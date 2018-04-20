
#include <AMReX_BLassert.H>
#include <AMReX_BLPgas.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Print.H>

//#include <unordered_map>
#include <map>

namespace amrex {

namespace
{
  // typedef std::unordered_multimap< upcxx::rank_t, BLPgas::SendInfo > pgas_send_info_map_t;
  typedef std::multimap< upcxx::rank_t, BLPgas::SendInfo > pgas_send_info_map_t;
  static pgas_send_info_map_t pgas_send_info_map;
}

namespace BLPgas {
  upcxx::event cp_send_event;
  upcxx::event cp_recv_event;
  int          cp_send_counter;

  upcxx::event fb_send_event;
  upcxx::event fb_recv_event;
  int          fb_send_counter;

  upcxx::event CollectData_recv_event;
}

/**
 * \brief send a message using PGAS one-sided communication
 *
 * BLPgas::Send is a drop-in replacement of MPI_Isend for non-blocking
 * communication. It uses an Active Receive (Sender Side Tag Matching)
 * protocol to perform point-to-point message passing as follows. The
 * receiver sends its recv request (recv buffer address, message size
 * and tag) to the sender using BLPgas::Request. The sender does a
 * simple tag matching with the recv request and then initiate an
 * one-sided put operation to transfer the message from the send
 * buffer to the recv buffer.
 *
 * @param src the source (send) buffer address
 * @param dst_rank the destination rank
 * @param nbytes number of bytes for the message
 * @param SeqNum the AMReX internal sequence number of message
 * @param done_event for notifying the sender when the data transfer is done
 * @param send_counter increment the counter only if the message is sent out
 */
void
BLPgas::Send(upcxx::global_ptr<void> src,
             upcxx::rank_t dst_rank,
             size_t nbytes,
             int SeqNum,
             upcxx::event *done_event,
             int *send_counter)
{
  BLPgas::Sendrecv(src, upcxx::global_ptr<void>(NULL, dst_rank),
                   nbytes, SeqNum, (upcxx::event*)NULL, done_event,
                   send_counter);
}

/**
 * \brief receive a message using PGAS one-sided communication
 *
 * BLPgas::Request is a drop-in replacement of MPI_Irecv for
 * non-blocking communication. It uses an Active Receive (Sender Side
 * Tag Matching) protocol to perform point-to-point message passing as
 * follows. The receiver sends its recv request (recv buffer address,
 * message size and tag) to the sender. The sender calls BLPgas::Send
 * to do a simple tag matching with the recv request and then initiate
 * an one-sided put operation to transfer the message from the send
 * buffer to the recv buffer.
 *
 * @param src_rank the source rank
 * @param dst the destination (recv) buffer address
 * @param nbytes number of bytes for the message
 * @param SeqNum the AMReX internal sequence number of message
 * @param signal_event for notifying the receiver when the data transfer is done
 */
void
BLPgas::Request(upcxx::rank_t src_rank,
                upcxx::global_ptr<void> dst,
                size_t nbytes,
                int SeqNum,
                upcxx::event *signal_event)
{
  upcxx::async(src_rank, NULL)(BLPgas::Sendrecv,
			 upcxx::global_ptr<void>(NULL, src_rank), dst,
			 nbytes, SeqNum, signal_event, (upcxx::event*)NULL,
			 (int *)NULL);
}

//! A common implementation for the above send and request functions
void
BLPgas::Sendrecv(upcxx::global_ptr<void> src,
                 upcxx::global_ptr<void> dst,
                 size_t nbytes,
                 int SeqNum,
                 upcxx::event *signal_event,
                 upcxx::event *done_event,
                 int *send_counter)
{
  // We use dst_rank as the key for the hash table
  std::pair <pgas_send_info_map_t::iterator, pgas_send_info_map_t::iterator> ret;
  ret = pgas_send_info_map.equal_range(dst.where());

  // try to match the SeqNum
  for (pgas_send_info_map_t::iterator it = ret.first;
       it != ret.second;
       ++it) {
    SendInfo& send_info = it->second;
    if (SeqNum == send_info.SeqNum) {
      // found the matched message
      // Check if data size matches
      assert(nbytes == send_info.nbytes);

      // Fire off the non-blocking one-sided communication (put)
      // If the receive request comes first, then the src_ptr in the existing
      // send_info entry should be NULL; otherwise, if the send request comes
      // first, then the dst_ptr must be NULL.  If neither is true, then there
      // must be something wrong!
      if (send_info.src_ptr == NULL) {
        // pgas_send request from Recv came earlier
        send_info.src_ptr = src;
        send_info.done_event = done_event;
        send_info.send_counter = send_counter;
      } else {
        // pgas_send request from Send came earlier
        assert(send_info.dst_ptr == NULL);
        send_info.dst_ptr = dst;
        send_info.signal_event = signal_event;
      }

      send_info.check();

      async_copy_and_signal(send_info.src_ptr,
                            send_info.dst_ptr,
                            send_info.nbytes,
                            send_info.signal_event,
                            send_info.done_event,
                            NULL);

      (*send_info.send_counter)++;
      // Delete the send_info from the map
      pgas_send_info_map.erase(it);
      return;
    }
  }

  // Can't find the send_info entry in the hash table
  // Create a new send_info and store the receiver part of it
  SendInfo send_info {src, dst, nbytes, SeqNum, signal_event, done_event, send_counter};
  pgas_send_info_map.insert(std::pair<upcxx::rank_t, SendInfo>(dst.where(), send_info));
}

//! Allocate _sz bytes in the global memory address space
void*
BLPgas::alloc (std::size_t _sz)
{
    if (_sz <= 0) {
	return nullptr;
    } else {
	auto p = upcxx::allocate(_sz);
	if (p == nullptr) {
	    amrex::Print(Print::AllProcs)
		<< "===== Proc. " << ParallelDescriptor::MyProc() << " =====\n"
		<< "   Failed to allocate " << _sz << " bytes global address space memory!\n"
		<< "   Please try to increase the GASNET_MAX_SEGSIZE environment variable.\n"
		<< "   For example, export GASNET_MAX_SEGSIZE=512MB\n"
		<< "   Total Bytes Allocated in Fabs: " << amrex::TotalBytesAllocatedInFabs()
		<< "   Highest Watermark in Fabs: " << amrex::TotalBytesAllocatedInFabsHWM()
		<< "\n";
	    amrex::Abort("BLPgas: upcxx::allocate failed");
	}
	return p;
    }
}

void
BLPgas::free (void* pt)
{
    upcxx::deallocate(pt);
}

}
