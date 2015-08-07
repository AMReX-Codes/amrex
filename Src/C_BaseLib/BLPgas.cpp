
#include <BLassert.H>
#include <BLPgas.H>

#ifdef BL_USE_CXX11
#include <unordered_map>
#else
#include <map>
#endif
namespace
{
#ifdef BL_USE_CXX11
  typedef std::unordered_multimap< upcxx::rank_t, BLPgas::SendInfo > pgas_send_info_map_t;
#else
  typedef std::multimap< upcxx::rank_t, BLPgas::SendInfo > pgas_send_info_map_t;
#endif
  static pgas_send_info_map_t pgas_send_info_map;
}

/**
 * \brief send a message using PGAS one-sided communication
 *
 * pgas_send is a drop-in replacement of MPI_Isend for non-blocking
 * communication.  It uses an Active Receive (Sender Side Tag Matching)
 * protocol to perform point-to-point message passing as follows.  The
 * receiver first sends its recv request (recv buffer address, message
 * size and tag) to the sender.  The sender does a simple tag matching with
 * the recv request and then initiate an one-sided put operation to transfer
 * the message from the send buffer to the recv buffer.
 *
 * @param src_ptr the source (send) buffer address
 * @param dst_ptr the destination (recv) buffer address
 * @param nbytes number of bytes for the message
 * @param SeqNum the BoxLib internal sequence number of message
 */
void
BLPgas::Send(upcxx::global_ptr<void> src_ptr,
             upcxx::global_ptr<void> dst_ptr,
             size_t nbytes,
             int SeqNum)
{

  // We use dst_rank as the key for the hash table
  std::pair <pgas_send_info_map_t::iterator, pgas_send_info_map_t::iterator> ret;
  ret = pgas_send_info_map.equal_range(dst_ptr.where());

  // try to match the SeqNum
  for (pgas_send_info_map_t::iterator it = ret.first;
       it != ret.second;
       ++it) {
    SendInfo& send_info = it->second;
    if (SeqNum == send_info.SeqNum) {
      // found the matched message
      // Check if data size matches
      BL_ASSERT(nbytes == send_info.nbytes);

      // Fire off the non-blocking one-sided communication (put)
      // If the receive request comes first, then the src_ptr in the existing
      // send_info entry should be NULL; otherwise, if the send request comes
      // first, then the dst_ptr must be NULL.  If neither is true, then there
      // must be something wrong!
      if (send_info.src_ptr == NULL) {
        upcxx::async_copy(src_ptr, send_info.dst_ptr, nbytes);
      } else {
        BL_ASSERT(send_info.dst_ptr == NULL);
        upcxx::async_copy(send_info.src_ptr, dst_ptr, nbytes);
      }

      // Delete the send_info from the map
      pgas_send_info_map.erase(it);
      return;
    }
  }

  // Can't find the send_info entry in the hash table
  // Create a new send_info and store the receiver part of it
  SendInfo send_info {src_ptr, dst_ptr, nbytes, SeqNum};
  pgas_send_info_map.insert(std::pair<upcxx::rank_t, SendInfo>(dst_ptr.where(), send_info));
}
