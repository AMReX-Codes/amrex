
#ifndef BL_PGAS_H
#define BL_PGAS_H

#include <BLassert.H>
#include <upcxx.h>

namespace BLPgas
{
  struct SendInfo {
    upcxx::global_ptr<void> src_ptr;
    upcxx::global_ptr<void> dst_ptr;
    size_t nbytes;
    int SeqNum;
    upcxx::event *signal_event;
    upcxx::event *done_event;
    int *send_counter;

    SendInfo(upcxx::global_ptr<void> src,
             upcxx::global_ptr<void> dst,
             size_t nb,
             int SN,
             upcxx::event *signal_e,
             upcxx::event *done_e,
             int *send_c
             )
	    : src_ptr (src),
	      dst_ptr (dst),
	      nbytes  (nb),
	      SeqNum  (SN),
        signal_event(signal_e),
        done_event(done_e),
        send_counter(send_c)
    {}

    inline void check()
    {
      BL_ASSERT(src_ptr != NULL);
      BL_ASSERT(dst_ptr != NULL);
      BL_ASSERT(nbytes != 0);
      BL_ASSERT(signal_event != NULL);
      // done_event can be NULL
      BL_ASSERT(send_counter != NULL);
    }
  };

  void Send(upcxx::global_ptr<void> src_ptr,
            upcxx::rank_t dst_rank,
            size_t nbytes,
            int SeqNum,
            upcxx::event *done_event,
            int *send_counter);

  void Request(upcxx::rank_t src_rank,
	       upcxx::global_ptr<void> dst_ptr,
               size_t nbytes,
               int SeqNum,
               upcxx::event *signal_event);

  void Sendrecv(upcxx::global_ptr<void> src_ptr,
                upcxx::global_ptr<void> dst_ptr,
                size_t nbytes,
                int SeqNum,
                upcxx::event *signal_event,
                upcxx::event *done_event,
                int *send_counter);

  extern upcxx::event cp_send_event;
  extern upcxx::event cp_recv_event;
  extern int          cp_send_counter;

  extern upcxx::event fb_send_event;
  extern upcxx::event fb_recv_event;
  extern int          fb_send_counter;

  extern upcxx::event CollectData_recv_event;

  void* alloc (std::size_t _sz);
  void free (void* pt);  
}

#endif
