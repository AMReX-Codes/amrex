
#ifndef BL_PGAS_H
#define BL_PGAS_H

#include <upcxx.h>

namespace BLPgas
{
    struct SendInfo {
	upcxx::global_ptr<void> src_ptr;
	upcxx::global_ptr<void> dst_ptr;
	size_t nbytes;
	int SeqNum;
	
	SendInfo(upcxx::global_ptr<void> src,
		 upcxx::global_ptr<void> dst,
		 size_t nb,
		 int SN)
	    : src_ptr (src),
	      dst_ptr (dst),
	      nbytes  (nb),
	      SeqNum  (SN)
	    { ; }
    };

    void Send(upcxx::global_ptr<void> src_ptr,
	      upcxx::global_ptr<void> dst_ptr,
	      size_t nbytes,
	      int SeqNum);
}

#endif
