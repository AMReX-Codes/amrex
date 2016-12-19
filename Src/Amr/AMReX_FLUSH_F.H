#ifndef _FLUSH_F_H_
#define _FLUSH_F_H_

#if defined(BL_AIX)
#define FORT_FLUSH flush_
#else
#define FORT_FLUSH flush
#endif

#endif

