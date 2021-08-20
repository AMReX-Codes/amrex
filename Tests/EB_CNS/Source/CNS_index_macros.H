#ifndef CNS_INDEX_H_
#define CNS_INDEX_H_

#if (AMREX_SPACEDIM == 2)

#define URHO  0
#define UMX   1
#define UMY   2
#define UEDEN 3
#define UEINT 4
#define UTEMP 5
#define NUM_STATE 6

#define NCONS 6
#define NEQNS 4

#define QRHO   0
#define QU     1
#define QV     2
#define QPRES  3
#define QCS    4
#define QEINT  5
#define QTEMP  6
#define NPRIM  7

#else

#define URHO  0
#define UMX   1
#define UMY   2
#define UMZ   3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define NUM_STATE 7

#define NCONS 7
#define NEQNS 5

#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QPRES  4
#define QCS    5
#define QEINT  6
#define QTEMP  7
#define NPRIM  8

#endif

#define CETA   0
#define CXI    1
#define CLAM   2
#define NCOEF  3

#endif
