//BL_COPYRIGHT_NOTICE

#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

//
// $Id: ParallelDescriptor.H,v 1.17 1998-03-24 20:18:52 car Exp $
//

#include <Assert.H>
#include <Boolean.H>
#include <REAL.H>
#include <Box.H>

#ifdef BL_USE_BSP
#include "bsp.h"

#ifdef FIXBSPLIBLEVEL1HEADER
extern "C"
{
  extern void bsp_fold(void (*)(void *, void *, void *, int *),
                       void*,void*,int);
  extern void bsp_fold_cpp(void (*)(void *, void *, void *, int *),
                           void*,void* , int, int, char *);
}
#endif /*FIXBSPLIBLEVEL1HEADER*/

#include "bsp_level1.h"

#define StartParallel(nprocs)    bsp_begin(nprocs)
#define StartParallelAllProcs()  bsp_begin(bsp_nprocs())
#define EndParallel()            bsp_end()

#else

#define StartParallel(nprocs)
#define StartParallelAllProcs() 
#define EndParallel()

#endif /*BL_USE_BSP*/

//
// Some of these can be collapsed.
//
struct FabComTag
{
    int fromProc;
    int toProc;
    int fabIndex;
    int fineIndex;
    int crseIndex;
    int srcComp;
    int destComp;
    int nComp;
    int face;
    int dataType;
    int fabArrayId;
    int fillBoxId;
    int procThatNeedsData;
    int procThatHasData;
    Box box;
    Box ovlpBox;
    Box subBox;
    //
    // A default constructor -- to quiet 3rd.
    //
    FabComTag ()
    {
        fromProc          = 0;
        toProc            = 0;
        fabIndex          = 0;
        fineIndex         = 0;
        crseIndex         = 0;
        srcComp           = 0;
        destComp          = 0;
        nComp             = 0;
        face              = 0;
        dataType          = 0;
        fabArrayId        = 0;
        fillBoxId         = 0;
        procThatNeedsData = 0;
        procThatHasData   = 0;
    }
};

#if defined(BL_USE_BSP)

class ParallelDescriptor
{
  public:

    static int MyProc ();
    static int NProcs ();
    static void Synchronize ();
    static void Synchronize (const char* msg);
    static bool IOProcessor ();
    static int  IOProcessorNumber ();
    static void Abort (const char* msg);
    static double second ();
    //
    // Reduction operations.
    //
    template<class T> static void ReduceMin (T& rvar);
    //
    // bool
    //
    static void ReduceBoolAnd (bool& rvar);
    static void ReduceBoolOr  (bool& rvar);
    //
    // Real
    //
    static void ReduceRealSum (Real& rvar);
    static void ReduceRealMax (Real& rvar);
    static void ReduceRealMin (Real& rvar);
    //
    // int
    //
    static void ReduceIntSum (int& rvar);
    static void ReduceIntMax (int& rvar);
    static void ReduceIntMin (int& rvar);
    //
    // long
    //
    static void ReduceLongSum (long& rvar);
    static void ReduceLongMax (long& rvar);
    static void ReduceLongMin (long& rvar);
    static void ReduceLongAnd (long& rvar);
    //
    // Data transfer functions.
    //
    static void ShareVar (const void* var,
                          int         bytes);
    static void UnshareVar (const void* var);
    static void WriteData (int         procnum,
                           const void* src,
                           void*       dest,
                           int         offset,
                           int         bytes);
    static void ReadData (int         procnum,
                          const void* src,
                          int         offset,
                          void*       dest,
                          int         bytes);
    static void SetMessageHeaderSize (int messageHeaderSize);
    static bool GetMessageHeader (int& dataSize,
                                  void* messageHeader);
    static bool MessageQueueEmpty ();

    static void SendData (int         toproc,
                          const void* messageHeader,
                          const void* data,
                          int         datasizeinbytes);
    static void ReceiveData (void* data,
                             int   datasizeinbytes);
    static void Broadcast (int fromproc, void*  src, void*  dest, int nbytes);

  private:

    enum { ioProcessor = 0 };
};

#elif defined(BL_USE_MPI)

#else


class ParallelDescriptor
{
public:

    static void Abort (const char* str)
    {
        BoxLib::Abort(str);
    }
    static int MyProc () { return 0; }
    static int NProcs () { return 1; }
    static void Synchronize () {}
    static void Synchronize (const char* msg)
    {
        cout << "----- " << 0 << " :  about to sync:  " << msg << endl;
    }
    static bool IOProcessor () { return true; }
    static int  IOProcessorNumber () { return 0; }
    static double second ();

    // reduction operations
    template<class T> static void ReduceMin (T &rvar) {}

    // bool
    static void ReduceBoolAnd (bool& rvar) {}
    static void ReduceBoolOr  (bool& rvar) {}
    // Real
    static void ReduceRealSum (Real& rvar) {}
    static void ReduceRealMax (Real& rvar) {}
    static void ReduceRealMin (Real& rvar) {}
    // int
    static void ReduceIntSum (int& rvar) {}
    static void ReduceIntMax (int& rvar) {}
    static void ReduceIntMin (int& rvar) {}
    // long
    static void ReduceLongSum (long& rvar) {}
    static void ReduceLongMax (long& rvar) {}
    static void ReduceLongMin (long& rvar) {}
    static void ReduceLongAnd (long& rvar) {}

    // data transfer functions
    static void ShareVar (const void* var,
                          int         bytes) {}
    static void UnshareVar (const void* var) {}
    static void WriteData (int         procnum,
                           const void* src,
                           void*       dest,
                           int         offset,
                           int         bytes) {}

    static void ReadData (int         procnum,
                          const void* src,
                          int         offset,
                          void*       dest,
                          int         bytes) {}

    static void SetMessageHeaderSize (int messageHeaderSize) {} 

    static bool GetMessageHeader (int&  dataSize,
                                  void* messageHeader)
    {
        return false;  // no messages waiting
    } 
    static bool MessageQueueEmpty ()
    {
        return true;  // no messages waiting
    } 
    static void SendData (int         toproc,
                          const void* messageHeader,
                          const void* data,
                          int         datasizeinbytes) {}

    static void ReceiveData (void* data,
                             int   datasizeinbytes) {}

    static void Broadcast (int fromproc, void*  src, void*  dest, size_t nbytes) {}

};

#endif /*BL_USE_BSP*/

#endif /*BL_PARALLELDESCRIPTOR_H*/
