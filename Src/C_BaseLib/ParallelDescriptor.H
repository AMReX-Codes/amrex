//BL_COPYRIGHT_NOTICE

#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

//
// $Id: ParallelDescriptor.H,v 1.2 1997-09-15 19:40:06 lijewski Exp $
//

#include <Assert.H>
#include <Boolean.H>
#include <REAL.H>
#include <Box.H>
#include <BaseFab.H>
#include <FArrayBox.H>

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

#include <Utility.H>

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
};

#ifndef BL_USE_BSP

class ParallelDescriptor
{
public:

    static void Abort (const char* str)
    {
	cerr << "Abort:  " << str << endl;
	abort();
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
    static double second ()  { return Utility::second(); }

    // reduction operations
    template<class T> static void ReduceMin (T &rvar) {}

    // bool
    static void ReduceBoolAnd (bool& rvar) {}
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

    static void SetMessageHeaderSize (int& messageHeaderSize) {} 

    static bool GetMessageHeader (int&  dataSize,
                                  void* messageHeader)
    {
        return false;  // no messages waiting
    } 

    static void SendData (int              toproc,
                          const FabComTag& messageHeader,
                          const void*      data,
                          int              datasizeinbytes) {}

    static void ReceiveData (void* data,
                             int   datasizeinbytes) {}

    template<class T> static void Broadcast (int toproc,
                                             T*  src,
                                             T*  dest) {}
};

#else

class ParallelDescriptor
{
  public:

    static int MyProc () { return bsp_pid();    }
    static int NProcs () { return bsp_nprocs(); }
    static void Synchronize () { bsp_sync(); }
    static void Synchronize (const char *msg)
    {
	  cout << "----- " << bsp_pid() << " :  about to sync:  " << msg
               << endl;
          bsp_sync();
    }
    static bool IOProcessor () { return (bsp_pid() == ioProcessor); }
    static int  IOProcessorNumber () { return (ioProcessor); }
    static void Abort (char *msg, ...)  { bsp_abort(msg); }
    static double second ()  { return bsp_time(); }

    // reduction operations
    template<class T> static void ReduceMin (T& rvar);

    // bool
    static void ReduceBoolAnd (bool& rvar);
    // Real
    static void ReduceRealSum (Real& rvar);
    static void ReduceRealMax (Real& rvar);
    static void ReduceRealMin (Real& rvar);
    // int
    static void ReduceIntSum (int& rvar);
    static void ReduceIntMax (int& rvar);
    static void ReduceIntMin (int& rvar);
    // long
    static void ReduceLongSum (long& rvar);
    static void ReduceLongMax (long& rvar);
    static void ReduceLongMin (long& rvar);

    // data transfer functions
    static void ShareVar (const void* var,
                          int         bytes)
    {
        bsp_pushregister(var, bytes);
    }
    static void UnshareVar (const void* var)
    {
        bsp_popregister(var);
    }
    static void WriteData (int         procnum,
                           const void* src,
                           void*       dest,
                           int         offset,
                           int         bytes)
    {
        bsp_put(procnum, src, dest, offset, bytes);
    }
    static void ReadData (int         procnum,
                          const void* src,
                          int         offset,
                          void*       dest,
                          int         bytes)
    {
        bsp_get(procnum, src, offset, dest, bytes);
    }
    static void SetMessageHeaderSize (int& messageHeaderSize)
    {
        bsp_set_tag_size(&messageHeaderSize);
    } 
    static bool GetMessageHeader (int& dataSize,
                                  void* messageHeader)
    {
        bsp_get_tag(&dataSize, messageHeader);

        return dataSize == -1
            ? false  // By bsp definition
            : true;  // A message is waiting
    } 
    static void SendData (int              toproc,
                          const FabComTag& messageHeader,
                          const void*      data,
                          int              datasizeinbytes)
    {
        bsp_send(toproc, &messageHeader, data, datasizeinbytes);
    }
    static void ReceiveData (void* data,
                             int   datasizeinbytes)
    {
        bsp_move(data, datasizeinbytes);
    }
    template<class T> static void Broadcast (int toproc,
                                             T*  src,
                                             T*  dest)
    {
        bsp_bcast(toproc, (void *) src, (void *) dest, sizeof(T));
    }

    //static void SendFab(int procnum, const FArrayBox &fab);

  private:

    enum { ioProcessor = 0 };
};

#endif /*BL_USE_BSP*/

#endif /*BL_PARALLELDESCRIPTOR_H*/
