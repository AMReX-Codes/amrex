//BL_COPYRIGHT_NOTICE

#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

//
// $Id: ParallelDescriptor.H,v 1.25 1998-04-16 22:36:16 lijewski Exp $
//

#include <Assert.H>
#include <Boolean.H>
#include <REAL.H>
#include <Box.H>

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

//
// Data structure used by a few routines when MPI is enabled.
//
// Used to communicate two integers and a box.
//
// We'll store all the info in a single array of integers.
//
struct CommData
{
    //
    // We encapsulate two `int's and a `Box' as an `int[3*BL_SPACEDIM+2]'.
    //
    enum { DIM = 3*BL_SPACEDIM+2 };

    int m_data[DIM];

    CommData ();
    CommData (int face, int fabindex, const Box& box);
    //
    // The number of integers.
    //
    int length () const { return DIM; }
    //
    // Pointer to the data.
    //
    int* dataPtr() { return &m_data[0]; }
    //
    // The face.
    //
    int face () const { return m_data[0]; }
    //
    // The fabindex.
    //
    int fabindex () const { return m_data[1]; }
    //
    // The contained box.
    //
    Box box () const
    {
        return Box(IntVect(&m_data[2]),
                   IntVect(&m_data[2+BL_SPACEDIM]),
                   IntVect(&m_data[2+2*BL_SPACEDIM]));
    }
};

class ParallelDescriptor
{
  public:
    static void StartParallel (int nprocs,
                               int* argc = 0,
                               char*** argv = 0);
    static void EndParallel ();

    static int    MyProc ();
    static int    NProcs ();
    static void   Barrier ();
    static void   Synchronize ();
    static bool   IOProcessor ();
    static int    IOProcessorNumber ();
    static void   Abort (const char* msg);
    static void   Abort (int errorcode);
    static double second ();
    //
    // Reduction operations.
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
    static void ShareVar (const void* var, int bytes);
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
    static bool GetMessageHeader (int&  dataSize,
                                  void* messageHeader);
    static bool MessageQueueEmpty ();

    static void SendData (int         toproc,
                          const void* messageHeader,
                          const void* data,
                          int         datasizeinbytes);

    static void ReceiveData (void* data,
                             int   datasizeinbytes);

    static void Broadcast (int   fromproc,
                           void* src,
                           void* dest,
                           int   nbytes);

    static void Gather (Real* sendbuf,
                        int   sendcount,
                        Real* recvbuf,
                        int   root);
  private:

    enum { ioProcessor = 0 };
};

#endif /*BL_PARALLELDESCRIPTOR_H*/
