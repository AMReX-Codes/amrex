//BL_COPYRIGHT_NOTICE

#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

//
// $Id: ParallelDescriptor.H,v 1.37 1998-08-21 22:20:20 car Exp $
//

#include <Assert.H>
#include <Boolean.H>
#include <REAL.H>
#include <Box.H>

//
// Used for collecting information used in communicating FABs.
//
struct FabComTag
{
    int fromProc;
    int toProc;
    int fabIndex;
    int fineIndex;
    int srcComp;
    int destComp;
    int nComp;
    int face;
    int fabArrayId;
    int fillBoxId;
    int procThatNeedsData;
    int procThatHasData;
    Box box;
    //
    // A default constructor -- to quiet 3rd.
    //
    FabComTag ()
    {
        fromProc          = 0;
        toProc            = 0;
        fabIndex          = 0;
        fineIndex         = 0;
        srcComp           = 0;
        destComp          = 0;
        nComp             = 0;
        face              = 0;
        fabArrayId        = 0;
        fillBoxId         = 0;
        procThatNeedsData = 0;
        procThatHasData   = 0;
    }
};

//
// Data structure used by a few routines when MPI is enabled.
//
// Used to communicate up to seven integers and a box.
//
// We'll store all the info in a single array of integers.
//
struct CommData
{
    //
    // We encapsulate seven `int's and a `Box' as an `int[3*BL_SPACEDIM+7]'.
    //
    enum { DIM = 3*BL_SPACEDIM+7 };

    int m_data[DIM];

    CommData ();
    CommData (int        face,
              int        fabindex,
              int        fromproc,
              int        id,
              int        ncomp,
              int        srccomp,
              int        fabarrayid,
              const Box& box);
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
    // The processor sending this data.
    //
    int fromproc () const { return m_data[2]; }
    //
    // The ID of this message.
    //
    // Meant to be used as the MPI tag in a send/receive of additional
    // data associated with this data.
    //
    int id () const { return m_data[3]; }
    //
    // The number of components.
    //
    int nComp () const { return m_data[4]; }
    //
    // The src component.
    //
    int srcComp () const { return m_data[5]; }
    //
    // The ID of the fab array.
    //
    int fabarrayid () const { return m_data[6]; }
    //
    // The contained box.
    //
    Box box () const
    {
        return Box(IntVect(&m_data[7]),
                   IntVect(&m_data[7+BL_SPACEDIM]),
                   IntVect(&m_data[7+2*BL_SPACEDIM]));
    }
};

//
//@Man:
//@Memo: Parallel functions.
/*@Doc:

  This class contains functions used for implementing parallelism.
*/

class ParallelDescriptor
{
  public:

    /*@ManDoc: Perform any needed parallel initialization.  This MUST be the
               first routine in this class called from within a program.
    */
    static void StartParallel (int     nprocs,
                               int*    argc = 0,
                               char*** argv = 0);

    /*@ManDoc: Perform any needed parallel finalization.  This MUST be the
               last routine in this class called from within a program.
    */
    static void EndParallel ();
    //
    //@ManDoc: Returns processor number of calling program.
    //
    static int MyProc ();
    //
    //@ManDoc: Returns number of CPUs involved in the computation.
    //
    static int NProcs ();
    //
    //@ManDoc: 
    //
    static void Barrier ();
    //
    //@ManDoc: Is this CPU the I/O Processor?
    //
    static bool IOProcessor ();
    //
    //@ManDoc: The CPU number of the I/O Processor.
    //
    static int IOProcessorNumber ();
    //
    //@ManDoc: Issue architecture specific Abort.
    //
    static void Abort ();
    //
    //@ManDoc: Abort with specified error code.
    //
    static void Abort (int errorcode);
    //
    //@ManDoc: ErrorString return string associated with error internal error condition
    //
    static const char* ErrorString(int errcode);
    //
    //@ManDoc: Returns wall-clock seconds since start of execution.
    //
    static double second ();
    //
    //@ManDoc: And-wise boolean reduction.
    //
    static void ReduceBoolAnd (bool& rvar);
    //
    //@ManDoc: Or-wise boolean reduction.
    //
    static void ReduceBoolOr  (bool& rvar);
    //
    //@ManDoc: Real sum reduction.
    //
    static void ReduceRealSum (Real& rvar);
    //
    //@ManDoc: Real max reduction.
    //
    static void ReduceRealMax (Real& rvar);
    //
    //@ManDoc: Real min reduction.
    //
    static void ReduceRealMin (Real& rvar);
    //
    //@ManDoc: Integer sum reduction.
    //
    static void ReduceIntSum (int& rvar);
    //
    //@ManDoc: Integer max reduction.
    //
    static void ReduceIntMax (int& rvar);
    //
    //@ManDoc: Integer min reduction.
    //
    static void ReduceIntMin (int& rvar);
    //
    //@ManDoc: Long sum reduction.
    //
    static void ReduceLongSum (long& rvar);
    //
    //@ManDoc: Long max reduction.
    //
    static void ReduceLongMax (long& rvar);
    //
    //@ManDoc: Long min reduction.
    //
    static void ReduceLongMin (long& rvar);
    //
    //@ManDoc: Long and-wise reduction.
    //
    static void ReduceLongAnd (long& rvar);
    //
    //@ManDoc: Parallel gather.
    //
    static void Gather (Real* sendbuf,
                        int   sendcount,
                        Real* recvbuf,
                        int   root);

    static void Broadcast (int   fromproc,
                           void* src,
                           void* dest,
                           int   nbytes);
  private:
    //
    // My processor ID.
    //
    static int m_MyId;
    //
    // The number of processors.
    //
    static int m_nProcs;

    enum { ioProcessor = 0 };
};

//
// Inlines.
//

inline
int
ParallelDescriptor::MyProc ()
{
    assert(m_MyId != -1);
    return m_MyId;
}

inline
int
ParallelDescriptor::NProcs ()
{
    assert(m_nProcs != -1);
    return m_nProcs;
}

inline
bool
ParallelDescriptor::IOProcessor ()
{
    assert(m_MyId != -1);
    return m_MyId == ioProcessor;
}

inline
int
ParallelDescriptor::IOProcessorNumber ()
{
    return ioProcessor;
}

#endif /*BL_PARALLELDESCRIPTOR_H*/
