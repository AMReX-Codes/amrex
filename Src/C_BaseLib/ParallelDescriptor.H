
#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

//
// $Id: ParallelDescriptor.H,v 1.51 2001-01-25 23:52:37 lijewski Exp $
//

#include <BLassert.H>
#include <REAL.H>
#include <Box.H>
#include <ccse-mpi.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

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

    CommData (const CommData& rhs);

    CommData& operator= (const CommData& rhs);
    //
    // Compare two CommData's.
    //
    bool operator== (const CommData& rhs) const;

    bool operator!= (const CommData& rhs) const { return !operator==(rhs); }
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
// Yes you can output CommData.
//
ostream& operator<< (ostream& os, const CommData& cd);

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
    static void StartParallel (int*    argc = 0,
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
    //@ManDoc: Returns number of CPUs to use in CFD portion of computation.
    //
    static int NProcsCFD ();
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
    //@ManDoc: And-wise boolean reduction to specified cpu.
    //
    static void ReduceBoolAnd (bool& rvar, int cpu);
    //
    //@ManDoc: Or-wise boolean reduction.
    //
    static void ReduceBoolOr  (bool& rvar);
    //
    //@ManDoc: Or-wise boolean reduction to specified cpu.
    //
    static void ReduceBoolOr  (bool& rvar, int cpu);
    //
    //@ManDoc: Real sum reduction.
    //
    static void ReduceRealSum (Real& rvar);
    //
    //@ManDoc: Real sum reduction to specified cpu.
    //
    static void ReduceRealSum (Real& rvar, int cpu);
    //
    //@ManDoc: Real max reduction.
    //
    static void ReduceRealMax (Real& rvar);
    //
    //@ManDoc: Real max reduction to specified cpu.
    //
    static void ReduceRealMax (Real& rvar, int cpu);
    //
    //@ManDoc: Real min reduction.
    //
    static void ReduceRealMin (Real& rvar);
    //
    //@ManDoc: Real min reduction to specified cpu.
    //
    static void ReduceRealMin (Real& rvar, int cpu);
    //
    //@ManDoc: Integer sum reduction.
    //
    static void ReduceIntSum (int& rvar);
    //
    //@ManDoc: Integer sum reduction to specified cpu.
    //
    static void ReduceIntSum (int& rvar, int cpu);
    //
    //@ManDoc: Integer max reduction.
    //
    static void ReduceIntMax (int& rvar);
    //
    //@ManDoc: Integer max reduction to specified cpu.
    //
    static void ReduceIntMax (int& rvar, int cpu);
    //
    //@ManDoc: Integer min reduction.
    //
    static void ReduceIntMin (int& rvar);
    //
    //@ManDoc: Integer min reduction to specified cpu.
    //
    static void ReduceIntMin (int& rvar, int cpu);
    //
    //@ManDoc: Long sum reduction.
    //
    static void ReduceLongSum (long& rvar);
    //
    //@ManDoc: Long sum reduction to specified cpu.
    //
    static void ReduceLongSum (long& rvar, int cpu);
    //
    //@ManDoc: Long max reduction.
    //
    static void ReduceLongMax (long& rvar);
    //
    //@ManDoc: Long max reduction to specified cpu.
    //
    static void ReduceLongMax (long& rvar, int cpu);
    //
    //@ManDoc: Long min reduction.
    //
    static void ReduceLongMin (long& rvar);
    //
    //@ManDoc: Long min reduction to specified cpu.
    //
    static void ReduceLongMin (long& rvar, int cpu);
    //
    //@ManDoc: Long and-wise reduction.
    //
    static void ReduceLongAnd (long& rvar);
    //
    //@ManDoc: Long and-wise reduction to specified cpu.
    //
    static void ReduceLongAnd (long& rvar, int cpu);
    //
    //@ManDoc: Parallel gather.
    //
    static void Gather (Real* sendbuf,
                        int   sendcount,
                        Real* recvbuf,
                        int   root);
    //
    //@ManDoc: Parallel broadcast of nbyte bytes from fromproc to all procs.
    //
    static void Broadcast (int   fromproc,
                           void* src,
                           void* dest,
                           int   nbytes);
    //
    //@ManDoc: BoxLib's Parallel Communicator, probably MPI_COMM_WORLD
    //
  static MPI_Comm Communicator ();

  private:
    //
    // Reduce helper functons.
    //
    static void DoAllReduceReal (Real& r, MPI_Op op);
    static void DoAllReduceLong (long& r, MPI_Op op);
    static void DoAllReduceInt (int& r, MPI_Op op);
    static void DoReduceReal (Real& r, MPI_Op op, int cpu);
    static void DoReduceLong (long& r, MPI_Op op, int cpu);
    static void DoReduceInt (int& r, MPI_Op op, int cpu);
    //
    // Sets number of CPUs to use in CFD portion of computation via ParmParse.
    //
    static void SetNProcsCFD ();
    //
    // My processor ID.
    //
    static int m_MyId;
    //
    // The number of processors.
    //
    static int m_nProcs;
    //
    // The number of processors in CFD part of computation.
    //
    static int m_nProcsCFD;
    //
    // BoxLib's Communicator
    //
    static MPI_Comm m_comm;

    enum { ioProcessor = 0 };
};

//
// Inlines.
//

inline
int
ParallelDescriptor::MyProc ()
{
    BL_ASSERT(m_MyId != -1);
    return m_MyId;
}

inline
int
ParallelDescriptor::NProcs ()
{
    BL_ASSERT(m_nProcs != -1);
    return m_nProcs;
}

inline
int
ParallelDescriptor::NProcsCFD ()
{
    if (m_nProcsCFD == -1)
        SetNProcsCFD();

    BL_ASSERT(m_nProcsCFD != -1);

    return m_nProcsCFD;
}

inline
bool
ParallelDescriptor::IOProcessor ()
{
    BL_ASSERT(m_MyId != -1);
    return m_MyId == ioProcessor;
}

inline
int
ParallelDescriptor::IOProcessorNumber ()
{
    return ioProcessor;
}

inline
MPI_Comm
ParallelDescriptor::Communicator ()
{
    return m_comm;
}

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_PARALLELDESCRIPTOR_H*/
