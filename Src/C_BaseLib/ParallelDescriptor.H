#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

//
// $Id: ParallelDescriptor.H,v 1.65 2001-07-23 04:16:52 car Exp $
//

#include <BLassert.H>
#include <REAL.H>
#include <Box.H>
#include <ccse-mpi.H>

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
std::ostream& operator<< (std::ostream& os, const CommData& cd);

std::ostream& operator<< (std::ostream& os, const Array<CommData>& cd);

//
//@Man:
//@Memo: Parallel functions.
/*@Doc:

  This class contains functions used for implementing parallelism.
*/

namespace ParallelDescriptor
{
    class Message
    {
    public:
	Message();
	Message(MPI_Request req_, MPI_Datatype type_);
	Message(MPI_Status stat_, MPI_Datatype type_);
	void wait();
	bool test();
	size_t count() const;
	int tag() const;
	int pid() const;
	MPI_Datatype type() const;
    private:
	bool m_finished;
	MPI_Datatype m_type;
	MPI_Request m_req;
	mutable MPI_Status m_stat;
    };

    /*@ManDoc: Perform any needed parallel initialization.  This MUST be the
               first routine in this class called from within a program.
    */
    void StartParallel (int*    argc = 0,
			char*** argv = 0);

    /*@ManDoc: Perform any needed parallel finalization.  This MUST be the
               last routine in this class called from within a program.
    */
    void EndParallel ();
    //
    //@ManDoc: Returns processor number of calling program.
    //
    int MyProc ();
    //
    //@ManDoc: Returns number of CPUs involved in the computation.
    //
    int NProcs ();
    //
    //@ManDoc: Returns number of CPUs to use in CFD portion of computation.
    //
    int NProcsCFD ();
    //
    //@ManDoc: 
    //
    void Barrier ();
    //
    //@ManDoc: Is this CPU the I/O Processor?
    //
    bool IOProcessor ();
    //
    //@ManDoc: The CPU number of the I/O Processor.
    //
    int IOProcessorNumber ();
    //
    //@ManDoc: Issue architecture specific Abort.
    //
    void Abort ();
    //
    //@ManDoc: Abort with specified error code.
    //
    void Abort (int errorcode);
    //
    //@ManDoc: ErrorString return string associated with error internal error condition
    //
    const char* ErrorString(int errcode);
    //
    //@ManDoc: Returns wall-clock seconds since start of execution.
    //
    double second ();
    //
    //@ManDoc: And-wise boolean reduction.
    //
    void ReduceBoolAnd (bool& rvar);
    //
    //@ManDoc: And-wise boolean reduction to specified cpu.
    //
    void ReduceBoolAnd (bool& rvar, int cpu);
    //
    //@ManDoc: Or-wise boolean reduction.
    //
    void ReduceBoolOr  (bool& rvar);
    //
    //@ManDoc: Or-wise boolean reduction to specified cpu.
    //
    void ReduceBoolOr  (bool& rvar, int cpu);
    //
    //@ManDoc: Real sum reduction.
    //
    void ReduceRealSum (Real& rvar);
    //
    //@ManDoc: Real sum reduction to specified cpu.
    //
    void ReduceRealSum (Real& rvar, int cpu);
    //
    //@ManDoc: Real max reduction.
    //
    void ReduceRealMax (Real& rvar);
    //
    //@ManDoc: Real max reduction to specified cpu.
    //
    void ReduceRealMax (Real& rvar, int cpu);
    //
    //@ManDoc: Real min reduction.
    //
    void ReduceRealMin (Real& rvar);
    //
    //@ManDoc: Real min reduction to specified cpu.
    //
    void ReduceRealMin (Real& rvar, int cpu);
    //
    //@ManDoc: Integer sum reduction.
    //
    void ReduceIntSum (int& rvar);
    //
    //@ManDoc: Integer sum reduction to specified cpu.
    //
    void ReduceIntSum (int& rvar, int cpu);
    //
    //@ManDoc: Integer max reduction.
    //
    void ReduceIntMax (int& rvar);
    //
    //@ManDoc: Integer max reduction to specified cpu.
    //
    void ReduceIntMax (int& rvar, int cpu);
    //
    //@ManDoc: Integer min reduction.
    //
    void ReduceIntMin (int& rvar);
    //
    //@ManDoc: Integer min reduction to specified cpu.
    //
    void ReduceIntMin (int& rvar, int cpu);
    //
    //@ManDoc: Long sum reduction.
    //
    void ReduceLongSum (long& rvar);
    //
    //@ManDoc: Long sum reduction to specified cpu.
    //
    void ReduceLongSum (long& rvar, int cpu);
    //
    //@ManDoc: Long max reduction.
    //
    void ReduceLongMax (long& rvar);
    //
    //@ManDoc: Long max reduction to specified cpu.
    //
    void ReduceLongMax (long& rvar, int cpu);
    //
    //@ManDoc: Long min reduction.
    //
    void ReduceLongMin (long& rvar);
    //
    //@ManDoc: Long min reduction to specified cpu.
    //
    void ReduceLongMin (long& rvar, int cpu);
    //
    //@ManDoc: Long and-wise reduction.
    //
    void ReduceLongAnd (long& rvar);
    //
    //@ManDoc: Long and-wise reduction to specified cpu.
    //
    void ReduceLongAnd (long& rvar, int cpu);
    //
    //@ManDoc: Parallel gather.
    //
    void Gather (Real* sendbuf,
                        int   sendcount,
                        Real* recvbuf,
                        int   root);
    //
    //@ManDoc: Parallel broadcast of nbyte bytes from fromproc to all procs.
    //
    void Broadcast (int   fromproc,
                           void* src,
                           void* dest,
                           int   nbytes);
    //
    //@ManDoc: BoxLib's Parallel Communicator, probably MPI_COMM_WORLD
    //
    MPI_Comm Communicator ();
    //
    //@ManDoc: Returns sequential message sequence numbers in range 1000-9000.
    //
    int SeqNum ();

    template <class T> Message Asend(const T*, size_t n, int pid, int tag);
    template <class T> Message Asend(const std::vector<T>& buf, int pid, int tag);

    template <class T> Message Arecv(T*, size_t n, int pid, int tag);
    template <class T> Message Arecv(std::vector<T>& buf, int pid, int tag);

    template <class T> Message Send(const T* buf, size_t n, int dst_pid, int tag);
    template <class T> Message Send(const std::vector<T>& buf, int dst_pid, int tag);

    template <class T> Message Recv(T*, size_t n, int pid, int tag);
    template <class T> Message Recv(std::vector<T>& t, int pid, int tag);

    template <class T> void Bcast(T*, size_t n, int root = 0);
    template <class T> T Bcast(const T&, int root = 0);

    template <class Op, class T> T Reduce(const T& t);
    template <class Op, class T> std::vector<T> Reduce(const std::vector<T>& t);

    template <class T, class T1> void Scatter(T*, size_t n, const T1*, size_t n1, int root);
    template <class T> T Scatter(const std::vector<T>& t, int root);

    template <class T, class T1> void Gather(T*, size_t n, const T1*, size_t n1, int root);
    template <class T> std::vector<T> Gather(const T&, int root);
    void MPI_Error(const char* file, int line, const char* msg, int rc);
}

#define BL_MPI_REQUIRE(x)						\
do									\
{									\
  if ( int l_status_ = (x) )						\
    {									\
      ParallelDescriptor::MPI_Error(__FILE__,__LINE__,#x, l_status_);   \
    }									\
}									\
while ( false )

#if BL_USE_MPI
template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend(const T* buf, size_t n, int dst_pid, int tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf), n, Mpi_typemap<T>::type(), dst_pid, tag, Communicator(), &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend(const std::vector<T>& buf, int dst_pid, int tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(&buf[0]), buf.size(), Mpi_typemap<T>::type(), dst_pid, tag, Communicator(), &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send(const T* buf, size_t n, int dst_pid, int tag)
{
    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(buf), n, Mpi_typemap<T>::type(), dst_pid, tag, Communicator()) );
    return Message();
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send(const std::vector<T>& buf, int dst_pid, int tag)
{
    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(&buf[0]), buf.size(), Mpi_typemap<T>::type(), dst_pid, tag, Communicator()) );
    return Message();
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv(T* buf, size_t n, int src_pid, int tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf, n, Mpi_typemap<T>::type(), src_pid, tag, Communicator(), &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv(std::vector<T>& buf, int src_pid, int tag)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(&buf[0], buf.size(), Mpi_typemap<T>::type(), src_pid, tag, Communicator(), &req) );
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv(T* buf, size_t n, int src_pid, int tag)
{
    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(buf, n, Mpi_typemap<T>::type(), src_pid, tag, Communicator(), &stat) );
    return Message(stat, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv(std::vector<T>& buf, int src_pid, int tag)
{
    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(&buf[0], buf.size(), Mpi_typemap<T>::type(), src_pid, tag, Communicator(), &stat) );
    return Message(stat, Mpi_typemap<T>::type());
}

template <class Op, class T>
T
ParallelDescriptor::Reduce(const T& t)
{
    T recv;
    BL_MPI_REQUIRE( MPI_Allreduce(const_cast<T*>(&t), &recv, 1, Mpi_typemap<T>::type(), Op::op(), Communicator()) );
    return recv;
}

template <class Op, class T>
std::vector<T>
ParallelDescriptor::Reduce(const std::vector<T>& t)
{
    std::vector<T> recv;
    BL_MPI_REQUIRE( MPI_Allreduce(const_cast<T*>(&t[0]), &recv[0], t.size(), Mpi_typemap<T>::type(), Op::op(), Communicator()) );
    return recv;
}

template <class T>
void
ParallelDescriptor::Bcast(T* t, size_t n, int root)
{
    BL_MPI_REQUIRE( MPI_Bcast(t, n, Mpi_typemap<T>::type(), root, Communicator()) );
}

template <class T>
T
ParallelDescriptor::Bcast(const T& t, int root)
{
    T resl;
    if ( MyId() == root ) resl = t;
    BL_MPI_REQUIRE( MPI_Bcast(&resl, 1, Mpi_typemap<T>::type(), root, Communicator()) );
    return resl;
}

template <class T, class T1>
void
ParallelDescriptor::Gather(T* t, size_t n, const T1* t1, size_t n1, int root)
{
    BL_MPI_REQUIRE( MPI_Gather(t1, n1, Mpi_typemap<T1>::type(),
			    t, n, Mpi_typemap<T>::type(),
			    root, Communicator()) );
}

template <class T>
std::vector<T>
ParallelDescriptor::Gather(const T& t, int root)
{
    std::vector<T> resl;
    if ( root == MyProc() ) resl.resize(NProcs());
    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(&t), 1, Mpi_typemap<T>::type(),
			    &resl[0], 1, Mpi_typemap<T>::type(),
			    root, Communicator()) );
    return resl;
}

template <class T, class T1>
void
ParallelDescriptor::Scatter(T* t, size_t n, const T1* t1, size_t n1, int root)
{
    BL_MPI_REQUIRE( MPI_Scatter(const_cast<T1*>(t1), n1, Mpi_typemap<T1>::type(),
			     t, n, Mpi_typemap<T>::type(),
			     root, Communicator()) );
}

template <class T>
T
ParallelDescriptor::Scatter(const std::vector<T>& t, int root)
{
    T resl;
    BL_MPI_REQUIRE( MPI_Scatter(&t[0], 1, Mpi_typemap<T>::type(),
			     &resl, 1, Mpi_typemap<T>::type(),
			     root, Communicator()) );
    return resl;
}
#else
namespace ParallelDescriptor
{
template <class T>
Message
Asend(const T* buf, size_t n, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Asend(const std::vector<T>& buf, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Send(const T* buf, size_t n, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Send(const std::vector<T>& buf, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Arecv(T* buf, size_t n, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Arecv(std::vector<T>& buf, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Recv(T* buf, size_t n, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Recv(std::vector<T>& buf, int src_pid, int tag)
{
    return Message();
}

template <class Op, class T>
T
Reduce(const T& t)
{
    return t;
}

template <class Op, class T>
std::vector<T>
Reduce(const std::vector<T>& t)
{
    return t;
}

template <class T>
void
Bcast(T* t, size_t n, int root)
{
}

template <class T>
T
Bcast(const T& t, int root)
{
    return t;
}

template <class T, class T1>
void
ParallelDescriptor::Gather(T* t, size_t n, const T1* t1, size_t n1, int root)
{
}

template <class T>
std::vector<T>
Gather(const T& t, int root)
{
    std::vector<T> resl(1);
    resl[0] = t;
    return resl;
}

template <class T, class T1>
void
ParallelDescriptor::Scatter(T* t, size_t n, const T1* t1, size_t n1, int root)
{
}

template <class T>
T
Scatter(const std::vector<T>& t, int root)
{
    return t[0];
}

}
#endif

#endif /*BL_PARALLELDESCRIPTOR_H*/
