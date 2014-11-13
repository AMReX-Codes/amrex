
#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

#include <Profiler.H>
#include <BLassert.H>
#include <REAL.H>
#include <Box.H>
#include <ccse-mpi.H>

#include <string>
#include <typeinfo>

namespace ParallelDescriptor
{
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
        // We encapsulate seven ints and a Box as an int[3*BL_SPACEDIM+7].
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

        void nComp (int ncomp) { m_data[4] = ncomp; }
        //
        // The src component.
        //
        int srcComp () const { return m_data[5]; }

        void srcComp (int scomp) { m_data[5] = scomp; }
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
    // Used as default argument to ParallelDescriptor::Barrier().
    //
    const std::string Unnamed("Unnamed");
}
//
// Yes you can output CommData.
//
std::ostream& operator<< (std::ostream& os, const ParallelDescriptor::CommData& cd);

std::ostream& operator<< (std::ostream& os, const Array<ParallelDescriptor::CommData>& cd);
//
// Functions used for implementing parallelism.
//

namespace ParallelDescriptor
{
    class Message
    {
    public:

	Message () :
            m_finished(true),
            m_type(MPI_DATATYPE_NULL),
            m_req(MPI_REQUEST_NULL) {}
	Message (MPI_Request req_, MPI_Datatype type_) :
            m_finished(false),
            m_type(type_),
            m_req(req_) {}
	Message (MPI_Status stat_, MPI_Datatype type_) :
            m_finished(true),
            m_type(type_),
            m_req(MPI_REQUEST_NULL), m_stat(stat_) {}
	void wait ();
	bool test ();
	size_t count () const;
	int tag () const;
	int pid () const;
	MPI_Datatype type () const { return m_type; }
	MPI_Request  req () const { return m_req; }

    private:

	bool               m_finished;
	MPI_Datatype       m_type;
	MPI_Request        m_req;
	mutable MPI_Status m_stat;
    };
    //
    // Perform any needed parallel initialization.  This MUST be the
    // first routine in this class called from within a program.
    //
    void StartParallel (int*    argc = 0,
			char*** argv = 0,
                        MPI_Comm mpi_comm = MPI_COMM_WORLD);
    //
    // Perform any needed parallel finalization.  This MUST be the
    // last routine in this class called from within a program.
    //
    void EndParallel ();
    //
    // Returns processor number of calling program.
    //
    extern int m_MyId_all;
    extern int m_MyId_comp;
    extern int m_MyId_perfmon;
    extern int m_MyId_all_perfmon;  // id of perfmon proc in all
    inline int
    MyProc ()  // use _comp here, almost always what the user wants
    {
        BL_ASSERT(m_MyId_comp != -1);
        return m_MyId_comp;
    }
    inline int
    MyProcAll ()
    {
        BL_ASSERT(m_MyId_all != -1);
        return m_MyId_all;
    }
    inline int
    MyProcComp ()
    {
        BL_ASSERT(m_MyId_comp != -1);
        return m_MyId_comp;
    }
    inline int
    MyProcPerfMon ()
    {
        BL_ASSERT(m_MyId_perfmon != -1);
        return m_MyId_perfmon;
    }
    inline int
    MyProcAllPerfMon ()
    {
        BL_ASSERT(m_MyId_all_perfmon != -1);
        return m_MyId_all_perfmon;
    }
    //
    // Returns number of CPUs involved in the computation.
    //
    extern int m_nProcs_all;
    extern int m_nProcs_comp;
    extern int m_nProcs_perfmon;
    extern int nPerfmonProcs;
    inline void
    SetNProcsPerfMon (const int npmp)
    {
        BL_ASSERT(npmp >= 0);
        nPerfmonProcs = npmp;
    }
    inline int
    NProcs ()  // use _comp here, almost always what the user wants
    {
        BL_ASSERT(m_nProcs_comp != -1);
        return m_nProcs_comp;
    }
    inline int
    NProcsAll ()
    {
        BL_ASSERT(m_nProcs_all != -1);
        return m_nProcs_all;
    }
    inline int
    NProcsComp ()
    {
        BL_ASSERT(m_nProcs_comp != -1);
        return m_nProcs_comp;
    }
    inline int
    NProcsPerfMon ()
    {
        BL_ASSERT(m_nProcs_perfmon != -1);
        return m_nProcs_perfmon;
    }
    //
    // The CPU number of the I/O Processor.
    //
    extern const int ioProcessor;
    inline int
    IOProcessorNumber ()
    {
        return ioProcessor;
    }
    //
    // Is this CPU the I/O Processor?
    //
    inline bool
    IOProcessor ()
    {
         return MyProc() == IOProcessorNumber();
    }
    //
    // BoxLib's Parallel Communicators
    //    m_comm_all     is for all procs, probably MPI_COMM_WORLD
    //    m_comm_comp    is for the nodes doing computations
    //    m_comm_perfmon is for the in-situ performance monitor
    //
    extern MPI_Comm m_comm_all;
    extern MPI_Comm m_comm_comp;
    extern MPI_Comm m_comm_perfmon;
    extern MPI_Comm m_comm_inter;
    inline MPI_Comm Communicator ()  // use _comp here, almost always what the user wants
    {
        return m_comm_comp;
    }
    inline MPI_Comm CommunicatorAll ()
    {
        return m_comm_all;
    }
    inline MPI_Comm CommunicatorComp ()
    {
        return m_comm_comp;
    }
    inline MPI_Comm CommunicatorPerfMon ()
    {
        return m_comm_perfmon;
    }
    inline MPI_Comm CommunicatorInter ()
    {
        return m_comm_inter;
    }

    void Barrier (const std::string& message = Unnamed);
    void Barrier (MPI_Comm comm, const std::string& message = Unnamed);

    void Test (MPI_Request& request, int& flag, MPI_Status& status);

    void Comm_dup (MPI_Comm comm, MPI_Comm& newcomm);
    //
    // Issue architecture specific Abort.
    //
    void Abort ();
    //
    // Abort with specified error code.
    //
    void Abort (int errorcode);
    //
    // ErrorString return string associated with error internal error condition
    //
    const char* ErrorString (int errcode);
    //
    // Returns wall-clock seconds since start of execution.
    //
    double second ();
    //
    // And-wise boolean reduction.
    //
    void ReduceBoolAnd (bool& rvar);
    //
    // And-wise boolean reduction to specified cpu.
    //
    void ReduceBoolAnd (bool& rvar, int cpu);
    //
    // Or-wise boolean reduction.
    //
    void ReduceBoolOr  (bool& rvar);
    //
    // Or-wise boolean reduction to specified cpu.
    //
    void ReduceBoolOr  (bool& rvar, int cpu);
    //
    // Real sum reduction.
    //
    void ReduceRealSum (Real& rvar);

    void ReduceRealSum (Real* rvar, int cnt);
    //
    // Real sum reduction to specified cpu.
    //
    void ReduceRealSum (Real& rvar, int cpu);

    void ReduceRealSum (Real* rvar, int cnt, int cpu);
    //
    // Real max reduction.
    //
    void ReduceRealMax (Real& rvar);

    void ReduceRealMax (Real* rvar, int cnt);
    //
    // Real max reduction to specified cpu.
    //
    void ReduceRealMax (Real& rvar, int cpu);

    void ReduceRealMax (Real* rvar, int cnt, int cpu);
    //
    // Real min reduction.
    //
    void ReduceRealMin (Real& rvar);

    void ReduceRealMin (Real* rvar, int cnt);
    //
    // Real min reduction to specified cpu.
    //
    void ReduceRealMin (Real& rvar, int cpu);

    void ReduceRealMin (Real* rvar, int cnt, int cpu);
    //
    // Integer sum reduction.
    //
    void ReduceIntSum (int& rvar);

    void ReduceIntSum (int* rvar, int cnt);
    //
    // Integer sum reduction to specified cpu.
    //
    void ReduceIntSum (int& rvar, int cpu);

    void ReduceIntSum (int* rvar, int cnt, int cpu);
    //
    // Integer max reduction.
    //
    void ReduceIntMax (int& rvar);

    void ReduceIntMax (int* rvar, int cnt);
    //
    // Integer max reduction to specified cpu.
    //
    void ReduceIntMax (int& rvar, int cpu);

    void ReduceIntMax (int* rvar, int cnt, int cpu);
    //
    // Integer min reduction.
    //
    void ReduceIntMin (int& rvar);

    void ReduceIntMin (int* rvar, int cnt);
    //
    // Integer min reduction to specified cpu.
    //
    void ReduceIntMin (int& rvar, int cpu);

    void ReduceIntMin (int* rvar, int cnt, int cpu);
    //
    // Long sum reduction.
    //
    void ReduceLongSum (long& rvar);

    void ReduceLongSum (long* rvar, int cnt);
    //
    // Long sum reduction to specified cpu.
    //
    void ReduceLongSum (long& rvar, int cpu);

    void ReduceLongSum (long* rvar, int cnt, int cpu);
    //
    // Long max reduction.
    //
    void ReduceLongMax (long& rvar);

    void ReduceLongMax (long* rvar, int cnt);
    //
    // Long max reduction to specified cpu.
    //
    void ReduceLongMax (long& rvar, int cpu);

    void ReduceLongMax (long* rvar, int cnt, int cpu);
    //
    // Long min reduction.
    //
    void ReduceLongMin (long& rvar);

    void ReduceLongMin (long* rvar, int cnt);
    //
    // Long min reduction to specified cpu.
    //
    void ReduceLongMin (long& rvar, int cpu);

    void ReduceLongMin (long* rvar, int cnt, int cpu);
    //
    // Long and-wise reduction.
    //
    void ReduceLongAnd (long& rvar);

    void ReduceLongAnd (long* rvar, int cnt);
    //
    // Long and-wise reduction to specified cpu.
    //
    void ReduceLongAnd (long& rvar, int cpu);

    void ReduceLongAnd (long* rvar, int cnt, int cpu);
    // 
    // long long sum reduction.
    //
    void ReduceLongLongSum (long long& rvar);

    void ReduceLongLongSum (long long* rvar, int cnt);
    //
    // Long long sum reduction to specified cpu.
    //
    void ReduceLongLongSum (long long& rvar, int cpu);

    void ReduceLongLongSum (long long* rvar, int cnt, int cpu);
    //
    // Parallel gather.
    //
    void Gather (Real* sendbuf,
                 int   sendcount,
                 Real* recvbuf,
                 int   root);
    //
    // Returns sequential message sequence numbers in range 1000-9000.
    //
    int SeqNum ();

    template <class T> Message Asend(const T*, size_t n, int pid, int tag);
    template <class T> Message Asend(const T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Asend(const std::vector<T>& buf, int pid, int tag);

    template <class T> Message Arecv(T*, size_t n, int pid, int tag);
    template <class T> Message Arecv(T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Arecv(std::vector<T>& buf, int pid, int tag);

    template <class T> Message Send(const T* buf, size_t n, int dst_pid, int tag);
    template <class T> Message Send(const std::vector<T>& buf, int dst_pid, int tag);

    template <class T> Message Recv(T*, size_t n, int pid, int tag);
    template <class T> Message Recv(std::vector<T>& t, int pid, int tag);

    template <class T> void Bcast(T*, size_t n, int root = 0);

    template <class Op, class T> T Reduce(const T& t);

    template <class T, class T1> void Scatter(T*, size_t n, const T1*, size_t n1, int root);

    template <class T, class T1> void Gather(const T*, size_t n, T1*, size_t n1, int root);
    template <class T> std::vector<T> Gather(const T&, int root);

    void Waitsome (Array<MPI_Request>&, int&, Array<int>&, Array<MPI_Status>&);

    void MPI_Error(const char* file, int line, const char* msg, int rc);

    void ReadAndBcastFile(const std::string &filename, Array<char> &charBuf,
                          const bool bExitOnError = true);
    void IProbe(int src_pid, int tag, int &mflag, MPI_Status &status);
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
ParallelDescriptor::Asend (const T* buf,
                           size_t   n,
                           int      dst_pid,
                           int      tag)
{
    BL_COMM_PROFILE(Profiler::AsendTsii, n * sizeof(T), dst_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf),
                              n,
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(Profiler::AsendTsii, Profiler::AfterCall(), dst_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const T* buf,
                           size_t   n,
                           int      dst_pid,
                           int      tag,
                           MPI_Comm comm)
{
    BL_COMM_PROFILE(Profiler::AsendTsiiM, n * sizeof(T), dst_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf),
                              n,
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              comm,
                              &req) );
    BL_COMM_PROFILE(Profiler::AsendTsiiM, Profiler::AfterCall(), dst_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const std::vector<T>& buf,
                           int                   dst_pid,
                           int                   tag)
{
    BL_COMM_PROFILE(Profiler::AsendvTii, buf.size() * sizeof(T), dst_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(&buf[0]),
                              buf.size(),
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(Profiler::AsendvTii, Profiler::AfterCall(), dst_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const T* buf,
                          size_t   n,
                          int      dst_pid,
                          int      tag)
{
    BL_COMM_PROFILE(Profiler::SendTsii, n * sizeof(T), dst_pid, tag);

    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(buf),
                             n,
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             Communicator()) );
    BL_COMM_PROFILE(Profiler::SendTsii, Profiler::AfterCall(), dst_pid, tag);
    return Message();
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const std::vector<T>& buf,
                          int                   dst_pid,
                          int                   tag)
{
    BL_COMM_PROFILE(Profiler::SendvTii, buf.size() * sizeof(T), dst_pid, tag);

    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(&buf[0]),
                             buf.size(),
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             Communicator()) );
    BL_COMM_PROFILE(Profiler::SendvTii, Profiler::AfterCall(), dst_pid, tag);
    return Message();
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (T*       buf,
                           size_t   n,
                           int      src_pid,
                           int      tag)
{
    BL_COMM_PROFILE(Profiler::ArecvTsii, n * sizeof(T), src_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf,
                              n,
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(Profiler::ArecvTsii, Profiler::AfterCall(), src_pid, tag);

    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (T*       buf,
                           size_t   n,
                           int      src_pid,
                           int      tag,
                           MPI_Comm comm)
{
    BL_COMM_PROFILE(Profiler::ArecvTsiiM, n * sizeof(T), src_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf,
                              n,
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              comm,
                              &req) );
    BL_COMM_PROFILE(Profiler::ArecvTsiiM, Profiler::AfterCall(), src_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (std::vector<T>& buf,
                           int             src_pid,
                           int             tag)
{
    BL_COMM_PROFILE(Profiler::ArecvvTii, buf.size() * sizeof(T), src_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(&buf[0],
                              buf.size(),
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(Profiler::ArecvvTii, Profiler::AfterCall(), src_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (T*     buf,
                          size_t n,
                          int    src_pid,
                          int    tag)
{
    BL_COMM_PROFILE(Profiler::RecvTsii, Profiler::BeforeCall(), src_pid, tag);

    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(buf,
                             n,
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             Communicator(),
                             &stat) );
    BL_COMM_PROFILE(Profiler::RecvTsii, n * sizeof(T), stat.MPI_SOURCE, stat.MPI_TAG);
    return Message(stat, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (std::vector<T>& buf,
                          int             src_pid,
                          int             tag)
{
    BL_COMM_PROFILE(Profiler::RecvvTii, Profiler::BeforeCall(), src_pid, tag);

    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(&buf[0],
                             buf.size(),
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             Communicator(),
                             &stat) );
    BL_COMM_PROFILE(Profiler::RecvvTii, buf.size() * sizeof(T), stat.MPI_SOURCE, stat.MPI_TAG);
    return Message(stat, Mpi_typemap<T>::type());
}

template <class Op, class T>
T
ParallelDescriptor::Reduce (const T& t)
{
    BL_COMM_PROFILE_ALLREDUCE(Profiler::AllReduceT, Profiler::BeforeCall(), true);

    T recv;
    BL_MPI_REQUIRE( MPI_Allreduce(const_cast<T*>(&t),
                                  &recv,
                                  1,
                                  Mpi_typemap<T>::type(),
                                  Op::op(),
                                  Communicator()) );
    BL_COMM_PROFILE_ALLREDUCE(Profiler::AllReduceT, sizeof(T), false);
    return recv;
}

template <class T>
void
ParallelDescriptor::Bcast (T*     t,
                           size_t n,
                           int    root)
{
    BL_COMM_PROFILE(Profiler::BCastTsi, Profiler::BeforeCall(), root, Profiler::NoTag());

    BL_MPI_REQUIRE( MPI_Bcast(t,
                              n,
                              Mpi_typemap<T>::type(),
                              root,
                              Communicator()) );
    BL_COMM_PROFILE(Profiler::BCastTsi, n * sizeof(T), root, Profiler::NoTag());
}

template <class T, class T1>
void
ParallelDescriptor::Gather (const T* t,
                            size_t   n,
                            T1*      t1,
                            size_t   n1,
                            int      root)
{
    BL_COMM_PROFILE(Profiler::GatherTsT1Si, Profiler::BeforeCall(), root, Profiler::NoTag());

    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(t),
                               n,
                               Mpi_typemap<T>::type(),
                               t1,
                               n1,
                               Mpi_typemap<T1>::type(),
                               root,
                               Communicator()) );
    BL_COMM_PROFILE(Profiler::GatherTsT1Si,  n * sizeof(T), root, Profiler::NoTag());
}

template <class T>
std::vector<T>
ParallelDescriptor::Gather (const T& t, int root)
{
    BL_COMM_PROFILE(Profiler::GatherTi, Profiler::BeforeCall(), root, Profiler::NoTag());

    std::vector<T> resl;
    if ( root == MyProc() ) resl.resize(NProcs());
    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(&t),
                               1,
                               Mpi_typemap<T>::type(),
                               &resl[0],
                               1,
                               Mpi_typemap<T>::type(),
                               root,
                               Communicator()) );
    BL_COMM_PROFILE(Profiler::GatherTi, sizeof(T), root, Profiler::NoTag());
    return resl;
}

template <class T, class T1>
void
ParallelDescriptor::Scatter (T*        t,
                             size_t    n,
                             const T1* t1,
                             size_t    n1,
                             int       root)
{
    BL_COMM_PROFILE(Profiler::ScatterTsT1si,  Profiler::BeforeCall(), root, Profiler::NoTag());

    BL_MPI_REQUIRE( MPI_Scatter(const_cast<T1*>(t1),
                                n1,
                                Mpi_typemap<T1>::type(),
                                t,
                                n,
                                Mpi_typemap<T>::type(),
                                root,
                                Communicator()) );
    BL_COMM_PROFILE(Profiler::ScatterTsT1si, n * sizeof(T), root, Profiler::NoTag());
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
Asend(const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm)
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
Arecv(T* buf, size_t n, int src_pid, int tag, MPI_Comm comm)
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

template <class T>
void
Bcast(T* t, size_t n, int root)
{}

template <class T, class T1>
void
Gather(const T* t, size_t n, T1* t1, size_t n1, int root)
{}

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
Scatter(T* t, size_t n, const T1* t1, size_t n1, int root)
{}

}
#endif

#endif /*BL_PARALLELDESCRIPTOR_H*/
