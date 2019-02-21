#ifndef AMREX_PARALLEL_REDUCE_H_
#define AMREX_PARALLEL_REDUCE_H_

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <type_traits>

namespace amrex {

namespace detail {

    enum ReduceOp : int {
        max = 0,
        min,
        sum,
        lor,
        land
    };

#ifdef BL_USE_MPI

    // NOTE: the order of these needs to match the order in the ReduceOp enum above
    const std::array<MPI_Op, 5> mpi_ops = {
        MPI_MAX,
        MPI_MIN,
        MPI_SUM,
        MPI_LOR,
        MPI_LAND
    };

    template<typename T>
    inline void Reduce (ReduceOp op, T* v, int cnt, int root, MPI_Comm comm)
    {
        auto mpi_op = mpi_ops[static_cast<int>(op)];
        Vector<T> tmp(v, v+cnt);
        if (root == -1) {
            // TODO: add BL_COMM_PROFILE commands
            MPI_Allreduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                          mpi_op, comm);
        } else {
            // TODO: add BL_COMM_PROFILE commands
            MPI_Reduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                       mpi_op, root, comm);
        }
    }

    template<typename T>
    inline void Reduce (ReduceOp op, T& v, int root, MPI_Comm comm) {
        Reduce(op, &v, 1, root, comm);
    }

    template<typename T>
    inline void Reduce (ReduceOp op, Vector<std::reference_wrapper<T> > const & v,
                        int root, MPI_Comm comm)
    {
        Vector<T> sndrcv(std::begin(v), std::end(v));
        Reduce(op, sndrcv.data(), v.size(), root, comm);
        for (int i = 0; i < v.size(); ++i) {
            v[i].get() = sndrcv[i];
        }
    }

    template<typename T>
    inline void Gather (const T* v, int cnt, T* vs, int root, MPI_Comm comm)
    {
        auto mpi_type = ParallelDescriptor::Mpi_typemap<T>::type();
        if (root == -1) {
            // TODO: check these BL_COMM_PROFILE commands
            BL_COMM_PROFILE(BLProfiler::Allgather, sizeof(T), BLProfiler::BeforeCall(),
                            BLProfiler::NoTag());
            // const_cast for MPI-2
            MPI_Allgather(const_cast<T*>(v), cnt, mpi_type, vs, cnt, mpi_type, comm);
            BL_COMM_PROFILE(BLProfiler::Allgather, sizeof(T), BLProfiler::AfterCall(),
                            BLProfiler::NoTag());
        } else {
            // TODO: add BL_COMM_PROFILE commands
            MPI_Gather(const_cast<T*>(v), cnt, mpi_type, vs, cnt, mpi_type, root, comm);
        }
    }

    template<typename T>
    inline void Gather (const T& v, T * vs, int root, MPI_Comm comm) {
        Gather(&v, 1, vs, root, comm);
    }

#else
    template<typename T> void Reduce (ReduceOp op, T* v, int cnt, int root, MPI_Comm comm) {}
    template<typename T> void Reduce (ReduceOp op, T& v, int root, MPI_Comm comm) {}
    template<typename T> void Reduce (ReduceOp op, Vector<std::reference_wrapper<T> > const & v, int root, MPI_Comm comm) {}

    template<typename T> void Gather (const T* v, int cnt, T* vs, int root, MPI_Comm comm) {}
    template<typename T> void Gather (const T& v, T * vs, int root, MPI_Comm comm) {}
#endif
}

namespace ParallelAllGather {
    template<typename T>
    void AllGather (const T* v, int cnt, T* vs, MPI_Comm comm) {
        detail::Gather(v, cnt, vs, -1, comm);
    }
    template<typename T>
    void AllGather (const T& v, T* vs, MPI_Comm comm) {
        detail::Gather(v, vs, -1, comm);
    }
}

namespace ParallelGather {
    template<typename T>
    void Gather (const T* v, int cnt, T* vs, int root, MPI_Comm comm) {
        detail::Gather(v, cnt, vs, root, comm);
    }
    template<typename T>
    void Gather (const T& v, T* vs, int root, MPI_Comm comm) {
        detail::Gather(v, vs, root, comm);
    }
}

namespace ParallelAllReduce {

    template<typename T>
    void Max (T& v, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::max, v, -1, comm);
    }
    template<typename T>
    void Max (T* v, int cnt, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::max, v, cnt, -1, comm);
    }
    template<typename T>
    void Max (Vector<std::reference_wrapper<T> > v, MPI_Comm comm) {
        detail::Reduce<T>(detail::ReduceOp::max, v, -1, comm);
    }

    template<typename T>
    void Min (T& v, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::min, v, -1, comm);
    }
    template<typename T>
    void Min (T* v, int cnt, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::min, v, cnt, -1, comm);
    }
    template<typename T>
    void Min (Vector<std::reference_wrapper<T> > v, MPI_Comm comm) {
        detail::Reduce<T>(detail::ReduceOp::min, v, -1, comm);
    }

    template<typename T>
    void Sum (T& v, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::sum, v, -1, comm);
    }
    template<typename T>
    void Sum (T* v, int cnt, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::sum, v, cnt, -1, comm);
    }
    template<typename T>
    void Sum (Vector<std::reference_wrapper<T> > v, MPI_Comm comm) {
        detail::Reduce<T>(detail::ReduceOp::sum, v, -1, comm);
    }

    inline void Or (bool & v, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        detail::Reduce(detail::ReduceOp::lor, iv, -1, comm);
        v = static_cast<bool>(iv);
    }

    inline void And (bool & v, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        detail::Reduce(detail::ReduceOp::land, iv, -1, comm);
        v = static_cast<bool>(iv);
    }
}

namespace ParallelReduce {

    template<typename T>
    void Max (T& v, int root, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::max, v, root, comm);
    }
    template<typename T>
    void Max (T* v, int cnt, int root, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::max, v, cnt, root, comm);
    }
    template<typename T>
    void Max (Vector<std::reference_wrapper<T> > v, int root, MPI_Comm comm) {
        detail::Reduce<T>(detail::ReduceOp::max, v, root, comm);
    }

    template<typename T>
    void Min (T& v, int root, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::min, v, root, comm);
    }
    template<typename T>
    void Min (T* v, int cnt, int root, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::min, v, cnt, root, comm);
    }
    template<typename T>
    void Min (Vector<std::reference_wrapper<T> > v, int root, MPI_Comm comm) {
        detail::Reduce<T>(detail::ReduceOp::min, v, root, comm);
    }

    template<typename T>
    void Sum (T& v, int root, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::sum, v, root, comm);
    }
    template<typename T>
    void Sum (T* v, int cnt, int root, MPI_Comm comm) {
        detail::Reduce(detail::ReduceOp::sum, v, cnt, root, comm);
    }
    template<typename T>
    void Sum (Vector<std::reference_wrapper<T> > v, int root, MPI_Comm comm) {
        detail::Reduce<T>(detail::ReduceOp::sum, v, root, comm);
    }

    inline void Or (bool & v, int root, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        detail::Reduce(detail::ReduceOp::lor, iv, root, comm);
        v = static_cast<bool>(iv);
    }

    inline void And (bool & v, int root, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        detail::Reduce(detail::ReduceOp::land, iv, root, comm);
        v = static_cast<bool>(iv);
    }
}

}

#endif
