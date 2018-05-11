#ifndef AMREX_PARALLEL_REDUCE_H_
#define AMREX_PARALLEL_REDUCE_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <type_traits>

namespace amrex {

namespace util {

#ifdef BL_USE_MPI

    inline void GccPlacater ()
    {
        std::allocator<int> foo;
        std::allocator<Real> bar;
    }

    template<typename T>
    inline void Reduce (MPI_Op op, T* v, int cnt, int root, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        if (root == -1) {
            MPI_Allreduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                          op, comm);
        } else {
            MPI_Reduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                       op, root, comm);
        }
    }

    template<typename T>
    inline void Reduce (MPI_Op op, T& v, int root, MPI_Comm comm) {
        Reduce(op, &v, 1, root, comm);
    }

    template<typename T>
    inline void Reduce (MPI_Op op, Vector<std::reference_wrapper<T> >&& v,
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
            MPI_Allgather(v, cnt, mpi_type, vs, cnt, mpi_type, comm);
        } else {
            MPI_Gather(v, cnt, mpi_type, vs, cnt, mpi_type, root, comm);
        }
    }

    template<typename T>
    inline void Gather (const T& v, T * vs, int root, MPI_Comm comm) {
        Gather(&v, 1, vs, root, comm);
    }

#else

    template<typename T> void Reduce (MPI_Op op, T* v, int cnt, int root, MPI_Comm comm) {}
    template<typename T> void Reduce (MPI_Op op, T& v, int root, MPI_Comm comm) {}
    template<typename T> void Reduce (MPI_Op op, Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {}

    template<typename T> void Gather (const T* v, int cnt, T* vs, int root, MPI_Comm comm) {}
    template<typename T> void Gather (const T& v, T * vs, int root, MPI_Comm comm) {}
#endif
}

namespace ParallelAllGather {
    template<typename T>
    void AllGather (const T* v, int cnt, T* vs, MPI_Comm comm) {
        util::Gather(v, cnt, vs, -1, comm);
    }
    template<typename T>
    void AllGather (const T& v, T* vs, MPI_Comm comm) {
        util::Gather(v, vs, -1, comm);
    }
}

namespace ParallelGather {
    template<typename T>
    void Gather (const T* v, int cnt, T* vs, int root, MPI_Comm comm) {
        util::Gather(v, cnt, vs, root, comm);
    }
    template<typename T>
    void Gather (const T& v, T* vs, int root, MPI_Comm comm) {
        util::Gather(v, vs, root, comm);
    }
}

namespace ParallelAllReduce {

    template<typename T>
    void Max (T& v, MPI_Comm comm) {
        util::Reduce(MPI_MAX, v, -1, comm);
    }
    template<typename T>
    void Max (T* v, int cnt, MPI_Comm comm) {
        util::Reduce(MPI_MAX, v, cnt, -1, comm);
    }
    template<typename T>
    void Max (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {
        util::Reduce(MPI_MAX, std::move(v), -1, comm);
    }

    template<typename T>
    void Min (T& v, MPI_Comm comm) {
        util::Reduce(MPI_MIN, v, -1, comm);
    }
    template<typename T>
    void Min (T* v, int cnt, MPI_Comm comm) {
        util::Reduce(MPI_MIN, v, cnt, -1, comm);
    }
    template<typename T>
    void Min (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {
        util::Reduce(MPI_MIN, std::move(v), -1, comm);
    }

    template<typename T>
    void Sum (T& v, MPI_Comm comm) {
        util::Reduce(MPI_SUM, v, -1, comm);
    }
    template<typename T>
    void Sum (T* v, int cnt, MPI_Comm comm) {
        util::Reduce(MPI_SUM, v, cnt, -1, comm);
    }
    template<typename T>
    void Sum (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {
        util::Reduce(MPI_SUM, std::move(v), -1, comm);
    }

    inline void Or (bool & v, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        util::Reduce(MPI_LOR, iv, -1, comm);
        v = static_cast<bool>(iv);
    }

    inline void And (bool & v, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        util::Reduce(MPI_LAND, iv, -1, comm);
        v = static_cast<bool>(iv);
    }
}

namespace ParallelReduce {

    template<typename T>
    void Max (T& v, int root, MPI_Comm comm) {
        util::Reduce(MPI_MAX, v, root, comm);
    }
    template<typename T>
    void Max (T* v, int cnt, int root, MPI_Comm comm) {
        util::Reduce(MPI_MAX, v, cnt, root, comm);
    }
    template<typename T>
    void Max (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {
        util::Reduce(MPI_MAX, std::move(v), root, comm);
    }

    template<typename T>
    void Min (T& v, int root, MPI_Comm comm) {
        util::Reduce(MPI_MIN, v, root, comm);
    }
    template<typename T>
    void Min (T* v, int cnt, int root, MPI_Comm comm) {
        util::Reduce(MPI_MIN, v, cnt, root, comm);
    }
    template<typename T>
    void Min (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {
        util::Reduce(MPI_MIN, std::move(v), root, comm);
    }

    template<typename T>
    void Sum (T& v, int root, MPI_Comm comm) {
        util::Reduce(MPI_SUM, v, root, comm);
    }
    template<typename T>
    void Sum (T* v, int cnt, int root, MPI_Comm comm) {
        util::Reduce(MPI_SUM, v, cnt, root, comm);
    }
    template<typename T>
    void Sum (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {
        util::Reduce(MPI_SUM, std::move(v), root, comm);
    }

    inline void Or (bool & v, int root, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        util::Reduce(MPI_LOR, iv, root, comm);
        v = static_cast<bool>(iv);
    }

    inline void And (bool & v, int root, MPI_Comm comm) {
        auto iv = static_cast<int>(v);
        util::Reduce(MPI_LAND, iv, root, comm);
        v = static_cast<bool>(iv);
    }
}

}

#endif
