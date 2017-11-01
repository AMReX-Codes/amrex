#ifndef AMREX_PARALLEL_REDUCE_H_
#define AMREX_PARALLEL_REDUCE_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

namespace amrex {

namespace ParallelAllReduce {

#ifdef BL_USE_MPI

    template<typename T>
    void Max (T& v, MPI_Comm comm)
    {
        const T tmp = v;
        MPI_Allreduce(&tmp, &v, 1, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MAX, comm);
    }

    template<typename T>
    void Max (T* v, int cnt, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        MPI_Allreduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MAX, comm);
    }

    template<typename T>
    void Max (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm)
    {
        const Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Allreduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MAX, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

    template<typename T>
    void Sum (T& v, MPI_Comm comm)
    {
        const T tmp = v;
        MPI_Allreduce(&tmp, &v, 1, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_SUM, comm);
    }

    template<typename T>
    void Sum (T* v, int cnt, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        MPI_Allreduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_SUM, comm);
    }

    template<typename T>
    void Sum (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm)
    {
        const Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Allreduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_SUM, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

#else

    template<typename T> void Max (T& rvar, MPI_Comm comm) {}
    template<typename T> void Max (T* rvar, int cnt, MPI_Comm comm) {}
    template<typename T> void Max (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {}

    template<typename T> void Sum (T& rvar, MPI_Comm comm) {}
    template<typename T> void Sum (T* rvar, int cnt, MPI_Comm comm) {}
    template<typename T> void Sum (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {}

#endif
}

}

#endif
