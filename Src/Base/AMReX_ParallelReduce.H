#ifndef AMREX_PARALLEL_REDUCE_H_
#define AMREX_PARALLEL_REDUCE_H_

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <type_traits>

namespace amrex {

namespace ParallelAllReduce {

#ifdef BL_USE_MPI

    template<typename T>
    void Max (T& v, MPI_Comm comm)
    {
        T tmp = v;
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
        Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Allreduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MAX, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

    template<typename T>
    void Min (T& v, MPI_Comm comm)
    {
        T tmp = v;
        MPI_Allreduce(&tmp, &v, 1, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MIN, comm);
    }

    template<typename T>
    void Min (T* v, int cnt, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        MPI_Allreduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MIN, comm);
    }

    template<typename T>
    void Min (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm)
    {
        Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Allreduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                      MPI_MIN, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

    template<typename T>
    void Sum (T& v, MPI_Comm comm)
    {
        T tmp = v;
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
        Vector<T> snd(std::begin(v), std::end(v));
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

    template<typename T> void Min (T& rvar, MPI_Comm comm) {}
    template<typename T> void Min (T* rvar, int cnt, MPI_Comm comm) {}
    template<typename T> void Min (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {}

    template<typename T> void Sum (T& rvar, MPI_Comm comm) {}
    template<typename T> void Sum (T* rvar, int cnt, MPI_Comm comm) {}
    template<typename T> void Sum (Vector<std::reference_wrapper<T> >&& v, MPI_Comm comm) {}

#endif
}

namespace ParallelReduce {

#ifdef BL_USE_MPI

    template<typename T>
    void Max (T& v, int root, MPI_Comm comm)
    {
        T tmp = v;
        MPI_Reduce(&tmp, &v, 1, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_MAX, root, comm);
    }

    template<typename T>
    void Max (T* v, int cnt, int root, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        MPI_Reduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_MAX, root, comm);
    }

    template<typename T>
    void Max (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm)
    {
        Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Reduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_MAX, root, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

    template<typename T>
    void Min (T& v, int root, MPI_Comm comm)
    {
        T tmp = v;
        MPI_Reduce(&tmp, &v, 1, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_MIN, root, comm);
    }

    template<typename T>
    void Min (T* v, int cnt, int root, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        MPI_Reduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_MIN, root, comm);
    }

    template<typename T>
    void Min (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm)
    {
        Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Reduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_MIN, root, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

    template<typename T>
    void Sum (T& v, int root, MPI_Comm comm)
    {
        T tmp = v;
        MPI_Reduce(&tmp, &v, 1, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_SUM, root, comm);
    }

    template<typename T>
    void Sum (T* v, int cnt, int root, MPI_Comm comm)
    {
        Vector<T> tmp(v, v+cnt);
        MPI_Reduce(tmp.data(), v, cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_SUM, root, comm);
    }

    template<typename T>
    void Sum (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm)
    {
        Vector<T> snd(std::begin(v), std::end(v));
        const int cnt = v.size();
        Vector<T> rcv(cnt);
        MPI_Reduce(snd.data(), rcv.data(), cnt, ParallelDescriptor::Mpi_typemap<T>::type(),
                   MPI_SUM, root, comm);
        for (int i = 0; i < cnt; ++i) {
            v[i].get() = rcv[i];
        }
    }

#else

    template<typename T> void Max (T& rvar, int root, MPI_Comm comm) {}
    template<typename T> void Max (T* rvar, int cnt, int root, MPI_Comm comm) {}
    template<typename T> void Max (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {}

    template<typename T> void Min (T& rvar, int root, MPI_Comm comm) {}
    template<typename T> void Min (T* rvar, int cnt, int root, MPI_Comm comm) {}
    template<typename T> void Min (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {}

    template<typename T> void Sum (T& rvar, int root, MPI_Comm comm) {}
    template<typename T> void Sum (T* rvar, int cnt, int root, MPI_Comm comm) {}
    template<typename T> void Sum (Vector<std::reference_wrapper<T> >&& v, int root, MPI_Comm comm) {}

#endif
}

}

#endif
