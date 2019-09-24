#ifndef AMREX_ARRAY4_H_
#define AMREX_ARRAY4_H_

#include <AMReX_IntVect.H>

namespace amrex {

    template <typename T>
    struct Array4
    {
        T* AMREX_RESTRICT p;
        long jstride;
        long kstride;
        long nstride;
        Dim3 begin;
        Dim3 end;  // end is hi + 1
        int  ncomp;

        AMREX_GPU_HOST_DEVICE
        constexpr Array4 () noexcept : p(nullptr) {}

        template <class U=T, class = typename std::enable_if<std::is_const<U>::value>::type >
        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (Array4<typename std::remove_const<T>::type> const& rhs) noexcept
            : p(rhs.p),
              jstride(rhs.jstride),
              kstride(rhs.kstride),
              nstride(rhs.nstride),
              begin(rhs.begin),
              end(rhs.end),
              ncomp(rhs.ncomp)
            {}

        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (T* a_p, Dim3 const& a_begin, Dim3 const& a_end, int a_ncomp) noexcept
            : p(a_p),
              jstride(a_end.x-a_begin.x),
              kstride(jstride*(a_end.y-a_begin.y)),
              nstride(kstride*(a_end.z-a_begin.z)),
              begin(a_begin),
              end(a_end),
              ncomp(a_ncomp)
            {}

        template <class U,
                  class = typename std::enable_if
                  <std::is_same<typename std::remove_const<T>::type,
                                typename std::remove_const<U>::type>::value>::type >
        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (Array4<U> const& rhs, int start_comp) noexcept
            : p((T*)(rhs.p+start_comp*rhs.nstride)),
              jstride(rhs.jstride),
              kstride(rhs.kstride),
              nstride(rhs.nstride),
              begin(rhs.begin),
              end(rhs.end),
              ncomp(rhs.ncomp-start_comp)
        {}

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,0);
#endif
            return p[(i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride];
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,n);
#endif
            return p[(i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride+n*nstride];
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,0);
#endif
            return p + ((i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride);
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,0);
#endif
            return p + ((i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride+n*nstride);
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv[0],0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv[0],iv[1],0);
#else
            return this->operator()(iv[0],iv[1],iv[2]);
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv[0],iv[1],0,n);
#else
            return this->operator()(iv[0],iv[1],iv[2],n);
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv[0],0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv[0],iv[1],0);
#else
            return this->ptr(iv[0],iv[1],iv[2]);
#endif
        }

        template <class U=T, class = typename std::enable_if<!std::is_void<U>::value>::type >
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv[0],iv[1],0,n);
#else
            return this->ptr(iv[0],iv[1],iv[2],n);
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* dataPtr () const noexcept {
            return this->p;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::size_t size () const noexcept {
            return this->nstride * this->ncomp;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nComp () const noexcept { return ncomp; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (int i, int j, int k) const noexcept {
            return (i>=begin.x and i<end.x and j>=begin.y and j<end.y and k>=begin.z and k<end.z);
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        AMREX_GPU_HOST_DEVICE inline
        void index_assert (int i, int j, int k, int n) const
        {
            if (i<begin.x || i>=end.x || j<begin.y || j>=end.y || k<begin.z || k>=end.z
                || n < 0 || n >= ncomp) {
#ifdef AMREX_DEVICE_COMPILE
                std::printf(" (%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,0:%d)\n",
                            i, j, k, n, begin.x, end.x-1, begin.y, end.y-1,
                            begin.z, end.z-1, ncomp-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << i << "," << j << "," << k << "," <<  n
                   << ") is out of bound ("
                   << begin.x << ":" << end.x-1 << ","
                   << begin.y << ":" << end.y-1 << ","
                   << begin.z << ":" << end.z-1 << ","
                   << "0:" << ncomp-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
        }
#endif
    };

    template <class Tto, class Tfrom>
    AMREX_GPU_HOST_DEVICE
    Array4<Tto> ToArray4 (Array4<Tfrom> const& a_in) noexcept
    {
        return Array4<Tto>((Tto*)(a_in.p), a_in.begin, a_in.end, a_in.ncomp);
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 lbound (Array4<T> const& a) noexcept
    {
        return a.begin;
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 ubound (Array4<T> const& a) noexcept
    {
        return Dim3{a.end.x-1,a.end.y-1,a.end.z-1};
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 length (Array4<T> const& a) noexcept
    {
        return Dim3{a.end.x-a.begin.x,a.end.y-a.begin.y,a.end.z-a.begin.z};
    }

    template <typename T>
    std::ostream& operator<< (std::ostream& os, const Array4<T>& a) {
        os << "((" << lbound(a) << ',' << ubound(a) << ")," << a.ncomp << ')';
        return os;
    }

}

#endif
