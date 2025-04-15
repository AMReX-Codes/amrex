#ifndef AMREX_ARRAY4_H_
#define AMREX_ARRAY4_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_GpuPrint.H>

#include <iostream>
#include <sstream>

namespace amrex {

    template <typename T>
    struct CellData // Data in a single cell
    {
        T* AMREX_RESTRICT p = nullptr;
        Long stride = 0;
        int ncomp = 0;

        AMREX_GPU_HOST_DEVICE
        constexpr CellData (T* a_p, Long a_stride, int a_ncomp)
            : p(a_p), stride(a_stride), ncomp(a_ncomp)
            {}

        template <class U=T,
                  std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr CellData (CellData<std::remove_const_t<T>> const& rhs) noexcept
            : p(rhs.p), stride(rhs.stride), ncomp(rhs.ncomp)
            {}

        AMREX_GPU_HOST_DEVICE
        explicit operator bool() const noexcept { return p != nullptr; }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE
        int nComp() const noexcept { return ncomp; }

        template <class U=T,
                  std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator[] (int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            if (n < 0 || n >= ncomp) {
                AMREX_IF_ON_DEVICE((
                    AMREX_DEVICE_PRINTF(" %d is out of bound (0:%d)", n, ncomp-1);
                ))
                AMREX_IF_ON_HOST((
                    std::stringstream ss;
                    ss << " " << n << " is out of bound: (0:" << ncomp-1 << ")";
                    amrex::Abort(ss.str());
                ))
            }
#endif
            return p[n*stride];
        }
    };

    template <typename T>
    struct Array4
    {
        T* AMREX_RESTRICT p;
        Dim3 begin{1,1,1};
        Dim3 len{0,0,0};
        int  ncomp=0;

        AMREX_GPU_HOST_DEVICE
        constexpr Array4 () noexcept : p(nullptr) {}

        template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (Array4<std::remove_const_t<T>> const& rhs) noexcept
            : p(rhs.p),
              begin(rhs.begin),
              len(rhs.len),
              ncomp(rhs.ncomp)
            {}

        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (T* a_p, Dim3 const& a_begin, Dim3 const& a_end, int a_ncomp) noexcept
            : p(a_p),
              begin(a_begin),
              len{a_end.x-a_begin.x, a_end.y-a_begin.y, a_end.z-a_begin.z},
              ncomp(a_ncomp)
            {}

        template <class U,
                  std::enable_if_t
                  <std::is_same_v<std::remove_const_t<T>,
                                std::remove_const_t<U>>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (Array4<U> const& rhs, int start_comp) noexcept
            : p((T*)(rhs.p+start_comp*rhs.nstride())),
              begin(rhs.begin),
              len(rhs.len),
              ncomp(rhs.ncomp-start_comp)
        {}

        template <class U,
                  std::enable_if_t
                  <std::is_same_v<std::remove_const_t<T>,
                                std::remove_const_t<U>>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (Array4<U> const& rhs, int start_comp, int num_comps) noexcept
            : p((T*)(rhs.p+start_comp*rhs.nstride())),
              begin(rhs.begin),
              len(rhs.len),
              ncomp(num_comps)
        {}

        AMREX_GPU_HOST_DEVICE
        explicit operator bool() const noexcept { return p != nullptr; }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Long jstride () const noexcept { return Long(len.x); }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Long kstride () const noexcept { return Long(len.x)*Long(len.y); }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Long nstride () const noexcept { return Long(len.x)*Long(len.y)*Long(len.z); }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,0);
#endif
            return p[(i-begin.x)+(j-begin.y)*jstride()+(k-begin.z)*kstride()];
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,n);
#endif
            return p[(i-begin.x)+Long(len.x)*((j-begin.y)+Long(len.y)*((k-begin.z)+Long(len.z)*n))];
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,0);
#endif
            return p + ((i-begin.x)+(j-begin.y)*jstride()+(k-begin.z)*kstride());
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,n);
#endif
            return p + ((i-begin.x)+(j-begin.y)*jstride()+(k-begin.z)*kstride()+n*nstride());
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv[0],0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv[0],iv[1],0);
#else
            return this->operator()(iv[0],iv[1],iv[2]);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv[0],iv[1],0,n);
#else
            return this->operator()(iv[0],iv[1],iv[2],n);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv[0],0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv[0],iv[1],0);
#else
            return this->ptr(iv[0],iv[1],iv[2]);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv[0],iv[1],0,n);
#else
            return this->ptr(iv[0],iv[1],iv[2],n);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z);
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell, int n) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z,n);
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (Dim3 const& cell) const noexcept {
            return this->ptr(cell.x,cell.y,cell.z);
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (Dim3 const& cell, int n) const noexcept {
            return this->ptr(cell.x,cell.y,cell.z,n);
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* dataPtr () const noexcept {
            return this->p;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::size_t size () const noexcept {
            return this->nstride() * this->ncomp;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nComp () const noexcept { return ncomp; }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (int i, int j, int k) const noexcept {
            return (i>=begin.x && i<(begin.x+len.x) && j>=begin.y && j<(begin.y+len.y) && k>=begin.z && k<(begin.z+len.z));
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (IntVect const& iv) const noexcept {
            return AMREX_D_TERM(   iv[0]>=begin.x && iv[0]<(begin.x+len.x),
                                && iv[1]>=begin.y && iv[1]<(begin.y+len.y),
                                && iv[2]>=begin.z && iv[2]<(begin.z+len.z));
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (Dim3 const& cell) const noexcept {
            return this->contains(cell.x,cell.y,cell.z);
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE inline
#endif
        void index_assert (int i, int j, int k, int n) const
        {
            if (i<begin.x || i>=(begin.x+len.x) || j<begin.y || j>=(begin.y+len.y) || k<begin.z || k>=(begin.z+len.z)
                || n < 0 || n >= ncomp) {
                AMREX_IF_ON_DEVICE((
                    AMREX_DEVICE_PRINTF(" (%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,0:%d)\n",
                                        i, j, k, n, begin.x, (begin.x+len.x)-1, begin.y, (begin.y+len.y)-1,
                                        begin.z, (begin.z+len.z)-1, ncomp-1);
                    amrex::Abort();
                ))
                AMREX_IF_ON_HOST((
                    std::stringstream ss;
                    ss << " (" << i << "," << j << "," << k << "," <<  n
                    << ") is out of bound ("
                    << begin.x << ":" << (begin.x+len.x)-1 << ","
                    << begin.y << ":" << (begin.y+len.y)-1 << ","
                    << begin.z << ":" << (begin.z+len.z)-1 << ","
                    << "0:" << ncomp-1 << ")";
                    amrex::Abort(ss.str());
                ))
            }
        }
#endif

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CellData<T> cellData (int i, int j, int k) const noexcept {
            return CellData<T>{this->ptr(i,j,k), nstride(), ncomp};
        }
    };

    template <class Tto, class Tfrom>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Array4<Tto> ToArray4 (Array4<Tfrom> const& a_in) noexcept
    {
        return Array4<Tto>((Tto*)(a_in.p), a_in.begin,
                           Dim3{a_in.begin.x + a_in.len.x,
                                a_in.begin.y + a_in.len.y,
                                a_in.begin.z + a_in.len.z},
                           a_in.ncomp);
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 lbound (Array4<T> const& a) noexcept
    {
        return a.begin;
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 ubound (Array4<T> const& a) noexcept
    {
        return Dim3{a.begin.x+a.len.x-1,a.begin.y+a.len.y-1,a.begin.z+a.len.z-1};
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 begin (Array4<T> const& a) noexcept
    {
        return a.begin;
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 end (Array4<T> const& a) noexcept
    {
        return Dim3{a.begin.x+a.len.x,a.begin.y+a.len.y,a.begin.z+a.len.z};
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 length (Array4<T> const& a) noexcept
    {
        return a.len;
    }

    template <typename T>
    std::ostream& operator<< (std::ostream& os, const Array4<T>& a) {
        os << "((" << lbound(a) << ',' << ubound(a) << ")," << a.ncomp << ')';
        return os;
    }

    //
    // Type traits for detecting if a class has a size() constexpr function.
    //
    template <class A, class Enable = void> struct HasMultiComp : std::false_type {};
    //
    template <class B>
    struct HasMultiComp<B, std::enable_if_t<B().size() >= 1>>
        : std::true_type {};

    //
    // PolymorphicArray4 can be used to access both AoS and SoA with
    // (i,j,k,n).  Here SoA refers multi-component BaseFab, and AoS refers
    // to single-component BaseFab of multi-component GpuArray.
    //
    template <typename T>
    struct PolymorphicArray4
        : public Array4<T>
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        PolymorphicArray4 (Array4<T> const& a)
            : Array4<T>{a} {}

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) const noexcept {
            return this->Array4<T>::operator()(i,j,k);
        }

        template <class U=T, std::enable_if_t< amrex::HasMultiComp<U>::value,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        typename U::reference_type
        operator() (int i, int j, int k, int n) const noexcept {
            return this->Array4<T>::operator()(i,j,k,0)[n];
        }

        template <class U=T, std::enable_if_t<!amrex::HasMultiComp<U>::value,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k, int n) const noexcept {
            return this->Array4<T>::operator()(i,j,k,n);
        }
    };

    template <typename T>
    [[nodiscard]] PolymorphicArray4<T>
    makePolymorphic (Array4<T> const& a)
    {
        return PolymorphicArray4<T>(a);
    }
}

#endif
