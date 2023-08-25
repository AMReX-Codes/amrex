#ifndef AMREX_TABLE_DATA_H_
#define AMREX_TABLE_DATA_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_DataAllocator.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuPrint.H>

#include <cstring>
#include <iostream>
#include <sstream>
#include <type_traits>

namespace amrex {

template <typename T>
struct Table1D
{
    T* AMREX_RESTRICT p = nullptr;
    int begin = 1;
    int end = 0;

    AMREX_GPU_HOST_DEVICE
    constexpr Table1D () noexcept {}

    template <class U=T, typename std::enable_if<std::is_const<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table1D (Table1D<typename std::remove_const<T>::type> const& rhs) noexcept
        : p(rhs.p),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table1D (T* a_p, int a_begin, int a_end) noexcept
        : p(a_p),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, typename std::enable_if<!std::is_void<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i);
#endif
        return p[i-begin];
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i) const
    {
        if (i < begin || i >= end) {
#if AMREX_DEVICE_COMPILE
            AMREX_DEVICE_PRINTF(" (%d) is out of bound (%d:%d)\n",
                                i, begin, end-1);
            amrex::Abort();
#else
            std::stringstream ss;
            ss << " (" << i << ") is out of bound ("
               << begin << ":" << end-1 << ")";
            amrex::Abort(ss.str());
#endif
        }
    }
#endif
};

template <typename T>
struct Table2D
{
    T* AMREX_RESTRICT p = nullptr;
    Long jstride = 0;
    GpuArray<int,2> begin{{1,1}};
    GpuArray<int,2> end{{0,0}};

    AMREX_GPU_HOST_DEVICE
    constexpr Table2D () noexcept {}

    template <class U=T, typename std::enable_if<std::is_const<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table2D (Table2D<typename std::remove_const<T>::type> const& rhs) noexcept
        : p(rhs.p),
          jstride(rhs.jstride),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table2D (T* a_p,
                       GpuArray<int,2> const& a_begin,
                       GpuArray<int,2> const& a_end) noexcept
        : p(a_p),
          jstride(a_end[0]-a_begin[0]),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, typename std::enable_if<!std::is_void<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i, int j) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i,j);
#endif
        return p[(i-begin[0])+(j-begin[1])*jstride];
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i, int j) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1]) {
#if AMREX_DEVICE_COMPILE
            AMREX_DEVICE_PRINTF(" (%d,%d) is out of bound (%d:%d,%d:%d)\n",
                                i, j, begin[0], end[0]-1, begin[1], end[1]-1);
            amrex::Abort();
#else
            std::stringstream ss;
            ss << " (" << i << "," << j << ") is out of bound ("
               << begin[0] << ":" << end[0]-1
               << "," << begin[1] << ":" << end[1]-1 << ")";
            amrex::Abort(ss.str());
#endif
        }
    }
#endif
};

template <typename T>
struct Table3D
{
    T* AMREX_RESTRICT p = nullptr;
    Long jstride = 0;
    Long kstride = 0;
    GpuArray<int,3> begin{{1,1,1}};
    GpuArray<int,3> end{{0,0,0}};

    AMREX_GPU_HOST_DEVICE
    constexpr Table3D () noexcept {}

    template <class U=T, typename std::enable_if<std::is_const<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table3D (Table3D<typename std::remove_const<T>::type> const& rhs) noexcept
        : p(rhs.p),
          jstride(rhs.jstride),
          kstride(rhs.kstride),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table3D (T* a_p,
                       GpuArray<int,3> const& a_begin,
                       GpuArray<int,3> const& a_end) noexcept
        : p(a_p),
          jstride(a_end[0]-a_begin[0]),
          kstride(jstride*(a_end[1]-a_begin[1])),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, typename std::enable_if<!std::is_void<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i,j,k);
#endif
        return p[(i-begin[0])+(j-begin[1])*jstride+(k-begin[2])*kstride];
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i, int j, int k) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1] ||
            k < begin[2] || k >= end[2]) {
#if AMREX_DEVICE_COMPILE
            AMREX_DEVICE_PRINTF(" (%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d)\n",
                                i, j, k, begin[0], end[0]-1, begin[1], end[1]-1,
                                begin[2], end[2]-1);
            amrex::Abort();
#else
            std::stringstream ss;
            ss << " (" << i << "," << j << "," << k << ") is out of bound ("
               << begin[0] << ":" << end[0]-1
               << "," << begin[1] << ":" << end[1]-1
               << "," << begin[2] << ":" << end[2]-1 << ")";
            amrex::Abort(ss.str());
#endif
        }
    }
#endif
};

template <typename T>
struct Table4D
{
    T* AMREX_RESTRICT p = nullptr;
    Long jstride = 0;
    Long kstride = 0;
    Long nstride = 0;
    GpuArray<int,4> begin{{1,1,1,1}};
    GpuArray<int,4> end{{0,0,0,0}};

    AMREX_GPU_HOST_DEVICE
    constexpr Table4D () noexcept {}

    template <class U=T, typename std::enable_if<std::is_const<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table4D (Table4D<typename std::remove_const<T>::type> const& rhs) noexcept
        : p(rhs.p),
          jstride(rhs.jstride),
          kstride(rhs.kstride),
          nstride(rhs.nstride),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table4D (T* a_p,
                       GpuArray<int,4> const& a_begin,
                       GpuArray<int,4> const& a_end) noexcept
        : p(a_p),
          jstride(a_end[0]-a_begin[0]),
          kstride(jstride*(a_end[1]-a_begin[1])),
          nstride(kstride*(a_end[2]-a_begin[2])),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, typename std::enable_if<!std::is_void<U>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i,j,k,n);
#endif
        return p[(i-begin[0])+(j-begin[1])*jstride+(k-begin[2])*kstride+(n-begin[3])*nstride];
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i, int j, int k, int n) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1] ||
            k < begin[2] || k >= end[2] ||
            n < begin[3] || n >= end[3]) {
#if AMREX_DEVICE_COMPILE
            AMREX_DEVICE_PRINTF(" (%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,%d:%d)\n",
                                i, j, k, n, begin[0], end[0]-1, begin[1], end[1]-1,
                                begin[2], end[2]-1, begin[3], end[3]-1);
            amrex::Abort();
#else
            std::stringstream ss;
            ss << " (" << i << "," << j << "," << k << "," << n << ") is out of bound ("
               << begin[0] << ":" << end[0]-1
               << "," << begin[1] << ":" << end[1]-1
               << "," << begin[2] << ":" << end[2]-1
               << "," << begin[3] << ":" << end[3]-1 << ")";
            amrex::Abort(ss.str());
#endif
        }
    }
#endif
};

/**
 * \brief Multi-dimensional array class.
 *
 * This class is somewhat similar to FArrayBox/BaseFab. The main difference
 * is the dimension of the array in this class can be 1, 2, 3, or 4, whereas
 * the dimension of FArrayBox/BaseFab is the spatial dimension
 * (AMREX_SPACEDIM) plus a component dimension.  Below is an example of
 * using it to store a 3D table of data that is initialized on CPU and is
 * read-only by all GPU threads on the device.
 *
 * \code
 *      Array<int,3> tlo{0,0,0}; // lower bounds
 *      Array<int,3> thi{100,100,100}; // upper bounds
 *      TableData<Real,3> table_data(tlo, thi);
 *  #ifdef AMREX_USE_GPU
 *      TableData<Real,3> h_table_data(tlo, thi, The_Pinned_Arena());
 *      auto const& h_table = h_table_data.table();
 *  #else
 *      auto const& h_table = table_data.table();
 *  #endif
 *      // Initialize data on the host
 *      for (int k = tlo[0]; k <= thi[0]; ++k) {
 *      for (int j = tlo[1]; j <= thi[1]; ++j) {
 *      for (int i = tlo[2]; i <= thi[2]; ++i) {
 *          h_table(i,j,k) = i + 1.e3*j + 1.e6*k;
 *      }}}
 *  #ifdef AMREX_USE_GPU
 *      // Copy data to GPU memory
 *      table_data.copy(h_table_data);
 *      Gpu::streamSynchronize();  // not needed if the kernel using it is on the same stream
 *  #endif
 *      auto const& table = table_data.const_table(); // const makes it read only
 *      // We can now use table in device lambda.
 * \endcode
 */
template <typename T, int N>
class TableData
    : public DataAllocator
{
public:

    template <class U, int M> friend class TableData;
    using value_type = T;
    using table_type = std::conditional_t<N==1, Table1D<T>,
                       std::conditional_t<N==2, Table2D<T>,
                       std::conditional_t<N==3, Table3D<T>,
                                                Table4D<T> > > >;
    using const_table_type = std::conditional_t<N==1, Table1D<T const>,
                             std::conditional_t<N==2, Table2D<T const>,
                             std::conditional_t<N==3, Table3D<T const>,
                                                      Table4D<T const> > > >;

    TableData () noexcept;

    explicit TableData (Arena* ar) noexcept;

    TableData (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar = nullptr);

    TableData (TableData<T,N> const&) = delete;
    TableData<T,N>& operator= (TableData<T,N> const&) = delete;

    TableData (TableData<T,N>&& rhs) noexcept;
    TableData<T,N>& operator= (TableData<T,N> && rhs) noexcept;

    ~TableData () noexcept;

    constexpr int dim () const noexcept { return N; }

    void resize (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar = nullptr);

    Long size () const noexcept;

    Array<int,N> const& lo () const noexcept { return m_lo; }

    Array<int,N> const& hi () const noexcept { return m_hi; }

    void clear () noexcept;

    void copy (TableData<T,N> const& rhs) noexcept;

    table_type table () noexcept;
    const_table_type table () const noexcept;
    const_table_type const_table () const noexcept;

private:

    void define ();

    T* m_dptr = nullptr;
    Array<int,N> m_lo;
    Array<int,N> m_hi;
    Long m_truesize = 0L;
    bool m_ptr_owner = false;
};

template <typename T, int N>
TableData<T,N>::TableData () noexcept {}

template <typename T, int N>
TableData<T,N>::TableData (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar)
    : DataAllocator{ar}, m_lo(lo), m_hi(hi)
{
    define();
}


template <typename T, int N>
TableData<T,N>::TableData (TableData<T,N>&& rhs) noexcept
    : DataAllocator{rhs.arena()},
      m_dptr(rhs.m_dptr),
      m_lo(rhs.m_lo),
      m_hi(rhs.m_hi),
      m_truesize(rhs.m_truesize),
      m_ptr_owner(rhs.m_ptr_owner)
{
    rhs.m_dptr = nullptr;
    rhs.m_ptr_owner = false;
}

template <typename T, int N>
TableData<T,N>&
TableData<T,N>::operator= (TableData<T,N> && rhs) noexcept
{
    m_arena     = rhs.m_arena;
    m_dptr      = rhs.m_dptr;
    m_lo        = rhs.m_lo;
    m_hi        = rhs.m_hi;
    m_truesize  = rhs.m_truesize;
    m_ptr_owner = rhs.m_ptr_owner;
    rhs.m_dptr = nullptr;
    rhs.m_ptr_owner = false;
    return *this;
}

template <typename T, int N>
TableData<T,N>::~TableData () noexcept
{
    static_assert(std::is_trivially_copyable<T>() &&
                  std::is_trivially_destructible<T>(),
                  "TableData<T,N>: T must be trivially copyable and trivially destructible");
    static_assert(N>=1 && N <=4, "TableData<T,N>: N must be in the range of [1,4]");
    clear();
}

template <typename T, int N>
void
TableData<T,N>::resize (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar)
{
    m_lo = lo;
    m_hi = hi;

    if (ar == nullptr) {
        ar = m_arena;
    }

    if (arena() != DataAllocator(ar).arena()) {
        clear();
        m_arena = ar;
        define();
    } else if (m_dptr == nullptr || !m_ptr_owner) {
        m_dptr = nullptr;
        define();
    } else if (size() > m_truesize) {
        clear();
        define();
    }
}

template <typename T, int N>
Long
TableData<T,N>::size () const noexcept
{
    Long r = 1;
    for (int i = 0; i < N; ++i) {
        r *= m_hi[i] - m_lo[i] + 1;
    }
    return r;
}

template <typename T, int N>
void
TableData<T,N>::clear () noexcept
{
    if (m_dptr) {
        if (m_ptr_owner) {
            this->free(m_dptr);
        }
        m_dptr = nullptr;
        m_truesize = 0;
    }
}

template <typename T, int N>
void
TableData<T,N>::define ()
{
    m_truesize = size();
    AMREX_ASSERT(m_truesize >= 0);
    if (m_truesize == 0) {
        return;
    } else {
        m_ptr_owner = true;
        m_dptr = static_cast<T*>(this->alloc(m_truesize*sizeof(T)));
    }
}

namespace detail {
    template <typename T>
    Table1D<T> make_table (T* p, Array<int,1> const& lo, Array<int,1> const& hi) {
        return Table1D<T>(p, lo[0], hi[0]+1);
    }
    template <typename T>
    Table2D<T> make_table (T* p, Array<int,2> const& lo, Array<int,2> const& hi) {
        return Table2D<T>(p, {lo[0],lo[1]}, {hi[0]+1,hi[1]+1});
    }
    template <typename T>
    Table3D<T> make_table (T* p, Array<int,3> const& lo, Array<int,3> const& hi) {
        return Table3D<T>(p, {lo[0],lo[1],lo[2]}, {hi[0]+1,hi[1]+1,hi[2]+1});
    }
    template <typename T>
    Table4D<T> make_table (T* p, Array<int,4> const& lo, Array<int,4> const& hi) {
        return Table4D<T>(p, {lo[0],lo[1],lo[2],lo[3]}, {hi[0]+1,hi[1]+1,hi[2]+1,hi[3]+1});
    }
}

template <typename T, int N>
typename TableData<T,N>::table_type
TableData<T,N>::table () noexcept
{
    return detail::make_table<T>(m_dptr, m_lo, m_hi);
}

template <typename T, int N>
typename TableData<T,N>::const_table_type
TableData<T,N>::table () const noexcept
{
    return detail::make_table<T const>(m_dptr, m_lo, m_hi);
}

template <typename T, int N>
typename TableData<T,N>::const_table_type
TableData<T,N>::const_table () const noexcept
{
    return detail::make_table<T const>(m_dptr, m_lo, m_hi);
}

template <typename T, int N>
void
TableData<T,N>::copy (TableData<T,N> const& rhs) noexcept
{
    std::size_t count = sizeof(T)*size();
#ifdef AMREX_USE_GPU
    bool this_on_device = arena()->isManaged() || arena()->isDevice();
    bool rhs_on_device = rhs.arena()->isManaged() || rhs.arena()->isDevice();
    if (this_on_device && rhs_on_device) {
        Gpu::dtod_memcpy_async(m_dptr, rhs.m_dptr, count);
    } else if (this_on_device && !rhs_on_device) {
        Gpu::htod_memcpy_async(m_dptr, rhs.m_dptr, count);
    } else if (!this_on_device && rhs_on_device) {
        Gpu::dtoh_memcpy_async(m_dptr, rhs.m_dptr, count);
    } else
#endif
    {
        std::memcpy(m_dptr, rhs.m_dptr, count);
    }
}

}

#endif
