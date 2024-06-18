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

    constexpr Table1D () noexcept = default;

    template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table1D (Table1D<std::remove_const_t<T>> const& rhs) noexcept
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

    template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
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
            AMREX_IF_ON_DEVICE((
                AMREX_DEVICE_PRINTF(" (%d) is out of bound (%d:%d)\n",
                                    i, begin, end-1);
                amrex::Abort();
            ))
            AMREX_IF_ON_HOST((
                std::stringstream ss;
                ss << " (" << i << ") is out of bound ("
                << begin << ":" << end-1 << ")";
                amrex::Abort(ss.str());
            ))
        }
    }
#endif
};

template <typename T, typename ORDER = Order::F>
struct Table2D
{
    T* AMREX_RESTRICT p = nullptr;
    Long stride1 = 0;
    GpuArray<int,2> begin{{1,1}};
    GpuArray<int,2> end{{0,0}};

    constexpr Table2D () noexcept = default;

    template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table2D (Table2D<std::remove_const_t<T>, ORDER> const& rhs) noexcept
        : p(rhs.p),
          stride1(rhs.stride1),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table2D (T* a_p,
                       GpuArray<int,2> const& a_begin,
                       GpuArray<int,2> const& a_end) noexcept
        : p(a_p),
          stride1(len0(a_begin,a_end)),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i, int j) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i,j);
#endif
        static_assert(std::is_same_v<ORDER,Order::F> ||
                      std::is_same_v<ORDER,Order::C>);
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return p[(i-begin[0])+(j-begin[1])*stride1];
        } else {
            return p[(i-begin[0])*stride1+(j-begin[1])];
        }
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i, int j) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1]) {
            AMREX_IF_ON_DEVICE((
                AMREX_DEVICE_PRINTF(" (%d,%d) is out of bound (%d:%d,%d:%d)\n",
                                    i, j, begin[0], end[0]-1, begin[1], end[1]-1);
                amrex::Abort();
            ))
            AMREX_IF_ON_HOST((
                std::stringstream ss;
                ss << " (" << i << "," << j << ") is out of bound ("
                << begin[0] << ":" << end[0]-1
                << "," << begin[1] << ":" << end[1]-1 << ")";
                amrex::Abort(ss.str());
            ))
        }
    }
#endif

private:

    static constexpr int len0 (GpuArray<int,2> const& a_begin,
                               GpuArray<int,2> const& a_end) noexcept
    {
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return a_end[0] - a_begin[0];
        } else {
            return a_end[1] - a_begin[1];
        }
    }
};

template <typename T, typename ORDER = Order::F>
struct Table3D
{
    T* AMREX_RESTRICT p = nullptr;
    Long stride1 = 0;
    Long stride2 = 0;
    GpuArray<int,3> begin{{1,1,1}};
    GpuArray<int,3> end{{0,0,0}};

    constexpr Table3D () noexcept = default;

    template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table3D (Table3D<std::remove_const_t<T>,ORDER> const& rhs) noexcept
        : p(rhs.p),
          stride1(rhs.stride1),
          stride2(rhs.stride2),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table3D (T* a_p,
                       GpuArray<int,3> const& a_begin,
                       GpuArray<int,3> const& a_end) noexcept
        : p(a_p),
          stride1(        len0(a_begin,a_end)),
          stride2(stride1*len1(a_begin,a_end)),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i,j,k);
#endif
        static_assert(std::is_same_v<ORDER,Order::F> ||
                      std::is_same_v<ORDER,Order::C>);
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return p[(i-begin[0])+(j-begin[1])*stride1+(k-begin[2])*stride2];
        } else {
            return p[(i-begin[0])*stride2+(j-begin[1])*stride1+(k-begin[2])];
        }
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i, int j, int k) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1] ||
            k < begin[2] || k >= end[2]) {
            AMREX_IF_ON_DEVICE((
                AMREX_DEVICE_PRINTF(" (%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d)\n",
                                    i, j, k, begin[0], end[0]-1, begin[1], end[1]-1,
                                    begin[2], end[2]-1);
                amrex::Abort();
            ))
            AMREX_IF_ON_HOST((
                std::stringstream ss;
                ss << " (" << i << "," << j << "," << k << ") is out of bound ("
                << begin[0] << ":" << end[0]-1
                << "," << begin[1] << ":" << end[1]-1
                << "," << begin[2] << ":" << end[2]-1 << ")";
                amrex::Abort(ss.str());
            ))
        }
    }
#endif

private:

    static constexpr int len0 (GpuArray<int,3> const& a_begin,
                               GpuArray<int,3> const& a_end) noexcept
    {
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return a_end[0] - a_begin[0];
        } else {
            return a_end[2] - a_begin[2];
        }
    }

    static constexpr int len1 (GpuArray<int,3> const& a_begin,
                               GpuArray<int,3> const& a_end) noexcept
    {
        return a_end[1] - a_begin[1];
    }
};

template <typename T, typename ORDER = Order::F>
struct Table4D
{
    T* AMREX_RESTRICT p = nullptr;
    Long stride1 = 0;
    Long stride2 = 0;
    Long stride3 = 0;
    GpuArray<int,4> begin{{1,1,1,1}};
    GpuArray<int,4> end{{0,0,0,0}};

    constexpr Table4D () noexcept = default;

    template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr Table4D (Table4D<std::remove_const_t<T>,ORDER> const& rhs) noexcept
        : p(rhs.p),
          stride1(rhs.stride1),
          stride2(rhs.stride2),
          stride3(rhs.stride3),
          begin(rhs.begin),
          end(rhs.end)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr Table4D (T* a_p,
                       GpuArray<int,4> const& a_begin,
                       GpuArray<int,4> const& a_end) noexcept
        : p(a_p),
          stride1(        len0(a_begin,a_end)),
          stride2(stride1*len1(a_begin,a_end)),
          stride3(stride2*len2(a_begin,a_end)),
          begin(a_begin),
          end(a_end)
        {}

    AMREX_GPU_HOST_DEVICE
    explicit operator bool () const noexcept { return p != nullptr; }

    template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        index_assert(i,j,k,n);
#endif
        static_assert(std::is_same_v<ORDER,Order::F> ||
                      std::is_same_v<ORDER,Order::C>);
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return p[(i-begin[0])+(j-begin[1])*stride1+(k-begin[2])*stride2+(n-begin[3])*stride3];
        } else {
            return p[(i-begin[0])*stride3+(j-begin[1])*stride2+(k-begin[2])*stride1+(n-begin[3])];
        }
    }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
    AMREX_GPU_HOST_DEVICE inline
    void index_assert (int i, int j, int k, int n) const
    {
        if (i < begin[0] || i >= end[0] ||
            j < begin[1] || j >= end[1] ||
            k < begin[2] || k >= end[2] ||
            n < begin[3] || n >= end[3]) {
            AMREX_IF_ON_DEVICE((
                AMREX_DEVICE_PRINTF(" (%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,%d:%d)\n",
                                    i, j, k, n, begin[0], end[0]-1, begin[1], end[1]-1,
                                    begin[2], end[2]-1, begin[3], end[3]-1);
                amrex::Abort();
            ))
            AMREX_IF_ON_HOST((
                std::stringstream ss;
                ss << " (" << i << "," << j << "," << k << "," << n << ") is out of bound ("
                << begin[0] << ":" << end[0]-1
                << "," << begin[1] << ":" << end[1]-1
                << "," << begin[2] << ":" << end[2]-1
                << "," << begin[3] << ":" << end[3]-1 << ")";
                amrex::Abort(ss.str());
            ))
        }
    }
#endif

private:

    static constexpr int len0 (GpuArray<int,4> const& a_begin,
                               GpuArray<int,4> const& a_end) noexcept
    {
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return a_end[0] - a_begin[0];
        } else {
            return a_end[3] - a_begin[3];
        }
    }

    static constexpr int len1 (GpuArray<int,4> const& a_begin,
                               GpuArray<int,4> const& a_end) noexcept
    {
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return a_end[1] - a_begin[1];
        } else {
            return a_end[2] - a_begin[2];
        }
    }

    static constexpr int len2 (GpuArray<int,4> const& a_begin,
                               GpuArray<int,4> const& a_end) noexcept
    {
        if constexpr (std::is_same_v<ORDER,Order::F>) {
            return a_end[2] - a_begin[2];
        } else {
            return a_end[1] - a_begin[1];
        }
    }
};

/**
 * \brief Multi-dimensional array class.
 *
 * This class is somewhat similar to FArrayBox/BaseFab. The main difference
 * is the dimension of the array in this class can be 1, 2, 3, or 4, whereas
 * the dimension of FArrayBox/BaseFab is the spatial dimension
 * (AMREX_SPACEDIM) plus a component dimension. Another difference is that
 * this class supports both column-major order (i.e., Fortran order) and
 * row-major order (i.e., C order), whereas FArrayBox/BaseFab is always
 * column-major. Below is an example of using it to store a 3D table of data
 * that is initialized on CPU and is read-only by all GPU threads on the
 * device.
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
template <typename T, int N, typename ORDER = Order::F>
class TableData
    : public DataAllocator
{
public:

    template <class U, int M, class O> friend class TableData;
    using value_type = T;
    using table_type = std::conditional_t<N==1, Table1D<T>,
                       std::conditional_t<N==2, Table2D<T, ORDER>,
                       std::conditional_t<N==3, Table3D<T, ORDER>,
                                                Table4D<T, ORDER> > > >;
    using const_table_type = std::conditional_t<N==1, Table1D<T const>,
                             std::conditional_t<N==2, Table2D<T const, ORDER>,
                             std::conditional_t<N==3, Table3D<T const, ORDER>,
                                                      Table4D<T const, ORDER> > > >;

    TableData () noexcept = default;

    explicit TableData (Arena* ar) noexcept;

    TableData (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar = nullptr);

    TableData (TableData<T,N,ORDER> const&) = delete;
    TableData<T,N,ORDER>& operator= (TableData<T,N,ORDER> const&) = delete;

    TableData (TableData<T,N,ORDER>&& rhs) noexcept;
    TableData<T,N,ORDER>& operator= (TableData<T,N,ORDER> && rhs) noexcept;

    ~TableData () noexcept;

    [[nodiscard]] constexpr int dim () const noexcept { return N; }

    void resize (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar = nullptr);

    [[nodiscard]] Long size () const noexcept;

    Array<int,N> const& lo () const noexcept { return m_lo; }

    Array<int,N> const& hi () const noexcept { return m_hi; }

    void clear () noexcept;

    void copy (TableData<T,N,ORDER> const& rhs) noexcept;

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

template <typename T, int N, typename ORDER>
TableData<T,N,ORDER>::TableData (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar)
    : DataAllocator{ar}, m_lo(lo), m_hi(hi)
{
    define();
}


template <typename T, int N, typename ORDER>
TableData<T,N,ORDER>::TableData (TableData<T,N,ORDER>&& rhs) noexcept
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

template <typename T, int N, typename ORDER>
TableData<T,N,ORDER>&
TableData<T,N,ORDER>::operator= (TableData<T,N,ORDER> && rhs) noexcept
{
    if (this != &rhs) {
        clear();
        m_arena     = rhs.m_arena;
        m_dptr      = rhs.m_dptr;
        m_lo        = rhs.m_lo;
        m_hi        = rhs.m_hi;
        m_truesize  = rhs.m_truesize;
        m_ptr_owner = rhs.m_ptr_owner;
        rhs.m_dptr = nullptr;
        rhs.m_ptr_owner = false;
    }
    return *this;
}

template <typename T, int N, typename ORDER>
TableData<T,N,ORDER>::~TableData () noexcept
{
    static_assert(std::is_trivially_copyable<T>() &&
                  std::is_trivially_destructible<T>(),
                  "TableData<T,N,ORDER>: T must be trivially copyable and trivially destructible");
    static_assert(N>=1 && N <=4, "TableData<T,N,ORDER>: N must be in the range of [1,4]");
    static_assert(std::is_same_v<ORDER,Order::F> ||
                  std::is_same_v<ORDER,Order::C>,
                  "TableDat<T,N,ORDER>: ORDER must be either Order::F or Order::C");
    clear();
}

template <typename T, int N, typename ORDER>
void
TableData<T,N,ORDER>::resize (Array<int,N> const& lo, Array<int,N> const& hi, Arena* ar)
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

template <typename T, int N, typename ORDER>
Long
TableData<T,N,ORDER>::size () const noexcept
{
    Long r = 1;
    for (int i = 0; i < N; ++i) {
        r *= m_hi[i] - m_lo[i] + 1;
    }
    return r;
}

template <typename T, int N, typename ORDER>
void
TableData<T,N,ORDER>::clear () noexcept
{
    if (m_dptr) {
        if (m_ptr_owner) {
            this->free(m_dptr);
        }
        m_dptr = nullptr;
        m_truesize = 0;
    }
}

template <typename T, int N, typename ORDER>
void
TableData<T,N,ORDER>::define ()
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
    template <typename T, typename>
    Table1D<T> make_table (T* p, Array<int,1> const& lo, Array<int,1> const& hi) {
        return Table1D<T>(p, lo[0], hi[0]+1);
    }
    template <typename T, typename ORDER>
    Table2D<T,ORDER> make_table (T* p, Array<int,2> const& lo, Array<int,2> const& hi) {
        return Table2D<T,ORDER>(p, {lo[0],lo[1]}, {hi[0]+1,hi[1]+1});
    }
    template <typename T, typename ORDER>
    Table3D<T> make_table (T* p, Array<int,3> const& lo, Array<int,3> const& hi) {
        return Table3D<T,ORDER>(p, {lo[0],lo[1],lo[2]}, {hi[0]+1,hi[1]+1,hi[2]+1});
    }
    template <typename T, typename ORDER>
    Table4D<T> make_table (T* p, Array<int,4> const& lo, Array<int,4> const& hi) {
        return Table4D<T,ORDER>(p, {lo[0],lo[1],lo[2],lo[3]}, {hi[0]+1,hi[1]+1,hi[2]+1,hi[3]+1});
    }
}

template <typename T, int N, typename ORDER>
typename TableData<T,N,ORDER>::table_type
TableData<T,N,ORDER>::table () noexcept
{
    return detail::make_table<T,ORDER>(m_dptr, m_lo, m_hi);
}

template <typename T, int N, typename ORDER>
typename TableData<T,N,ORDER>::const_table_type
TableData<T,N,ORDER>::table () const noexcept
{
    return detail::make_table<T const, ORDER>(m_dptr, m_lo, m_hi);
}

template <typename T, int N, typename ORDER>
typename TableData<T,N,ORDER>::const_table_type
TableData<T,N,ORDER>::const_table () const noexcept
{
    return detail::make_table<T const, ORDER>(m_dptr, m_lo, m_hi);
}

template <typename T, int N, typename ORDER>
void
TableData<T,N,ORDER>::copy (TableData<T,N,ORDER> const& rhs) noexcept
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
