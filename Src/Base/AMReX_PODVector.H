#ifndef AMREX_PODVECTOR_H_
#define AMREX_PODVECTOR_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MemPool.H>
#include <AMReX_TypeTraits.H>

#include <iterator>
#include <type_traits>
#include <utility>
#include <memory>
#include <cstring>

namespace amrex
{
    namespace detail
    {
        template <typename T, typename Size, template<class> class Allocator>
        FatPtr<T> allocate_in_place ([[maybe_unused]] T* p, [[maybe_unused]] Size nmin, Size nmax,
                                     Allocator<T>& allocator)
        {
            if constexpr (IsArenaAllocator<Allocator<T>>::value) {
                return allocator.allocate_in_place(p, nmin, nmax);
            } else {
                T* pnew = allocator.allocate(nmax);
                return {pnew, nmax};
            }
        }

        template <typename T, typename Size, template<class> class Allocator>
        T* shrink_in_place ([[maybe_unused]] T* p, Size n, Allocator<T>& allocator)
        {
            if constexpr (IsArenaAllocator<Allocator<T>>::value) {
                return allocator.shrink_in_place(p, n);
            } else {
                return allocator.allocate(n);
            }
        }

        template <typename T, typename Size, template<class> class Allocator>
        void uninitializedFillNImpl (T* data, Size count, const T& value,
                                     [[maybe_unused]] Allocator<T> const& allocator)
        {
#ifdef AMREX_USE_GPU
            if constexpr (RunOnGpu<Allocator<T>>::value)
            {
                amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (Size i) noexcept {
                    data[i] = value;
                });
                Gpu::streamSynchronize();
                return;
            }
            else if constexpr (IsPolymorphicArenaAllocator<Allocator<T>>::value)
            {
                if (allocator.arena()->isManaged() ||
                    allocator.arena()->isDevice())
                {
                    amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (Size i) noexcept
                    {
                        data[i] = value;
                    });
                    Gpu::streamSynchronize();
                    return;
                }
            }
#endif
            std::uninitialized_fill_n(data, count, value);
        }

        template <typename T, template<class> class Allocator>
        void initFromListImpl (T* data, std::initializer_list<T> const& list,
                               [[maybe_unused]] Allocator<T> const & allocator)
        {
            auto count = list.size() * sizeof(T);
#ifdef AMREX_USE_GPU
            if constexpr (RunOnGpu<Allocator<T>>::value)
            {
                Gpu::htod_memcpy_async(data, std::data(list), count);
                Gpu::streamSynchronize();
                return;
            }
            else if constexpr (IsPolymorphicArenaAllocator<Allocator<T>>::value)
            {
                if (allocator.arena()->isManaged() ||
                    allocator.arena()->isDevice())
                {
                    Gpu::htod_memcpy_async(data, std::data(list), count);
                    Gpu::streamSynchronize();
                    return;
                }
            }
#endif
            std::memcpy(data, std::data(list), count);
        }

        template <typename T, typename Size, template<class> class Allocator>
        void fillValuesImpl (T* dst, T const* src, Size count,
                             [[maybe_unused]] Allocator<T> const& allocator)
        {
#ifdef AMREX_USE_GPU
            if constexpr (RunOnGpu<Allocator<T>>::value)
            {
                amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (Size i) noexcept {
                    dst[i] = src[i];
                });
                Gpu::Device::streamSynchronize();
                return;
            }
            else if constexpr (IsPolymorphicArenaAllocator<Allocator<T>>::value)
            {
                if (allocator.arena()->isManaged() ||
                    allocator.arena()->isDevice())
                {
                    amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (Size i) noexcept
                    {
                        dst[i] = src[i];
                    });
                    Gpu::streamSynchronize();
                    return;
                }
            }
#else
            static_assert(RunOnGpu<Allocator<T>>::value == false);
#endif
            if constexpr (! RunOnGpu<Allocator<T>>::value) {
                for (Size i = 0; i < count; ++i) { dst[i] = src[i]; }
            }
        }

        template <typename Allocator>
        void memCopyImpl (void* dst, const void* src, std::size_t count,
                          [[maybe_unused]] Allocator const& dst_allocator,
                          [[maybe_unused]] Allocator const& src_allocator,
                          [[maybe_unused]] bool sync = true)
        {
#ifdef AMREX_USE_GPU
            if constexpr (RunOnGpu<Allocator>::value)
            {
                Gpu::dtod_memcpy_async(dst, src, count);
                if (sync) { Gpu::streamSynchronize(); }
                return;
            }
            else if constexpr (IsPolymorphicArenaAllocator<Allocator>::value)
            {
                bool dst_on_device = dst_allocator.arena()->isManaged() ||
                                     dst_allocator.arena()->isDevice();
                bool src_on_device = src_allocator.arena()->isManaged() ||
                                     src_allocator.arena()->isDevice();
                if (dst_on_device || src_on_device)
                {
                    if (dst_on_device && src_on_device) {
                        Gpu::dtod_memcpy_async(dst, src, count);
                    } else if (dst_on_device) {
                        Gpu::htod_memcpy_async(dst, src, count);
                    } else {
                        Gpu::dtoh_memcpy_async(dst, src, count);
                    }
                    if (sync) { Gpu::streamSynchronize(); }
                    return;
                }
            }
#endif
            std::memcpy(dst, src, count);
        }

        template <typename Allocator>
        void memMoveImpl (void* dst, const void* src, std::size_t count,
                          [[maybe_unused]] Allocator const& allocator)
        {
#ifdef AMREX_USE_GPU
            if constexpr (RunOnGpu<Allocator>::value)
            {
                auto* tmp = The_Arena()->alloc(count);
                Gpu::dtod_memcpy_async(tmp, src, count);
                Gpu::dtod_memcpy_async(dst, tmp, count);
                Gpu::streamSynchronize();
                The_Arena()->free(tmp);
                return;
            }
            else if constexpr (IsPolymorphicArenaAllocator<Allocator>::value)
            {
                if (allocator.arena()->isManaged() ||
                    allocator.arena()->isDevice())
                {
                    auto* tmp = The_Arena()->alloc(count);
                    Gpu::dtod_memcpy_async(tmp, src, count);
                    Gpu::dtod_memcpy_async(dst, tmp, count);
                    Gpu::streamSynchronize();
                    The_Arena()->free(tmp);
                    return;
                }
            }
#endif
            std::memmove(dst, src, count);
        }

        template <typename T, typename Size, template<class> class Allocator>
        void maybe_init_snan (T* data, Size count, Allocator<T> const& allocator)
        {
            amrex::ignore_unused(data, count, allocator);
            if constexpr (std::is_same_v<float,  std::remove_cv_t<T>> ||
                          std::is_same_v<double, std::remove_cv_t<T>>) {
                if (amrex::InitSNaN()) {
#ifdef AMREX_USE_GPU
                   if constexpr (RunOnGpu<Allocator<T>>::value) {
                       amrex::fill_snan<RunOn::Device>(data, count);
                       Gpu::streamSynchronize();
                       return;
                   } else if constexpr (IsPolymorphicArenaAllocator<Allocator<T>>::value) {
                       if (allocator.arena()->isManaged() ||
                           allocator.arena()->isDevice())
                       {
                           amrex::fill_snan<RunOn::Device>(data, count);
                           Gpu::streamSynchronize();
                           return;
                       }
                   }
#endif
                   amrex::fill_snan<RunOn::Host>(data, count);
                }
            }
        }
    }

    namespace VectorGrowthStrategy
    {
        extern AMREX_EXPORT Real growth_factor;
        inline Real GetGrowthFactor () { return growth_factor; }
        inline void SetGrowthFactor (Real a_factor);

        namespace detail
        {
            void ValidateUserInput ();
        }

        void Initialize ();
    }

    template <class T, class Allocator = std::allocator<T> >
    class PODVector : public Allocator
    {
        //        static_assert(std::is_standard_layout<T>(), "PODVector can only hold standard layout types");
        static_assert(std::is_trivially_copyable<T>(), "PODVector can only hold trivially copyable types");
        //        static_assert(std::is_trivially_default_constructible<T>(), "PODVector can only hold trivial dc types");

        using Allocator::allocate;
        using Allocator::deallocate;

    public:
        using value_type      = T;
        using allocator_type  = Allocator;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        using reference        = T&;
        using pointer          = T*;
        using iterator         = T*;
        using reverse_iterator = std::reverse_iterator<iterator>;

        using const_reference        = const T&;
        using const_pointer          = const T*;
        using const_iterator         = const T*;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    private:
        pointer m_data = nullptr;
        size_type m_size{0}, m_capacity{0};

    public:
        constexpr PODVector () noexcept = default;

        constexpr explicit PODVector (const allocator_type& a_allocator) noexcept
            : Allocator(a_allocator)
        {}

        explicit PODVector (size_type a_size)
            : m_size(a_size), m_capacity(a_size)
        {
            if (a_size != 0) {
                m_data = allocate(m_size);
                detail::maybe_init_snan(m_data, m_size, (Allocator const&)(*this));
            }
        }

        PODVector (size_type a_size, const value_type& a_value,
                   const allocator_type& a_allocator = Allocator())
            : Allocator(a_allocator), m_size(a_size), m_capacity(a_size)
        {
            if (a_size != 0) {
                m_data = allocate(m_size);
                detail::uninitializedFillNImpl(m_data, a_size, a_value,
                                               (Allocator const&)(*this));
            }
        }

        PODVector (std::initializer_list<T> a_initializer_list,
                   const allocator_type& a_allocator = Allocator())
            : Allocator(a_allocator),
              m_size    (a_initializer_list.size()),
              m_capacity(a_initializer_list.size())
        {
            if (a_initializer_list.size() != 0) {
                m_data = allocate(m_size);
                detail::initFromListImpl(m_data, a_initializer_list,
                                         (Allocator const&)(*this));
            }
        }

        PODVector (const PODVector<T, Allocator>& a_vector)
            : Allocator(a_vector),
              m_size    (a_vector.size()),
              m_capacity(a_vector.size())
        {
            if (a_vector.size() != 0) {
                m_data = allocate(m_size);
                detail::memCopyImpl(m_data, a_vector.m_data, a_vector.nBytes(),
                                    (Allocator const&)(*this),
                                    (Allocator const&)a_vector);
            }
        }

        PODVector (PODVector<T, Allocator>&& a_vector) noexcept
            : Allocator(static_cast<Allocator&&>(a_vector)),
              m_data(a_vector.m_data),
              m_size(a_vector.m_size),
              m_capacity(a_vector.m_capacity)
        {
            a_vector.m_data = nullptr;
            a_vector.m_size = 0;
            a_vector.m_capacity = 0;
        }

        ~PODVector ()
        {
            // let's not worry about other allocators
            static_assert(std::is_same<Allocator,std::allocator<T>>::value ||
                          IsArenaAllocator<Allocator>::value);
            if (m_data != nullptr) {
                deallocate(m_data, capacity());
            }
        }

        PODVector& operator= (const PODVector<T, Allocator>& a_vector)
        {
            if (this == &a_vector) { return *this; }

            if ((Allocator const&)(*this) != (Allocator const&)a_vector) {
                if (m_data != nullptr) {
                    deallocate(m_data, m_capacity);
                    m_data = nullptr;
                    m_size = 0;
                    m_capacity = 0;
                }
                (Allocator&)(*this) = (Allocator const&)a_vector;
            }

            const auto other_size = a_vector.size();
            if ( other_size > m_capacity ) {
                clear();
                reserve(other_size);
            }

            m_size = other_size;
            if (m_size > 0) {
                detail::memCopyImpl(m_data, a_vector.m_data, nBytes(),
                                    (Allocator const&)(*this),
                                    (Allocator const&)a_vector);
            }
            return *this;
        }

        PODVector& operator= (PODVector<T, Allocator>&& a_vector) noexcept
        {
            if (this == &a_vector) { return *this; }

            if (static_cast<Allocator const&>(a_vector) ==
                static_cast<Allocator const&>(*this))
            {
                if (m_data != nullptr) {
                    deallocate(m_data, m_capacity);
                }

                m_data = a_vector.m_data;
                m_size = a_vector.m_size;
                m_capacity = a_vector.m_capacity;

                a_vector.m_data = nullptr;
                a_vector.m_size = 0;
                a_vector.m_capacity = 0;
            }
            else
            {
                // if the allocators are not the same we give up and copy
                *this = a_vector; // must copy instead of move
            }

            return *this;
        }

        iterator erase (const_iterator a_pos)
        {
            auto* pos = const_cast<iterator>(a_pos);
            --m_size;
            detail::memMoveImpl(pos, a_pos+1, (end() - pos)*sizeof(T),
                                (Allocator const&)(*this));
            return pos;
        }

        iterator erase (const_iterator a_first, const_iterator a_last)
        {
            size_type num_to_erase = a_last - a_first;
            auto* first = const_cast<iterator>(a_first);
            if (num_to_erase > 0) {
                m_size -= num_to_erase;
                detail::memMoveImpl(first, a_last, (end() - first)*sizeof(T),
                                    (Allocator const&)(*this));
            }
            return first;
        }

        iterator insert (const_iterator a_pos, const T& a_item)
        {
            return insert(a_pos, 1, a_item);
        }

        iterator insert (const_iterator a_pos, size_type a_count, const T& a_value)
        {
            auto* pos = const_cast<iterator>(a_pos);
            if (a_count > 0) {
                if (m_capacity < m_size + a_count)
                {
                    std::size_t insert_index = std::distance(m_data, pos);
                    AllocateBufferForInsert(insert_index, a_count);
                    pos = m_data + insert_index;
                }
                else
                {
                    detail::memMoveImpl(pos+a_count, a_pos, (end() - pos) * sizeof(T),
                                        (Allocator const&)(*this));
                    m_size += a_count;
                }
                detail::uninitializedFillNImpl(pos, a_count, a_value,
                                               (Allocator const&)(*this));
            }
            return pos;
        }

        iterator insert (const_iterator a_pos, T&& a_item)
        {
            // This is *POD* vector after all
            return insert(a_pos, 1, std::move(a_item));
        }

        iterator insert (const_iterator a_pos,
                         std::initializer_list<T> a_initializer_list)
        {
            auto* pos = const_cast<iterator>(a_pos);
            size_type count = a_initializer_list.size();
            if (count > 0) {
                if (m_capacity < m_size + count)
                {
                    std::size_t insert_index = std::distance(m_data, pos);
                    AllocateBufferForInsert(insert_index, count);
                    pos = m_data + insert_index;
                }
                else
                {
                    detail::memMoveImpl(pos+count, a_pos, (end() - pos) * sizeof(T),
                                        (Allocator const&)(*this));
                    m_size += count;
                }
                detail::initFromListImpl(pos, a_initializer_list,
                                         (Allocator const&)(*this));
            }
            return pos;
        }

        template <class InputIt, class bar = typename std::iterator_traits<InputIt>::difference_type>
        iterator insert (const_iterator a_pos, InputIt a_first, InputIt a_last)
        {
            auto* pos = const_cast<iterator>(a_pos);
            size_type count = std::distance(a_first, a_last);
            if (count > 0) {
                if (m_capacity < m_size + count)
                {
                    std::size_t insert_index = std::distance(m_data, pos);
                    AllocateBufferForInsert(insert_index, count);
                    pos = m_data + insert_index;
                }
                else
                {
                    detail::memMoveImpl(pos+count, a_pos, (end() - pos) * sizeof(T),
                                        (Allocator const&)(*this));
                    m_size += count;
                }
                // Unfortunately we don't know whether InputIt points
                // GPU or CPU memory. We will assume it's the same as
                // the vector.
                detail::fillValuesImpl(pos, a_first, count,
                                       (Allocator const&)(*this));
            }
            return pos;
        }

        void assign (size_type a_count, const T& a_value)
        {
            if ( a_count > m_capacity ) {
                clear();
                reserve(a_count);
            }
            m_size = a_count;
            detail::uninitializedFillNImpl(m_data, a_count, a_value,
                                           (Allocator const&)(*this));
        }

        void assign (std::initializer_list<T> a_initializer_list)
        {
            if (a_initializer_list.size() > m_capacity) {
                clear();
                reserve(a_initializer_list.size());
            }
            m_size = a_initializer_list.size();
            detail::initFromListImpl(m_data, a_initializer_list,
                                     (Allocator const&)(*this));
        }

        template <class InputIt, class bar = typename std::iterator_traits<InputIt>::difference_type>
        void assign (InputIt a_first, InputIt a_last)
        {
            std::size_t count = std::distance(a_first, a_last);
            if (count > m_capacity) {
                clear();
                reserve(count);
            }
            m_size = count;
            detail::fillValuesImpl(m_data, a_first, count,
                                   (Allocator const&)(*this));
        }

        /** Set the same value to every element of the vector
         *
         * @param a_value the value to assign
         */
        void assign (const T& a_value)
        {
            assign(m_size, a_value);
        }

        [[nodiscard]] allocator_type get_allocator () const noexcept { return *this; }

        void push_back (const T& a_value)
        {
            if (m_size == m_capacity) {
                auto new_capacity = GetNewCapacityForPush();
                AllocateBufferForPush(new_capacity);
            }
            detail::uninitializedFillNImpl(m_data+m_size, 1, a_value,
                                           (Allocator const&)(*this));
            ++m_size;
        }

        // Because T is trivial, there is no need for push_back(T&&)

        // Don't have the emplace methods, but not sure how often we use those.

        void pop_back () noexcept { --m_size; }

        void clear () noexcept { m_size = 0; }

        [[nodiscard]] size_type size () const noexcept { return m_size; }

        [[nodiscard]] size_type capacity () const noexcept { return m_capacity; }

        [[nodiscard]] bool empty () const noexcept { return m_size == 0; }

        [[nodiscard]] T& operator[] (size_type a_index) noexcept { return m_data[a_index]; }

        [[nodiscard]] const T& operator[] (size_type a_index) const noexcept { return m_data[a_index]; }

        [[nodiscard]] T& front () noexcept { return *m_data; }

        [[nodiscard]] const T& front () const noexcept { return *m_data; }

        [[nodiscard]] T& back () noexcept { return *(m_data + m_size - 1); }

        [[nodiscard]] const T& back () const noexcept { return *(m_data + m_size - 1); }

        [[nodiscard]] T* data () noexcept { return m_data; }

        [[nodiscard]] const T* data () const noexcept { return m_data; }

        [[nodiscard]] T* dataPtr () noexcept { return m_data; }

        [[nodiscard]] const T* dataPtr () const noexcept { return m_data; }

        [[nodiscard]] iterator begin () noexcept { return m_data; }

        [[nodiscard]] const_iterator begin () const noexcept { return m_data; }

        [[nodiscard]] iterator end () noexcept { return m_data + m_size; }

        [[nodiscard]] const_iterator end () const noexcept { return m_data + m_size; }

        [[nodiscard]] reverse_iterator rbegin () noexcept { return reverse_iterator(end()); }

        [[nodiscard]] const_reverse_iterator rbegin () const noexcept { return const_reverse_iterator(end()); }

        [[nodiscard]] reverse_iterator rend () noexcept { return reverse_iterator(begin()); }

        [[nodiscard]] const_reverse_iterator rend () const noexcept { return const_reverse_iterator(begin()); }

        [[nodiscard]] const_iterator cbegin () const noexcept { return m_data; }

        [[nodiscard]] const_iterator cend () const noexcept { return m_data + m_size; }

        [[nodiscard]] const_reverse_iterator crbegin () const noexcept { return const_reverse_iterator(end()); }

        [[nodiscard]] const_reverse_iterator crend () const noexcept { return const_reverse_iterator(begin()); }

        void resize (size_type a_new_size)
        {
            auto old_size = m_size;
            resize_without_init_snan(a_new_size);
            if (old_size < a_new_size) {
                detail::maybe_init_snan(m_data + old_size,
                                        m_size - old_size, (Allocator const&)(*this));
            }
        }

        void resize (size_type a_new_size, const T& a_val)
        {
            size_type old_size = m_size;
            resize_without_init_snan(a_new_size);
            if (old_size < a_new_size)
            {
                detail::uninitializedFillNImpl(m_data + old_size,
                                               m_size - old_size, a_val,
                                               (Allocator const&)(*this));
            }
        }

        void reserve (size_type a_capacity)
        {
            if (m_capacity < a_capacity) {
                auto fp = detail::allocate_in_place(m_data, a_capacity, a_capacity,
                                                    (Allocator&)(*this));
                UpdateDataPtr(fp);
            }
        }

        void shrink_to_fit ()
        {
            if (m_data != nullptr) {
                if (m_size == 0) {
                    deallocate(m_data, m_capacity);
                    m_data = nullptr;
                    m_capacity = 0;
                } else if (m_size < m_capacity) {
                    auto* new_data = detail::shrink_in_place(m_data, m_size,
                                                             (Allocator&)(*this));
                    if (new_data != m_data) {
                        detail::memCopyImpl(new_data, m_data, nBytes(),
                                            (Allocator const&)(*this),
                                            (Allocator const&)(*this));
                        deallocate(m_data, m_capacity);
                    }
                    m_capacity = m_size;
                }
            }
        }

        void swap (PODVector<T, Allocator>& a_vector) noexcept
        {
            std::swap(m_data, a_vector.m_data);
            std::swap(m_size, a_vector.m_size);
            std::swap(m_capacity, a_vector.m_capacity);
            std::swap(static_cast<Allocator&>(a_vector), static_cast<Allocator&>(*this));
        }

    private:

        [[nodiscard]] size_type nBytes () const noexcept
        {
            return m_size*sizeof(T);
        }

        // this is where we would change the growth strategy for push_back
        [[nodiscard]] size_type GetNewCapacityForPush () const noexcept
        {
            if (m_capacity == 0) {
                return std::max(64/sizeof(T), size_type(1));
            } else {
                Real const gf = VectorGrowthStrategy::GetGrowthFactor();
                if (amrex::almostEqual(gf, Real(1.5))) {
                    return (m_capacity*3+1)/2;
                } else {
                    return size_type(gf*Real(m_capacity+1));
                }
            }
        }

        void UpdateDataPtr (FatPtr<T> const& fp)
        {
            auto* new_data = fp.ptr();
            auto new_capacity = fp.size();
            if (m_data != nullptr && m_data != new_data) {
                if (m_size > 0) {
                    detail::memCopyImpl(new_data, m_data, nBytes(),
                                        (Allocator const&)(*this),
                                        (Allocator const&)(*this));
                }
                deallocate(m_data, capacity());
            }
            m_data = new_data;
            m_capacity = new_capacity;
        }

        // This is where we play games with the allocator. This function
        // updates m_data and m_capacity, but not m_size.
        void AllocateBufferForPush (size_type target_capacity)
        {
            auto fp = detail::allocate_in_place(m_data, m_size+1, target_capacity,
                                                (Allocator&)(*this));
            UpdateDataPtr(fp);
        }

        // This is where we play games with the allocator and the growth
        // strategy for insert. This function updates m_data, m_size and
        // m_capacity.
        void AllocateBufferForInsert (size_type a_index, size_type a_count)
        {
            size_type new_size = m_size + a_count;
            size_type new_capacity = std::max(new_size, GetNewCapacityForPush());
            auto fp = detail::allocate_in_place(m_data, new_size, new_capacity,
                                                (Allocator&)(*this));
            auto* new_data = fp.ptr();
            new_capacity = fp.size();

            if (m_data != nullptr) {
                if (m_data == new_data) {
                    if (m_size > a_index) {
                        detail::memMoveImpl(m_data+a_index+a_count, m_data+a_index,
                                            (m_size-a_index)*sizeof(T),
                                            (Allocator const&)(*this));
                    }
                } else {
                    if (m_size > 0) {
                        if (a_index > 0) {
                            detail::memCopyImpl(new_data, m_data, a_index*sizeof(T),
                                                (Allocator const&)(*this),
                                                (Allocator const&)(*this), false);
                        }
                        if (m_size > a_index) {
                            detail::memCopyImpl(new_data+a_index+a_count, m_data+a_index,
                                                (m_size-a_index)*sizeof(T),
                                                (Allocator const&)(*this),
                                                (Allocator const&)(*this), false);
                        }
                        Gpu::streamSynchronize();
                    }
                    deallocate(m_data, m_capacity);
                }
            }
            m_data = new_data;
            m_size = new_size;
            m_capacity = new_capacity;
        }

        void resize_without_init_snan (size_type a_new_size)
        {
            if (m_capacity < a_new_size) {
                reserve(a_new_size);
            }
            m_size = a_new_size;
        }
    };
}

#endif
