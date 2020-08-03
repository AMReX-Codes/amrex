#ifndef AMREX_PODVECTOR_H_
#define AMREX_PODVECTOR_H_

#include <AMReX_TypeTraits.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>

#include <type_traits>
#include <utility>
#include <memory>
#include <cstring>

namespace amrex
{
    namespace detail
    {
        template <typename T, typename U, typename Size, typename Value>
        typename std::enable_if<RunOnGpu<T>::value>::type
        uninitializedFillNImpl (U* data, Size count, const Value& value)
        {
            amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (Size i) noexcept {
                data[i] = value;
            });
            Gpu::Device::streamSynchronize();
        }

        template <typename T, typename U, typename Size, typename Value>
        typename std::enable_if<!RunOnGpu<T>::value>::type
        uninitializedFillNImpl (U* data, Size count, const Value& value)
        {
            std::uninitialized_fill_n<U*, U>(data, count, value);
        }

        template <typename T>
        typename std::enable_if<RunOnGpu<T>::value>::type
        memCopyImpl (void* dst, const void* src, std::size_t count)
        {
#ifdef AMREX_USE_GPU
            Gpu::dtod_memcpy(dst, src, count);
#else
            std::memcpy(dst, src, count);
#endif
        }

        template <typename T>
        typename std::enable_if<!RunOnGpu<T>::value>::type
        memCopyImpl (void* dst, const void* src, std::size_t count)
        {
            std::memcpy(dst, src, count);
        }

        template <typename T>
        typename std::enable_if<RunOnGpu<T>::value>::type
        memMoveImpl (void* dst, const void* src, std::size_t count, T& allocator)
        {
#ifdef AMREX_USE_GPU
            AMREX_ASSERT(count % sizeof(typename T::value_type) == 0);

            auto N = count / sizeof(typename T::value_type);
            auto tmp = allocator.allocate(N);

            Gpu::dtod_memcpy(tmp, src, count);
            Gpu::dtod_memcpy(dst, tmp, count);

            allocator.deallocate(tmp, N);
#else
            amrex::ignore_unused(allocator);
            std::memmove(dst, src, count);
#endif
        }

        template <typename T>
        typename std::enable_if<!RunOnGpu<T>::value>::type
        memMoveImpl (void* dst, const void* src, std::size_t count, T& /*allocator*/)
        {
            std::memmove(dst, src, count);
        }

        template <typename T, typename U>
        typename std::enable_if<RunOnGpu<T>::value>::type
        initFromListImpl (U* data, std::initializer_list<U> list)
        {
#ifdef AMREX_USE_GPU
            Gpu::htod_memcpy(data, &(*list.begin()), list.size() * sizeof(U));
#else
            std::memcpy(data, &(*list.begin()), list.size() * sizeof(U));
#endif
        }

        template <typename T, typename U>
        typename std::enable_if<!RunOnGpu<T>::value>::type
        initFromListImpl (U* data, std::initializer_list<U> list)
        {
            for (const auto& value : list) { *(data++) = value; }
        }
    }
    
    template <class T, class Allocator = std::allocator<T> >
    class PODVector : public Allocator
    {			
        static_assert(std::is_standard_layout<T>(), "PODVector can only hold standard layout types");
        static_assert(AMREX_IS_TRIVIALLY_COPYABLE(T), "PODVector can only hold trivially copyable types");
            
        using Allocator::allocate;
        using Allocator::deallocate;
            
    public:            
        typedef T value_type;
        typedef Allocator allocator_type;
        typedef std::size_t size_t;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        
        typedef T& reference;
        typedef T* pointer;
        typedef T* iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;
        
        typedef const T& const_reference;
        typedef const T* const_pointer;
        typedef const T* const_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        
    private:        
        pointer m_data;
        size_type m_size, m_capacity;

    public:
        explicit PODVector () noexcept :
        Allocator(), m_data(nullptr), m_size(0), m_capacity(0)
        {}
	
        explicit PODVector (const allocator_type& a_allocator) noexcept
            : Allocator(a_allocator), m_data(nullptr), m_size(0), m_capacity(0)
        {}
        
        explicit PODVector (size_type a_size) noexcept
            : m_data(nullptr), m_size(0), m_capacity(0) 
        {
            AllocateBuffer(GetNewCapacity(a_size));
            m_size = a_size;
        }
	
        PODVector (size_type a_size, const value_type& a_value,
                   const allocator_type& a_allocator = Allocator()) noexcept
            : Allocator(a_allocator), m_data(nullptr), m_size(0), m_capacity(0) 
        {
            AllocateBuffer(GetNewCapacity(a_size));
            m_size = a_size;
            detail::uninitializedFillNImpl<Allocator>(m_data, a_size, a_value);
        }
	
        PODVector (std::initializer_list<T> a_initializer_list,
                   const allocator_type& a_allocator = Allocator()) noexcept
            : Allocator(a_allocator), m_data(nullptr), m_size(0), m_capacity(0) 
        {
            size_type size = a_initializer_list.size();
            AllocateBuffer(GetNewCapacity(size));
            m_size = size;            
            detail::initFromListImpl<Allocator>(m_data, a_initializer_list);
        }
	
        PODVector (const PODVector<T, Allocator>& a_vector) noexcept
            : Allocator(a_vector), m_data(nullptr), m_size(0), m_capacity(0)
        {
            using namespace detail;
            AllocateBuffer(a_vector.capacity());
            m_size = a_vector.size();
            memCopyImpl<Allocator>(m_data, a_vector.m_data, a_vector.size() * sizeof(T));
        }
        
        PODVector (PODVector<T, Allocator>&& a_vector) noexcept
            : Allocator(std::move(a_vector)),
              m_data(a_vector.m_data),
              m_size(a_vector.m_size),
              m_capacity(a_vector.m_capacity)
        {
            a_vector.m_data = nullptr;
            a_vector.m_size = 0;
            a_vector.m_capacity = 0;
        }
	
        PODVector (PODVector<T, Allocator>&& a_vector, const allocator_type& a_allocator) noexcept
            : Allocator(a_allocator),
              m_data(a_vector.m_data),
              m_size(a_vector.m_size),
              m_capacity(a_vector.m_capacity)
        {
            a_vector.m_data = nullptr;
            a_vector.m_size = 0;
            a_vector.m_capacity = 0;
        }
	
        ~PODVector () noexcept { deallocate(m_data, capacity()); }
	
        PODVector& operator= (const PODVector<T, Allocator>& a_vector) noexcept
        {
            return CopyAssignmentDoIt(a_vector, typename std::allocator_traits<Allocator>::
                                      propagate_on_container_copy_assignment());
        }
	
        PODVector& operator= (PODVector<T, Allocator>&& a_vector) noexcept
            (std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value)
        {
            return MoveAssignmentDoIt(std::move(a_vector), typename std::allocator_traits<Allocator>::
                                      propagate_on_container_move_assignment());
        }
        
        iterator erase (const_iterator a_pos) noexcept
        {
            --m_size;
            detail::memMoveImpl<Allocator>(const_cast<iterator>(a_pos), a_pos+1, (end() - a_pos)*sizeof(T), *this);
            return const_cast<iterator>(a_pos);
        }
        
        iterator erase (const_iterator a_first, const_iterator a_last) noexcept
        {
            size_type num_to_erase = a_last - a_first;
            if (num_to_erase == 0) return const_cast<iterator>(a_first);
            m_size -= num_to_erase;
            detail::memMoveImpl<Allocator>(const_cast<iterator>(a_first), a_first+num_to_erase, (end() - a_first)*sizeof(T), *this);
            return const_cast<iterator>(a_first);
        }
        
        iterator insert (const_iterator a_pos, const T& a_item) noexcept
        {
            if(m_size == m_capacity)
            {
                size_t insert_index = std::distance(m_data, const_cast<iterator>(a_pos));
                AllocateBufferForInsert(GetNewCapacity(1), insert_index, 1);
                a_pos = m_data;
                std::advance(a_pos, insert_index);
            }
            else
            {
                detail::memMoveImpl<Allocator>(const_cast<iterator>(a_pos)+1, a_pos, (end() - a_pos) * sizeof(T), *this);
                ++m_size;
            }
            
            *const_cast<iterator>(a_pos) = a_item;
            
            return const_cast<iterator>(a_pos);
        }
        
        iterator insert (const_iterator a_pos, size_type a_count, const T& a_value) noexcept
        {
            if (a_count == 0) return const_cast<iterator>(a_pos);
            if( capacity() < size() + a_count)
            {
                size_t insert_index = std::distance(m_data, const_cast<iterator>(a_pos));
                AllocateBufferForInsert(GetNewCapacity(a_count), insert_index, a_count);
                a_pos = m_data;
                std::advance(a_pos, insert_index);
            }
            else
            {
                detail::memMoveImpl<Allocator>(const_cast<iterator>(a_pos)+a_count, a_pos, (end() - a_pos) * sizeof(T), *this);
                m_size += a_count;
            }

            detail::uninitializedFillNImpl<Allocator>(const_cast<iterator>(a_pos), a_count, a_value);
            
            return const_cast<iterator>(a_pos);
        }
        
        iterator insert (const_iterator a_pos, T&& a_item) noexcept
        {
            if(m_size == m_capacity)
            {
                size_t insert_index = std::distance(m_data, const_cast<iterator>(a_pos));
                AllocateBufferForInsert(GetNewCapacity(1), insert_index, 1);
                a_pos = m_data;
                std::advance(a_pos, insert_index);
            }
            else
            {
                detail::memMoveImpl<Allocator>(const_cast<iterator>(a_pos)+1, a_pos, (end() - a_pos) * sizeof(T), *this);
                ++m_size;
            }
            
            *const_cast<iterator>(a_pos) = a_item;
            
            return const_cast<iterator>(a_pos);
        }
        
        iterator insert (const_iterator a_pos, std::initializer_list<T> a_initializer_list) noexcept
        {
            size_type count = a_initializer_list.size();
            if( capacity() < size() + count)
            {
                size_t insert_index = std::distance(m_data, const_cast<iterator>(a_pos));
                AllocateBufferForInsert(GetNewCapacity(count), insert_index, count);
                a_pos = m_data;
                std::advance(a_pos, insert_index);
            }
            else
            {
                detail::memMoveImpl<Allocator>(const_cast<iterator>(a_pos)+count, a_pos, (end() - a_pos) * sizeof(T), *this);
                m_size += count;
            }        // copy the initializer list
            pointer dst = const_cast<iterator>(a_pos);
            for (const auto& value : a_initializer_list) { *(dst++) = value; }        
            return const_cast<iterator>(a_pos);
        }    
        
        template <class InputIt, class bar = typename std::iterator_traits<InputIt>::difference_type>
        iterator insert (const_iterator a_pos, InputIt a_first, InputIt a_last) noexcept
        {
            size_type count = std::distance(a_first, a_last);
            if (count == 0) return const_cast<iterator>(a_pos);
            if( capacity() < size() + count)
            {
                size_t insert_index = std::distance(m_data, const_cast<iterator>(a_pos));
                AllocateBufferForInsert(GetNewCapacity(count), insert_index, count);
                a_pos = m_data;
                std::advance(a_pos, insert_index);
            }
            else
            {
                detail::memMoveImpl<Allocator>(const_cast<iterator>(a_pos)+count, a_pos, (end() - a_pos) * sizeof(T), *this);
                m_size += count;
            }        pointer dst = const_cast<iterator>(a_pos);
            while(a_first != a_last) { *(dst++) = *(a_first++); }
            return const_cast<iterator>(a_pos);
        }
        
        void assign (size_type a_count, const T& a_value) noexcept
        {
            if ( a_count > capacity() ) AllocateBuffer(GetNewCapacity(a_count));
            m_size = a_count;
            detail::uninitializedFillNImpl<Allocator>(m_data, a_count, a_value);
        }    
        
        void assign (std::initializer_list<T> a_initializer_list) noexcept
        {
            if(a_initializer_list.size() > capacity())
                AllocateBuffer(GetNewCapacity(a_initializer_list.size()));
            m_size = a_initializer_list.size();        // copy the initializer list
            pointer dst = const_cast<iterator>(m_data);
            for (const auto& value : a_initializer_list) { *(dst++) = value; }
        }    
        
        template <class InputIt, class bar = typename std::iterator_traits<InputIt>::difference_type>
        void assign (InputIt a_first, InputIt a_last) noexcept
        {
            size_t count = std::distance(a_first, a_last);
            if (count > capacity()) AllocateBuffer(GetNewCapacity(count));
            m_size = count;        pointer dst = const_cast<iterator>(m_data);
            while(a_first != a_last) { *(dst++) = *(a_first++); }
        }
        
        // don't have the emplace methods, but not sure how often we use those.	
        allocator_type get_allocator () const noexcept { return *this; }
	
        void push_back (const T& a_value) noexcept
        {
            if (m_size == m_capacity) AllocateBuffer(GetNewCapacity(1));
            m_data[m_size] = a_value;
            ++m_size;
        }
        
        void push_back (T&& a_value) noexcept
        {
            if (m_size == m_capacity) AllocateBuffer(GetNewCapacity(1));
            m_data[m_size] = std::move(a_value);
            ++m_size;
        }
	
        void pop_back () noexcept { --m_size; }
        
        void clear () noexcept { m_size = 0; }
	
        size_type size () const noexcept { return m_size; }
        
        size_type max_size () const noexcept { return Allocator::max_size(); }
	
        size_type capacity () const noexcept { return m_capacity; }
        
        bool empty () const noexcept { return m_size == 0; }
	
        T& operator[] (size_type a_index) noexcept { return m_data[a_index]; }
	
        const T& operator[] (size_type a_index) const noexcept { return m_data[a_index]; }
        
        T& front () noexcept { return *m_data; }
	
        const T& front () const noexcept { return *m_data; }
	
        T& back () noexcept { return *(m_data + m_size - 1); }
	
        const T& back () const noexcept { return *(m_data + m_size - 1); }
	
        T* data () noexcept { return m_data; }
	
        const T* data () const noexcept { return m_data; }
        
        T* dataPtr () noexcept { return m_data; }
	
        const T* dataPtr () const noexcept { return m_data; }
        
        iterator begin () noexcept { return m_data; }
	
        const_iterator begin () const noexcept { return m_data; }
	
        iterator end () noexcept { return m_data + m_size; }
	
        const_iterator end () const noexcept { return m_data + m_size; }
	
        reverse_iterator rbegin () noexcept { return reverse_iterator(end()); }
	
        const_reverse_iterator rbegin () const noexcept { return const_reverse_iterator(end()); }
	
        reverse_iterator rend () noexcept { return reverse_iterator(begin()); }
	
        const_reverse_iterator rend () const noexcept { return const_reverse_iterator(begin()); }
	
        const_iterator cbegin () const noexcept { return m_data; }
	
        const_iterator cend () const noexcept { return m_data + m_size; }
	
        const_reverse_iterator crbegin () const noexcept { return const_reverse_iterator(end()); }
	
        const_reverse_iterator crend () const noexcept { return const_reverse_iterator(begin()); }
        
        void resize (size_type a_new_size) noexcept
        {
            if (capacity() < a_new_size) AllocateBuffer(GetNewCapacity(a_new_size - capacity()));
            m_size = a_new_size;
        }
	
        void resize (size_type a_new_size, const T& a_val) noexcept
        {
            size_type old_size = size();
            resize(a_new_size);
            if (old_size < a_new_size)
            {
                detail::uninitializedFillNImpl<Allocator>(m_data + old_size,
                                                          size() - old_size,
                                                          a_val);
            }
        }
        
        void reserve (size_type a_capacity) noexcept 
        {
            if(capacity() < a_capacity) AllocateBuffer(a_capacity);
        }
        
        void shrink_to_fit () noexcept
        {
            const size_type current_size = size();
            if ( current_size == 0 )
            {
                deallocate(m_data, capacity());
                m_data = nullptr;
                m_size = 0;
                m_capacity = 0;
            }
            else if ( current_size < capacity() )
            {
                AllocateBuffer(current_size);
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
        
        // this is where we would change the growth strategy
        size_type GetNewCapacity (size_type a_num_to_be_added) const noexcept
        {
            size_type new_capacity = capacity();
            
            if (capacity() == 0)
            {
                new_capacity = std::max(64 / sizeof(T), size_type(1));
            }
            
            while (new_capacity < (capacity() + a_num_to_be_added))
            {
                new_capacity = (3 * new_capacity + 1)/2;
            }
            
            return new_capacity;
        }
        
        // this is where we would play games with the allocator
        void AllocateBuffer (size_type a_capacity) noexcept
        {
            pointer new_data = allocate(a_capacity);
            if (m_data) detail::memCopyImpl<Allocator>(new_data, m_data, size() * sizeof(T));
            deallocate(m_data, capacity());
            m_data = new_data;
            m_capacity = a_capacity;
        }
        
        // need to have this version too
        void AllocateBufferForInsert(size_type a_capacity, size_type a_index, size_type a_count) noexcept
        {
            using namespace detail;
            pointer new_data = allocate(a_capacity);
            if (m_data)
            {
                memCopyImpl<Allocator>(new_data, m_data, a_index * sizeof(T));
                memCopyImpl<Allocator>(new_data + a_index + a_count, m_data + a_index,
                                       (size() - a_index)*sizeof(T));
            }
            deallocate(m_data, capacity());
            m_data = new_data;
            m_size = size() + a_count;
            m_capacity = a_capacity;
        }
        
        PODVector& CopyAssignmentDoIt (const PODVector<T, Allocator>& a_vector, std::true_type) noexcept
        {
            const size_t other_size = a_vector.size();            
            if ( other_size > capacity() ) { AllocateBuffer(other_size); }            
            m_size = other_size;            
            detail::memCopyImpl<Allocator>(m_data, a_vector.m_data, size() * sizeof(T));
            return *this;
        }
        
        PODVector& CopyAssignmentDoIt (const PODVector<T, Allocator>& a_vector, std::false_type) noexcept
        {
            if (static_cast<Allocator const&>(a_vector) == static_cast<Allocator&>(*this))
            {
                return CopyAssignmentDoIt(a_vector, std::true_type());
            }
            else
            {
                const size_t other_size = a_vector.size();                
                if ( other_size > capacity() ) { AllocateBuffer(other_size); }                
                m_size = other_size;                    
                detail::memCopyImpl<Allocator>(m_data, a_vector.m_data, size() * sizeof(T));
                Allocator::operator=(static_cast<Allocator const&>(a_vector));
                return *this;
            }
        }
        
        PODVector& MoveAssignmentDoIt (PODVector<T, Allocator>&& a_vector, std::false_type) noexcept
        {
            if(static_cast<Allocator&>(a_vector) == static_cast<Allocator&>(*this))
            {
                deallocate(m_data, capacity());
                
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
                CopyAssignmentDoIt(a_vector, std::false_type());
            }
            
            return *this;
        }
        
        PODVector& MoveAssignmentDoIt (PODVector<T, Allocator>&& a_vector, std::true_type) noexcept
        {
            deallocate(m_data, capacity());
            
            Allocator::operator=(std::move(static_cast<Allocator const&>(a_vector)));
            
            m_data = a_vector.m_data;
            m_size = a_vector.m_size;
            m_capacity = a_vector.m_capacity;
            
            a_vector.m_data = nullptr;
            a_vector.m_size = 0;
            a_vector.m_capacity = 0;
            
            return *this;
        }
    };
}

#endif
