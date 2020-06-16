#ifndef AMREX_GPUALLOCATORS_H_
#define AMREX_GPUALLOCATORS_H_

#include <map>
#include <memory>
#include <limits>
#include <memory>

#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_GpuDevice.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#endif // AMREX_USE_CUDA

namespace amrex {

    template <typename T>
    struct RunOnGpu : std::false_type {};

    template<typename T>
    class ArenaAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Arena()->alloc(n * sizeof(T));
	    return result;
        }

        inline void deallocate(value_type* ptr, std::size_t)
        {
            The_Arena()->free(ptr);
        }
    };

    template<typename T>
    class DeviceArenaAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Device_Arena()->alloc(n * sizeof(T));
	    return result;
        }

        inline void deallocate(value_type* ptr, std::size_t)
        {
            The_Device_Arena()->free(ptr);
        }
    };

    template<typename T>
    class PinnedArenaAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Pinned_Arena()->alloc(n * sizeof(T));
	    return result;
        }

        inline void deallocate(value_type* ptr, std::size_t)
        {
            The_Pinned_Arena()->free(ptr);
        }
    };

    template<typename T>
    class ManagedArenaAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Managed_Arena()->alloc(n * sizeof(T));
	    return result;
        }

        inline void deallocate(value_type* ptr, std::size_t)
        {
            The_Managed_Arena()->free(ptr);
        }
    };

    template<typename T>
    class PolymorphicAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
	    value_type* result = nullptr;
            if (ParallelDescriptor::UseGpuAwareMpi())
            {
                result = (value_type*) The_Device_Arena()->alloc(n * sizeof(T));
            }
            else
            {
                result = (value_type*) The_Pinned_Arena()->alloc(n * sizeof(T));
            }
	    return result;
        }

        inline void deallocate(value_type* ptr, std::size_t)
        {
            if (ParallelDescriptor::UseGpuAwareMpi())
            {
                The_Device_Arena()->free(ptr);
            }
            else
            {
                The_Pinned_Arena()->free(ptr);
            }
        }
    };

    template <template <typename> class Allocator, class T, class U>
    bool
    operator==(Allocator<T> const&, Allocator<U> const&) noexcept
    {
        return true;
    }

    template <template <typename> class Allocator, class T, class U>
    bool
    operator!=(Allocator<T> const& x, Allocator<U> const& y) noexcept
    {
        return !(x == y);
    }

#ifdef AMREX_USE_GPU
    template <class T>
    using DefaultAllocator = amrex::ArenaAllocator<T>;

    template <typename T>
    struct RunOnGpu<ArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<DeviceArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<ManagedArenaAllocator<T> > : std::true_type {};

#endif // AMREX_USE_GPU

    template <class T>
    using DefaultAllocator = std::allocator<T>;

} // namespace amrex

#endif // AMREX_GPUALLOCATORS_H_
