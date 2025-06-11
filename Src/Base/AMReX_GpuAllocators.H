#ifndef AMREX_GPUALLOCATORS_H_
#define AMREX_GPUALLOCATORS_H_
#include <AMReX_Config.H>

#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_GpuDevice.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#endif // AMREX_USE_CUDA

#include <limits>
#include <map>
#include <memory>
#include <type_traits>

namespace amrex {

    template <typename T>
    struct FatPtr
    {
        T* m_ptr = nullptr;
        std::size_t m_size = 0;
        [[nodiscard]] constexpr T* ptr () const noexcept { return m_ptr; }
        [[nodiscard]] constexpr std::size_t size () const noexcept { return m_size; }
    };

    template <class T, class AR>
    struct ArenaAllocatorBase
    {
        using value_type = T;
        using arena_wrapper_type = AR;

        constexpr ArenaAllocatorBase () = default;
        explicit constexpr ArenaAllocatorBase (AR a_ar) : m_ar(a_ar) {}

        [[nodiscard]] T* allocate (std::size_t n)
        {
            return (T*) arena()->alloc(n * sizeof(T));
        }

        [[nodiscard]] FatPtr<T>
        allocate_in_place (T* p, std::size_t nmin, std::size_t nmax)
        {
            auto pn = arena()->alloc_in_place(p, nmin*sizeof(T), nmax*sizeof(T));
            return FatPtr<T>{(T*)pn.first, pn.second/sizeof(T)};
        }

        [[nodiscard]] T*
        shrink_in_place (T* p, std::size_t n)
        {
            return (T*) arena()->shrink_in_place(p,n*sizeof(T));
        }

        void deallocate (T* ptr, std::size_t)
        {
            if (ptr != nullptr) { arena()->free(ptr); }
        }

        [[nodiscard]] Arena* arena () const noexcept {
            return m_ar.arena();
        }

    private:
        AR m_ar{};
    };

    struct ArenaWrapper {
        [[nodiscard]] static Arena* arena () noexcept {
            return The_Arena();
        }
    };

    struct DeviceArenaWrapper {
        [[nodiscard]] static Arena* arena () noexcept {
            return The_Device_Arena();
        }
    };

    struct PinnedArenaWrapper {
        [[nodiscard]] static Arena* arena () noexcept {
            return The_Pinned_Arena();
        }
    };

    struct ManagedArenaWrapper {
        [[nodiscard]] static Arena* arena () noexcept {
            return The_Managed_Arena();
        }
    };

    struct AsyncArenaWrapper {
        [[nodiscard]] static Arena* arena () noexcept {
            return The_Async_Arena();
        }
    };

    struct PolymorphicArenaWrapper {
        constexpr PolymorphicArenaWrapper () = default;
        explicit constexpr PolymorphicArenaWrapper (Arena* a_arena)
            : m_arena(a_arena) {}
        [[nodiscard]] Arena* arena () const noexcept {
            return (m_arena) ? m_arena : The_Arena();
        }
        Arena* m_arena = nullptr;
    };

    template<typename T>
    class ArenaAllocator
        : public ArenaAllocatorBase<T,ArenaWrapper>
    {
    };

    template<typename T>
    class DeviceArenaAllocator
        : public ArenaAllocatorBase<T,DeviceArenaWrapper>
    {
    };

    template<typename T>
    class PinnedArenaAllocator
        : public ArenaAllocatorBase<T,PinnedArenaWrapper>
    {
    };

    template<typename T>
    class ManagedArenaAllocator
        : public ArenaAllocatorBase<T,ManagedArenaWrapper>
    {
    };

    template<typename T>
    class AsyncArenaAllocator
        : public ArenaAllocatorBase<T,AsyncArenaWrapper>
    {
    };

    template<typename T>
    class PolymorphicArenaAllocator
        : public ArenaAllocatorBase<T,PolymorphicArenaWrapper>
    {
    public :
        constexpr PolymorphicArenaAllocator () = default;
        explicit constexpr PolymorphicArenaAllocator (Arena* a_arena)
            : ArenaAllocatorBase<T,PolymorphicArenaWrapper>
                (PolymorphicArenaWrapper(a_arena))
            {}
        void setArena (Arena* a_ar) noexcept
        {
            *this = PolymorphicArenaAllocator<T>(a_ar);
        }
    };

    template <typename T>
    struct RunOnGpu : std::false_type {};

    template <class T, class Enable = void>
    struct IsArenaAllocator : std::false_type {};
    //
    template <class T>
    struct IsArenaAllocator
                <T,std::enable_if_t<std::is_base_of_v
                                    <ArenaAllocatorBase<typename T::value_type,
                                                        typename T::arena_wrapper_type>,
                                     T>>>
        : std::true_type {};

    template <typename T>
    struct IsPolymorphicArenaAllocator : std::false_type {};

#ifdef AMREX_USE_GPU
    template <typename T>
    struct RunOnGpu<ArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<DeviceArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<ManagedArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<AsyncArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct IsPolymorphicArenaAllocator<PolymorphicArenaAllocator<T> > : std::true_type {};

#endif // AMREX_USE_GPU

#ifdef AMREX_USE_GPU
    template <class T>
    using DefaultAllocator = amrex::ArenaAllocator<T>;
#else
    template <class T>
    using DefaultAllocator = std::allocator<T>;
#endif // AMREX_USE_GPU

    template <typename A1, typename A2,
              std::enable_if_t<IsArenaAllocator<A1>::value &&
                               IsArenaAllocator<A2>::value, int> = 0>
    bool operator== (A1 const& a1, A2 const& a2)
    {
        return a1.arena() == a2.arena();
    }

    template <typename A1, typename A2,
              std::enable_if_t<IsArenaAllocator<A1>::value &&
                               IsArenaAllocator<A2>::value, int> = 0>
    bool operator!= (A1 const& a1, A2 const& a2)
    {
        return a1.arena() != a2.arena();
    }

} // namespace amrex

#endif // AMREX_GPUALLOCATORS_H_
