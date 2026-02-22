#ifndef AMREX_GPUALLOCATORS_H_
#define AMREX_GPUALLOCATORS_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuAllocators.H
 * \brief Arena-backed STL-compatible allocators for GPU and host memory.
 */

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

/**
 * \struct amrex::FatPtr
 * \ingroup amrex_memory
 * \brief A pointer bundled with its allocation size in elements.
 *
 * Returned by allocate_in_place to convey both the (possibly reused)
 * pointer and the actual number of elements that were allocated.
 *
 * \tparam T Pointed-to type.
 */
template <typename T>
struct FatPtr
{
    T* m_ptr = nullptr;         ///< Pointer to the allocation.
    std::size_t m_size = 0;     ///< Number of elements available at \p m_ptr.
    //! Return the pointer.
    [[nodiscard]] constexpr T* ptr () const noexcept { return m_ptr; }
    //! Return the number of elements.
    [[nodiscard]] constexpr std::size_t size () const noexcept { return m_size; }
};

/**
 * \class amrex::ArenaAllocatorBase
 * \ingroup amrex_memory
 * \brief CRTP base providing STL allocator operations backed by an Arena.
 *
 * Implements allocate(), deallocate(), allocate_in_place(), and
 * shrink_in_place() through an Arena obtained from the \p AR wrapper.
 *
 * \tparam T  Value type.
 * \tparam AR Arena-wrapper type that supplies an \c arena() accessor.
 */
template <class T, class AR>
struct ArenaAllocatorBase
{
    using value_type = T;               ///< Element type (STL requirement).
    using arena_wrapper_type = AR;      ///< Arena-wrapper type.

    //! Default-construct; uses the default arena supplied by \p AR.
    constexpr ArenaAllocatorBase () = default;
    //! Construct with an explicit arena wrapper \p a_ar.
    explicit constexpr ArenaAllocatorBase (AR a_ar) : m_ar(a_ar) {}

    //! Allocate storage for \p n objects of type \p T.
    [[nodiscard]] T* allocate (std::size_t n)
    {
        return (T*) arena()->alloc(n * sizeof(T));
    }

    //! Try to grow or shrink the allocation at \p p in-place to between \p nmin and \p nmax elements; returns the resulting pointer/size pair.
    [[nodiscard]] FatPtr<T>
    allocate_in_place (T* p, std::size_t nmin, std::size_t nmax)
    {
        auto pn = arena()->alloc_in_place(p, nmin*sizeof(T), nmax*sizeof(T));
        return FatPtr<T>{(T*)pn.first, pn.second/sizeof(T)};
    }

    //! Shrink the allocation at \p p to exactly \p n elements; returns the new pointer.
    [[nodiscard]] T*
    shrink_in_place (T* p, std::size_t n)
    {
        return (T*) arena()->shrink_in_place(p,n*sizeof(T));
    }

    //! Deallocate the storage pointed to by \p ptr.
    void deallocate (T* ptr, std::size_t)
    {
        if (ptr != nullptr) { arena()->free(ptr); }
    }

    //! Return a pointer to the underlying Arena.
    [[nodiscard]] Arena* arena () const noexcept {
        return m_ar.arena();
    }

private:
    AR m_ar{};
};

/**
 * \struct amrex::ArenaWrapper
 * \ingroup amrex_memory
 * \brief Arena wrapper that returns The_Arena() (the default GPU/CPU arena).
 */
struct ArenaWrapper {
    //! Return the default AMReX arena.
    [[nodiscard]] static Arena* arena () noexcept {
        return The_Arena();
    }
};

/**
 * \struct amrex::DeviceArenaWrapper
 * \ingroup amrex_memory
 * \brief Arena wrapper that returns The_Device_Arena() (device-only memory).
 */
struct DeviceArenaWrapper {
    //! Return the device-only arena.
    [[nodiscard]] static Arena* arena () noexcept {
        return The_Device_Arena();
    }
};

/**
 * \struct amrex::PinnedArenaWrapper
 * \ingroup amrex_memory
 * \brief Arena wrapper that returns The_Pinned_Arena() (page-locked host memory).
 */
struct PinnedArenaWrapper {
    //! Return the pinned (page-locked) host arena.
    [[nodiscard]] static Arena* arena () noexcept {
        return The_Pinned_Arena();
    }
};

/**
 * \struct amrex::ManagedArenaWrapper
 * \ingroup amrex_memory
 * \brief Arena wrapper that returns The_Managed_Arena() (unified/managed GPU memory).
 */
struct ManagedArenaWrapper {
    //! Return the managed (unified) GPU arena.
    [[nodiscard]] static Arena* arena () noexcept {
        return The_Managed_Arena();
    }
};

/**
 * \struct amrex::AsyncArenaWrapper
 * \ingroup amrex_memory
 * \brief Arena wrapper that returns The_Async_Arena() (stream-ordered GPU memory).
 */
struct AsyncArenaWrapper {
    //! Return the asynchronous (stream-ordered) GPU arena.
    [[nodiscard]] static Arena* arena () noexcept {
        return The_Async_Arena();
    }
};

/**
 * \struct amrex::PolymorphicArenaWrapper
 * \ingroup amrex_memory
 * \brief Arena wrapper that holds a runtime-selectable Arena pointer.
 *
 * Falls back to The_Arena() when no arena is specified.
 */
struct PolymorphicArenaWrapper {
    //! Default-construct; arena() will return The_Arena().
    constexpr PolymorphicArenaWrapper () = default;
    //! Construct with an explicit arena pointer \p a_arena.
    explicit constexpr PolymorphicArenaWrapper (Arena* a_arena)
        : m_arena(a_arena) {}
    //! Return the stored arena, or The_Arena() if none was set.
    [[nodiscard]] Arena* arena () const noexcept {
        return (m_arena) ? m_arena : The_Arena();
    }
    Arena* m_arena = nullptr; ///< Stored arena pointer (may be null).
};

/**
 * \class amrex::ArenaAllocator
 * \ingroup amrex_memory
 * \brief STL-compatible allocator backed by The_Arena() (default GPU/CPU arena).
 * \tparam T Element type.
 */
template<typename T>
class ArenaAllocator
    : public ArenaAllocatorBase<T,ArenaWrapper>
{
};

/**
 * \class amrex::DeviceArenaAllocator
 * \ingroup amrex_memory
 * \brief STL-compatible allocator backed by The_Device_Arena() (device-only memory).
 * \tparam T Element type.
 */
template<typename T>
class DeviceArenaAllocator
    : public ArenaAllocatorBase<T,DeviceArenaWrapper>
{
};

/**
 * \class amrex::PinnedArenaAllocator
 * \ingroup amrex_memory
 * \brief STL-compatible allocator backed by The_Pinned_Arena() (page-locked host memory).
 * \tparam T Element type.
 */
template<typename T>
class PinnedArenaAllocator
    : public ArenaAllocatorBase<T,PinnedArenaWrapper>
{
};

/**
 * \class amrex::ManagedArenaAllocator
 * \ingroup amrex_memory
 * \brief STL-compatible allocator backed by The_Managed_Arena() (unified GPU memory).
 * \tparam T Element type.
 */
template<typename T>
class ManagedArenaAllocator
    : public ArenaAllocatorBase<T,ManagedArenaWrapper>
{
};

/**
 * \class amrex::AsyncArenaAllocator
 * \ingroup amrex_memory
 * \brief STL-compatible allocator backed by The_Async_Arena() (stream-ordered GPU memory).
 * \tparam T Element type.
 */
template<typename T>
class AsyncArenaAllocator
    : public ArenaAllocatorBase<T,AsyncArenaWrapper>
{
};

/**
 * \class amrex::PolymorphicArenaAllocator
 * \ingroup amrex_memory
 * \brief STL-compatible allocator backed by a runtime-selectable Arena.
 *
 * Falls back to The_Arena() when no arena is provided at construction.
 *
 * \tparam T Element type.
 */
template<typename T>
class PolymorphicArenaAllocator
    : public ArenaAllocatorBase<T,PolymorphicArenaWrapper>
{
public :
    //! Default-construct; uses The_Arena().
    constexpr PolymorphicArenaAllocator () = default;
    //! Construct with an explicit arena pointer \p a_arena.
    explicit constexpr PolymorphicArenaAllocator (Arena* a_arena)
        : ArenaAllocatorBase<T,PolymorphicArenaWrapper>
            (PolymorphicArenaWrapper(a_arena))
        {}
    //! Switch to a different arena \p a_ar.
    void setArena (Arena* a_ar) noexcept
    {
        *this = PolymorphicArenaAllocator<T>(a_ar);
    }
};

/**
 * \struct amrex::RunOnGpu
 * \ingroup amrex_memory
 * \brief Type trait that is \c std::true_type when \p T is a GPU-capable arena allocator.
 *
 * Specializations for ArenaAllocator, DeviceArenaAllocator, ManagedArenaAllocator,
 * and AsyncArenaAllocator are enabled when AMREX_USE_GPU is defined.
 *
 * \tparam T Allocator type to test.
 */
template <typename T>
struct RunOnGpu : std::false_type {};

/**
 * \struct amrex::IsArenaAllocator
 * \ingroup amrex_memory
 * \brief Type trait that is \c std::true_type when \p T derives from ArenaAllocatorBase.
 * \tparam T Type to test.
 */
template <class T, class Enable = void>
struct IsArenaAllocator : std::false_type {};
/// \cond DOXYGEN_IGNORE
//
template <class T>
struct IsArenaAllocator
            <T,std::enable_if_t<std::is_base_of_v
                                <ArenaAllocatorBase<typename T::value_type,
                                                    typename T::arena_wrapper_type>,
                                 T>>>
    : std::true_type {};
/// \endcond

/**
 * \struct amrex::IsPolymorphicArenaAllocator
 * \ingroup amrex_memory
 * \brief Type trait that is \c std::true_type when \p T is a PolymorphicArenaAllocator.
 * \tparam T Type to test.
 */
template <typename T>
struct IsPolymorphicArenaAllocator : std::false_type {};
/// \cond DOXYGEN_IGNORE
template <typename T>
struct IsPolymorphicArenaAllocator<PolymorphicArenaAllocator<T> > : std::true_type {};
/// \endcond

#ifdef AMREX_USE_GPU
/// \cond DOXYGEN_IGNORE
    template <typename T>
    struct RunOnGpu<ArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<DeviceArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<ManagedArenaAllocator<T> > : std::true_type {};

    template <typename T>
    struct RunOnGpu<AsyncArenaAllocator<T> > : std::true_type {};
/// \endcond
#endif // AMREX_USE_GPU

#ifdef AMREX_USE_GPU
    /// \ingroup amrex_memory
    //! The default AMReX allocator type; ArenaAllocator when GPU is enabled, std::allocator otherwise.
    template <class T>
    using DefaultAllocator = amrex::ArenaAllocator<T>;
/// \cond DOXYGEN_IGNORE
#else
    template <class T>
    using DefaultAllocator = std::allocator<T>;
#endif // AMREX_USE_GPU
/// \endcond

/// \ingroup amrex_memory
//! Return true if two arena allocators refer to the same arena.
template <typename A1, typename A2,
          std::enable_if_t<IsArenaAllocator<A1>::value &&
                           IsArenaAllocator<A2>::value, int> = 0>
bool operator== (A1 const& a1, A2 const& a2)
{
    return a1.arena() == a2.arena();
}

/// \ingroup amrex_memory
//! Return true if two arena allocators refer to different arenas.
template <typename A1, typename A2,
          std::enable_if_t<IsArenaAllocator<A1>::value &&
                           IsArenaAllocator<A2>::value, int> = 0>
bool operator!= (A1 const& a1, A2 const& a2)
{
    return a1.arena() != a2.arena();
}

} // namespace amrex

#endif // AMREX_GPUALLOCATORS_H_
