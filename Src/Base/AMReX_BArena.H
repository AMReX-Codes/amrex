#ifndef AMREX_BARENA_H_
#define AMREX_BARENA_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>

/**
 * \file AMReX_BArena.H
 * \brief CPU heap-backed implementation of the Arena interface.
 */

namespace amrex {
/**
 * \class amrex::BArena
 * \brief Basic Arena that allocates on the host heap with `std::malloc`.
 * \ingroup amrex_memory
 */
class BArena
    :
    public Arena
{
public:
    /**
     * \brief Allocate host memory.
     * \param sz Number of bytes to allocate.
     * \return Pointer to the allocated buffer.
     */
    /// \return Host pointer allocated with std::malloc.
    [[nodiscard]] void* alloc (std::size_t sz) final;
    /**
     * \brief Free memory previously returned by alloc().
     * \param pt Pointer to the allocation to release.
     */
    void free (void* pt) final;

    /** \brief Return true; BArena allocations are host accessible. */
    [[nodiscard]] bool isDeviceAccessible () const final;
    /** \brief Return true; BArena allocations reside on the host. */
    [[nodiscard]] bool isHostAccessible () const final;

    /** \brief Return false; BArena memory is not managed/unified. */
    [[nodiscard]] bool isManaged () const final;
    /** \brief Return false; BArena operates on the host. */
    [[nodiscard]] bool isDevice () const final;
    /** \brief Return false; BArena does not allocate pinned memory. */
    [[nodiscard]] bool isPinned () const final;
};

}

#endif
