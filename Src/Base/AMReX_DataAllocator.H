#ifndef AMREX_DATA_ALLOCATOR_H_
#define AMREX_DATA_ALLOCATOR_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_DataAllocator.H
 * \brief Lightweight Arena-backed allocator and deleter helpers.
 */

#include <AMReX_Arena.H>

namespace amrex {

/**
 * \struct amrex::DataAllocator
 * \ingroup amrex_memory
 * \brief Lightweight allocator/deallocator backed by an Arena.
 *
 * Falls back to The_Arena() when no arena is supplied.
 */
struct DataAllocator {
    Arena* m_arena = nullptr; ///< Stored arena pointer (may be null).

    //! Default-construct; alloc/free will use The_Arena().
    DataAllocator ()  noexcept  = default;

    //! Construct with an explicit arena \p ar.
    DataAllocator (Arena* ar) noexcept : m_arena(ar) {}

    //! Allocate \p sz bytes from the arena.
    [[nodiscard]] void* alloc (std::size_t sz) const noexcept {
        return arena()->alloc(sz);
    }

    //! Return the allocation at \p pt to the arena.
    void free (void* pt) const noexcept {
        arena()->free(pt);
    }

    //! Return the stored arena, or The_Arena() if none was set.
    [[nodiscard]] Arena* arena () const noexcept {
        return (m_arena) ? m_arena : The_Arena();
    }
};

/**
 * \struct amrex::DataDeleter
 * \ingroup amrex_memory
 * \brief Custom deleter for \c std::unique_ptr that returns memory to an Arena.
 *
 * Falls back to The_Arena() when no arena is supplied.
 */
struct DataDeleter {
    Arena* m_arena = nullptr; ///< Stored arena pointer (may be null).
    //! Default-construct; operator() will use The_Arena().
    DataDeleter () noexcept  = default;
    //! Construct with an explicit arena \p ar.
    DataDeleter (Arena* ar) noexcept : m_arena(ar) {}
    //! Return the stored arena, or The_Arena() if none was set.
    [[nodiscard]] Arena* arena () const noexcept { return (m_arena) ? m_arena : The_Arena(); }
    //! Return the allocation at \p pt to the arena.
    void operator() (void* pt) const noexcept { arena()->free(pt); }
};

}

#endif
