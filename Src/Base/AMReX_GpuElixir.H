#ifndef AMREX_GPUELIXIR_H_
#define AMREX_GPUELIXIR_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuElixir.H
 * \brief RAII handle for deferring GPU memory deallocation until stream completion.
 */

#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <utility>

namespace amrex {
namespace Gpu {

/**
 * \class amrex::Elixir
 * \ingroup amrex_memory
 * \brief RAII handle that frees GPU memory when it goes out of scope.
 *
 * An Elixir holds a set of (pointer, arena) pairs. When the Elixir is
 * destroyed it calls arena->free(pointer) for each pair. This is used to
 * keep GPU memory alive until a kernel or stream operation that uses it
 * has completed, at which point the Elixir is allowed to go out of scope.
 *
 * Elixirs are move-only; they cannot be copied.
 */
class Elixir
{
public:

    //! Construct an empty Elixir holding no memory.
    Elixir () noexcept = default;

    //! Construct an Elixir that will free \p p from \p arena when destroyed.
    Elixir (void* p, Arena* arena) noexcept : m_pa({std::make_pair(p,arena)}) {}

    Elixir (Elixir const&) = delete;
    Elixir& operator= (Elixir const&) = delete;

    //! Move constructor: transfer ownership of all held pointers.
    Elixir (Elixir && rhs) noexcept
        : m_pa(std::move(rhs.m_pa))
    {
        rhs.m_pa.clear();
    }

    //! Move assignment: release current pointers then take ownership from \p rhs.
    Elixir& operator= (Elixir && rhs) noexcept
    {
        clear();
        m_pa = std::move(rhs.m_pa);
        rhs.m_pa.clear();
        return *this;
    }

    //! Transfer all (pointer, arena) pairs from \p rhs into this Elixir.
    void append (Elixir rhs) noexcept
    {
        m_pa.insert(m_pa.end(), rhs.m_pa.begin(), rhs.m_pa.end());
        rhs.m_pa.clear();
    }

    ~Elixir () { clear(); }

    //! Free all held pointers via their respective arenas.
    void clear () noexcept;

private:
    Vector<std::pair<void*,Arena*> > m_pa;
};

}

/// \ingroup amrex_memory
using Gpu::Elixir;
}

#endif
