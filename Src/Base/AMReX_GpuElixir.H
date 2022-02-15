#ifndef AMREX_GPUELIXIR_H_
#define AMREX_GPUELIXIR_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <utility>

namespace amrex {
namespace Gpu {

class Elixir
{
public:

    Elixir () noexcept {}

    Elixir (void* p, Arena* arena) noexcept : m_pa({std::make_pair(p,arena)}) {}

    Elixir (Elixir const&) = delete;
    void operator= (Elixir const&) = delete;

    Elixir (Elixir && rhs) noexcept
        : m_pa(std::move(rhs.m_pa))
    {
        rhs.m_pa.clear();
    }

    void operator= (Elixir && rhs) noexcept
    {
        clear();
        m_pa = std::move(rhs.m_pa);
        rhs.m_pa.clear();
    }

    void append (Elixir && rhs) noexcept
    {
        m_pa.insert(m_pa.end(), rhs.m_pa.begin(), rhs.m_pa.end());
        rhs.m_pa.clear();
    }

    ~Elixir () { clear(); }

    void clear () noexcept;

private:
    Vector<std::pair<void*,Arena*> > m_pa;
};

}

using Gpu::Elixir;
}

#endif
