#ifndef AMREX_CUDAELIXIR_H_
#define AMREX_CUDAELIXIR_H_

#include <AMReX_Arena.H>

namespace amrex {
namespace Cuda {

class Elixir
{
public:

    Elixir (void* p, Arena* arena) : m_p(p), m_arena(arena) {}

    Elixir (Elixir const&) = delete;
    void operator= (Elixir const&) = delete;
    void operator= (Elixir &&) = delete;

    Elixir (Elixir && rhs)
        : m_p(rhs.m_p), m_arena(rhs.m_arena)
    {
        rhs.m_p = nullptr;
        rhs.m_arena = nullptr;
    }

    ~Elixir () { clear(); }

    void clear ();

private:
    void* m_p;
    Arena* m_arena;
};

}

using Cuda::Elixir;
}

#endif
