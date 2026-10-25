#include <AMReX_GpuElixir.H>
#include <AMReX_GpuDevice.H>

namespace amrex::Gpu {

void
Elixir::clear ()
{
    for (auto const& pa : m_pa) {
        if (pa.first != nullptr) { Gpu::freeAsync(pa.second, pa.first); }
    }
    m_pa.clear();
}

}
