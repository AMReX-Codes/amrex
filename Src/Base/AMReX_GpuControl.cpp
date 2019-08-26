
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>

namespace amrex {
namespace Gpu {

#if defined(AMREX_USE_GPU)
bool in_launch_region = true;
bool in_graph_region = false;

ScopedDefaultStream::ScopedDefaultStream () noexcept
    : m_prev_stream(Device::resetStream())
{}

ScopedDefaultStream::~ScopedDefaultStream ()
{
    Device::setStream(m_prev_stream);
}

#endif

}
}
