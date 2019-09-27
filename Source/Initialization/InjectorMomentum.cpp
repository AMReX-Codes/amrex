#include <InjectorMomentum.H>
#include <PlasmaInjector.H>

using namespace amrex;

InjectorMomentum::~InjectorMomentum ()
{
    switch (type)
    {
    case Type::parser:
    {
        object.parser.m_ux_parser.clear();
        object.parser.m_uy_parser.clear();
        object.parser.m_uz_parser.clear();
        break;
    }
    case Type::custom:
    {
        object.custom.clear();
        break;
    }
    }
}

// Compute the amount of memory needed in GPU Shared Memory.
std::size_t
InjectorMomentum::sharedMemoryNeeded () const noexcept
{
    switch (type)
    {
    case Type::parser:
    {
        // For parser injector, the 3D position of each particle
        // is stored in shared memory.
        return amrex::Gpu::numThreadsPerBlockParallelFor() * sizeof(double) * 3;
    }
    default:
        return 0;
    }
}

