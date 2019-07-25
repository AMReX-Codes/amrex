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

std::size_t
InjectorMomentum::sharedMemoryNeeded () const noexcept
{
    switch (type)
    {
    case Type::parser:
    {
        return amrex::Gpu::numThreadsPerBlockParallelFor() * sizeof(double);
    }
    default:
        return 0;
    }
}

bool
InjectorMomentum::useRandom () const noexcept
{
    switch (type)
    {
    case Type::gaussian:
    {
        return true;
    }
    default:
        return false;
    }
}
