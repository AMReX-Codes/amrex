#include <PlasmaInjector.H>

using namespace amrex;

bool
InjectorPosition::useRandom () const noexcept
{
    switch (type)
    {
    case Type::random:
    {
        return true;
    }
    default:
    {
        return false;
    }
    };
}
