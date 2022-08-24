#include "AMReX_InSituUtils.H"

#include "Error.h"
#include "SVTKUtils.h"

namespace amrex {
namespace InSituUtils {

// --------------------------------------------------------------------------
int StateMap::GetIndex(const std::string &name, int centering,
    int &desc, int &comp)
{
    auto cit = this->Map.find(centering);

    if (cit == this->Map.end())
    {
        SENSEI_ERROR("No " << sensei::SVTKUtils::GetAttributesName(centering)
          << " arrays")
        return -1;
    }

    auto nit = cit->second.find(name);
    if (nit == cit->second.end())
    {
        SENSEI_ERROR("No array named \"" << name  << "\" in "
            << sensei::SVTKUtils::GetAttributesName(centering)
            << " centered data")
        return -1;
    }

    desc = nit->second.first;
    comp = nit->second.second;

    return 0;
}

// --------------------------------------------------------------------------
int StateMap::GetName(int centering, int id, std::string &name)
{
    auto cit = this->Map.find(centering);

    if (cit == this->Map.end())
    {
        SENSEI_ERROR("No " << sensei::SVTKUtils::GetAttributesName(centering)
          << " arrays")
        return -1;
    }

    if (size_t(id) >= cit->second.size())
    {
        SENSEI_ERROR("Array index " << id << " out of bounds " << cit->second.size())
        return -1;
    }

    auto nit = cit->second.begin();
    while (id--)
        ++nit;

    name = nit->first;

    return 0;
}

}
}
