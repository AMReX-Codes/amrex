/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef SMART_UTILS_H_
#define SMART_UTILS_H_

#include "DefaultInitialization.H"

using NameMap = std::map<std::string, int>;
using PolicyVec = amrex::Gpu::DeviceVector<InitializationPolicy>;

struct SmartCopyTag
{
    std::vector<std::string> common_names;
    amrex::Gpu::DeviceVector<int> src_comps;
    amrex::Gpu::DeviceVector<int> dst_comps;

    int size () const noexcept { return common_names.size(); }
};

PolicyVec getPolicies (const NameMap& names) noexcept;

SmartCopyTag getSmartCopyTag (const NameMap& src, const NameMap& dst) noexcept;

/**
 * \brief Sets the ids of newly created particles to the next values.
 *
 * \tparam PTile the particle tile type
 *
 * \param ptile the particle tile
 * \param old_size the index of the first new particle
 * \param num_added the number of particles to set the ids for.
 */
template <typename PTile>
void setNewParticleIDs (PTile& ptile, int old_size, int num_added)
{
    int pid;
#ifdef _OPENMP
#pragma omp critical (ionization_nextid)
#endif
    {
        pid = PTile::ParticleType::NextID();
        PTile::ParticleType::NextID(pid + num_added);
    }

    const int cpuid = amrex::ParallelDescriptor::MyProc();
    auto pp = ptile.GetArrayOfStructs()().data() + old_size;
    amrex::ParallelFor(num_added, [=] AMREX_GPU_DEVICE (int ip) noexcept
    {
        auto& p = pp[ip];
        p.id() = pid+ip;
        p.cpu() = cpuid;
    });
}

#endif //SMART_UTILS_H_
