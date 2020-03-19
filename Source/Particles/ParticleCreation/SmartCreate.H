/* Copyright 2019-2020 Luca Fedeli, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef SMART_CREATE_H_
#define SMART_CREATE_H_

#include "DefaultInitialization.H"

#include <AMReX_AmrCore.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

/**
 * \brief This is a functor for performing a "smart create" that works
 * in both host and device code.
 *
 * A "smart" create does the following:
 * First, it initializes the position, cpuid and id of
 * the particle (in most cases IDs should be initialized using setNewParticleIDs
 * after their creation).
 * Then, it initializes all the other components according to initialization policies.
 *
 * You don't create this directly - use the SmartCreateFactory object below.
 *
 */

struct SmartCreate
{
    const InitializationPolicy* m_policy_real;
    const InitializationPolicy* m_policy_int;
    const int m_weight_index;

    template <typename PartData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (
        PartData& prt, const int i_prt,
        const amrex::Real x = 0.0,
#if (AMREX_SPACEDIM == 3)
        const amrex::Real y = 0.0,
#endif
        const amrex::Real z = 0.0,
        const int cpu = 0,
        const int id = 0) const noexcept
    {
        prt.m_aos[i_prt].pos(0) = x;
#if (AMREX_SPACEDIM == 3)
        prt.m_aos[i_prt].pos(1) = y;
#endif
        prt.m_aos[i_prt].pos(2) = z;

        prt.m_aos[i_prt].cpu() = cpu;
        prt.m_aos[i_prt].id() = id;

         // initialize the real components
         for (int j = 0; j < PartData::NAR; ++j)
             prt.m_rdata[j][i_prt] = initializeRealValue(m_policy_real[j]);
         for (int j = 0; j < prt.m_num_runtime_real; ++j)
             prt.m_runtime_rdata[j][i_prt] = initializeRealValue(m_policy_real[j+PartData::NAR]);

         // initialize the int components
         for (int j = 0; j < PartData::NAI; ++j)
             prt.m_idata[j][i_prt] = initializeIntValue(m_policy_int[j]);
         for (int j = 0; j < prt.m_num_runtime_int; ++j)
             prt.m_runtime_idata[j][i_prt] = initializeIntValue(m_policy_int[j+PartData::NAI]);
    }
};

/**
 * \brief A factory for creating SmartCreate functors.
 *
 * Given a particle container, this can create a functor
 * that will perform the smart create operation on a tile
 * of that particle container
 */
class SmartCreateFactory
{
    PolicyVec m_policy_real;
    PolicyVec m_policy_int;
    bool m_defined;

public:
    template <class PartTileData>
    SmartCreateFactory (const PartTileData& part) noexcept
        : m_defined(false)
    {
        m_policy_real = getPolicies(part.getParticleComps());
        m_policy_int  = getPolicies(part.getParticleiComps());
        m_defined = true;
    }

    SmartCreate getSmartCreate () const noexcept
    {
        AMREX_ASSERT(m_defined);
        return SmartCreate{m_policy_real.dataPtr(),
                         m_policy_int.dataPtr()};
    }

    bool isDefined () const noexcept { return m_defined; }
};

#endif //SMART_CREATE_H_
