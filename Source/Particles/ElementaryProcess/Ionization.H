/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef IONIZATION_H_
#define IONIZATION_H_

#include "Utils/WarpXConst.H"
#include "Particles/WarpXParticleContainer.H"

struct IonizationFilterFunc
{
    const amrex::Real* const AMREX_RESTRICT m_ionization_energies;
    const amrex::Real* const AMREX_RESTRICT m_adk_prefactor;
    const amrex::Real* const AMREX_RESTRICT m_adk_exp_prefactor;
    const amrex::Real* const AMREX_RESTRICT m_adk_power;

    int comp;
    int m_atomic_number;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int i) const noexcept
    {
        const int ion_lev = ptd.m_runtime_idata[comp][i];
        if (ion_lev < m_atomic_number)
        {
            constexpr amrex::Real c = PhysConst::c;
            constexpr amrex::Real c2_inv = 1./c/c;

            // Compute electric field amplitude in the particle's frame of
            // reference (particularly important when in boosted frame).
            amrex::ParticleReal ux = ptd.m_rdata[PIdx::ux][i];
            amrex::ParticleReal uy = ptd.m_rdata[PIdx::uy][i];
            amrex::ParticleReal uz = ptd.m_rdata[PIdx::uz][i];
            amrex::ParticleReal ex = ptd.m_rdata[PIdx::Ex][i];
            amrex::ParticleReal ey = ptd.m_rdata[PIdx::Ey][i];
            amrex::ParticleReal ez = ptd.m_rdata[PIdx::Ez][i];
            amrex::ParticleReal bx = ptd.m_rdata[PIdx::Bx][i];
            amrex::ParticleReal by = ptd.m_rdata[PIdx::By][i];
            amrex::ParticleReal bz = ptd.m_rdata[PIdx::Bz][i];

            amrex::Real ga = std::sqrt(1. + (ux*ux + uy*uy + uz*uz) * c2_inv);
            amrex::Real E = std::sqrt(
                               - ( ux*ex + uy*ey + uz*ez ) * ( ux*ex + uy*ey + uz*ez ) * c2_inv
                               + ( ga   *ex + uy*bz - uz*by ) * ( ga   *ex + uy*bz - uz*by )
                               + ( ga   *ey + uz*bx - ux*bz ) * ( ga   *ey + uz*bx - ux*bz )
                               + ( ga   *ez + ux*by - uy*bx ) * ( ga   *ez + ux*by - uy*bx )
                               );

            // Compute probability of ionization p
            amrex::Real w_dtau = 1./ ga * m_adk_prefactor[ion_lev] *
                std::pow(E, m_adk_power[ion_lev]) *
                std::exp( m_adk_exp_prefactor[ion_lev]/E );
            amrex::Real p = 1. - std::exp( - w_dtau );

            amrex::Real random_draw = amrex::Random();
            if (random_draw < p)
            {
                return true;
            }
        }
        return false;
    }
};

struct IonizationTransformFunc
{
    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, SrcData& src, int i_src, int i_dst) const noexcept
    {
        src.m_runtime_idata[0][i_src] += 1;
    }
};

#endif
