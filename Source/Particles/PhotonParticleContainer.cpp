/* Copyright 2019 David Grote, Luca Fedeli, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PhotonParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

// Import low-level single-particle kernels
#include "Particles/Pusher/UpdatePositionPhoton.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#ifdef _OPENMP
#   include <omp.h>
#endif

#include <limits>
#include <sstream>
#include <algorithm>


using namespace amrex;

PhotonParticleContainer::PhotonParticleContainer (AmrCore* amr_core, int ispecies,
                                                  const std::string& name)
    : PhysicalParticleContainer(amr_core, ispecies, name)
{
    ParmParse pp(species_name);

#ifdef WARPX_QED
        //IF m_do_qed is enabled, find out if Breit Wheeler process is enabled
        if(m_do_qed)
            pp.query("do_qed_breit_wheeler", m_do_qed_breit_wheeler);

        //If Breit Wheeler process is enabled, look for the target electron and positron
        //species
        if(m_do_qed_breit_wheeler){
            pp.get("qed_breit_wheeler_ele_product_species", m_qed_breit_wheeler_ele_product_name);
            pp.get("qed_breit_wheeler_pos_product_species", m_qed_breit_wheeler_pos_product_name);
        }

        //Check for processes which do not make sense for photons
        bool test_quantum_sync = false;
        pp.query("do_qed_quantum_sync", test_quantum_sync);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        test_quantum_sync == 0,
        "ERROR: do_qed_quantum_sync can be 1 for species NOT listed in particles.photon_species only!");
        //_________________________________________________________
#endif

}

void PhotonParticleContainer::InitData()
{
    AddParticles(0); // Note - add on level 0

    Redistribute();  // We then redistribute

}

void
PhotonParticleContainer::PushPX(WarpXParIter& pti, Real dt, DtType /*a_dt_type*/)
{

    // This wraps the momentum and position advance so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();

    // Extract pointers to the different particle quantities
    ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ex = attribs[PIdx::Ex].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ey = attribs[PIdx::Ey].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ez = attribs[PIdx::Ez].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bx = attribs[PIdx::Bx].dataPtr();
    const ParticleReal* const AMREX_RESTRICT By = attribs[PIdx::By].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bz = attribs[PIdx::Bz].dataPtr();

    if (WarpX::do_back_transformed_diagnostics && do_back_transformed_diagnostics)
    {
        copy_attribs(pti);
    }

    const auto GetPosition = GetParticlePosition(pti);
          auto SetPosition = SetParticlePosition(pti);

    amrex::ParallelFor(
        pti.numParticles(),
        [=] AMREX_GPU_DEVICE (long i) {
            ParticleReal x, y, z;
            GetPosition(i, x, y, z);
            UpdatePositionPhoton( x, y, z, ux[i], uy[i], uz[i], dt );
            SetPosition(i, x, y, z);
        }
    );
}

void
PhotonParticleContainer::Evolve (int lev,
                                 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                 MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                 MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                 MultiFab* rho, MultiFab* crho,
                                 const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                 const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                 Real t, Real dt, DtType /*a_dt_type*/)
{
    // This does gather, push and depose.
    // Push and depose have been re-written for photon,
    // so they do not do anything.
    PhysicalParticleContainer::Evolve (lev,
                                       Ex, Ey, Ez,
                                       Bx, By, Bz,
                                       jx, jy, jz,
                                       cjx, cjy, cjz,
                                       rho, crho,
                                       cEx, cEy, cEz,
                                       cBx, cBy, cBz,
                                       t, dt);

}

#ifdef WARPX_QED

void
PhotonParticleContainer::EvolveOpticalDepth(
    WarpXParIter& pti,amrex::Real dt)
{
    if(!has_breit_wheeler())
        return;

    auto& attribs = pti.GetAttribs();
    ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ex = attribs[PIdx::Ex].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ey = attribs[PIdx::Ey].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ez = attribs[PIdx::Ez].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bx = attribs[PIdx::Bx].dataPtr();
    const ParticleReal* const AMREX_RESTRICT By = attribs[PIdx::By].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bz = attribs[PIdx::Bz].dataPtr();

    BreitWheelerEvolveOpticalDepth evolve_opt =
        m_shr_p_bw_engine->build_evolve_functor();

    amrex::Real* AMREX_RESTRICT p_optical_depth_BW =
        pti.GetAttribs(particle_comps["optical_depth_BW"]).dataPtr();

    const auto me = PhysConst::m_e;

    amrex::ParallelFor(
        pti.numParticles(),
        [=] AMREX_GPU_DEVICE (long i) {
            const ParticleReal px = me * ux[i];
            const ParticleReal py = me * uy[i];
            const ParticleReal pz = me * uz[i];

            bool has_event_happened = evolve_opt(
                px, py, pz,
                Ex[i], Ey[i], Ez[i],
                Bx[i], By[i], Bz[i],
                dt, p_optical_depth_BW[i]);
        }
        );
}

#endif
