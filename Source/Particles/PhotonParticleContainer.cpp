#include <limits>
#include <sstream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <PhotonParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXConst.H>


// Import low-level single-particle kernels
#include <UpdatePositionPhoton.H>


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
PhotonParticleContainer::PushPX(WarpXParIter& pti,
                                Gpu::ManagedDeviceVector<ParticleReal>& xp,
                                Gpu::ManagedDeviceVector<ParticleReal>& yp,
                                Gpu::ManagedDeviceVector<ParticleReal>& zp,
                                Real dt, DtType a_dt_type)
{

    // This wraps the momentum and position advance so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();
    // Extract pointers to the different particle quantities
    ParticleReal* const AMREX_RESTRICT x = xp.dataPtr();
    ParticleReal* const AMREX_RESTRICT y = yp.dataPtr();
    ParticleReal* const AMREX_RESTRICT z = zp.dataPtr();
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
        copy_attribs(pti, x, y, z);
    }

    amrex::ParallelFor(
        pti.numParticles(),
        [=] AMREX_GPU_DEVICE (long i) {

            UpdatePositionPhoton( x[i], y[i], z[i],
                            ux[i], uy[i], uz[i], dt );
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
                                 Real t, Real dt, DtType a_dt_type)
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

    amrex::Real* AMREX_RESTRICT p_tau =
        pti.GetAttribs(particle_comps["tau"]).dataPtr();

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
                dt, p_tau[i]);
        }
    );
}
#endif
