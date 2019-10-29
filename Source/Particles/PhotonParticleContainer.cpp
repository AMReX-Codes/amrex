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
        //IF do_qed is enabled, find out if Breit Wheeler process is enabled
        if(do_qed)
            pp.query("do_qed_breit_wheeler", do_qed_breit_wheeler);

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
                                Cuda::ManagedDeviceVector<ParticleReal>& xp,
                                Cuda::ManagedDeviceVector<ParticleReal>& yp,
                                Cuda::ManagedDeviceVector<ParticleReal>& zp,
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

#ifdef WARPX_QED
     if(has_breit_wheeler()) DoBreitWheelerPti(pti, dt);
#endif

}

void
PhotonParticleContainer::PushP (int lev,
                                amrex::Real dt,
                                const amrex::MultiFab& Ex,
                                const amrex::MultiFab& Ey,
                                const amrex::MultiFab& Ez,
                                const amrex::MultiFab& Bx,
                                const amrex::MultiFab& By,
                                const amrex::MultiFab& Bz)
{
    BL_PROFILE("PhotonParticleContainer::PushP");
    if (do_not_push) return;

#ifdef WARPX_QED
    if(has_breit_wheeler()) DoBreitWheeler(lev,dt, Ex,Ey,Ez,Bx,By,Bz);
#endif
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
PhotonParticleContainer::DoBreitWheeler(int lev,
                                        amrex::Real dt,
                                        const amrex::MultiFab& Ex,
                                        const amrex::MultiFab& Ey,
                                        const amrex::MultiFab& Ez,
                                        const amrex::MultiFab& Bx,
                                        const amrex::MultiFab& By,
                                        const amrex::MultiFab& Bz)
{
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
            auto& Ezp = attribs[PIdx::Ez];
            auto& Bxp = attribs[PIdx::Bx];
            auto& Byp = attribs[PIdx::By];
            auto& Bzp = attribs[PIdx::Bz];

            const long np = pti.numParticles();

            // Data on the grid
            const FArrayBox& exfab = Ex[pti];
            const FArrayBox& eyfab = Ey[pti];
            const FArrayBox& ezfab = Ez[pti];
            const FArrayBox& bxfab = Bx[pti];
            const FArrayBox& byfab = By[pti];
            const FArrayBox& bzfab = Bz[pti];

            Exp.assign(np,WarpX::E_external_particle[0]);
            Eyp.assign(np,WarpX::E_external_particle[1]);
            Ezp.assign(np,WarpX::E_external_particle[2]);

            Bxp.assign(np,WarpX::B_external_particle[0]);
            Byp.assign(np,WarpX::B_external_particle[1]);
            Bzp.assign(np,WarpX::B_external_particle[2]);

            //
            // copy data from particle container to temp arrays
            //
            pti.GetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);

            int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();
            FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                        &exfab, &eyfab, &ezfab, &bxfab, &byfab, &bzfab,
                        Ex.nGrow(), e_is_nodal,
                        0, np, thread_num, lev, lev);

            // This wraps the momentum advance so that inheritors can modify the call.
            // Extract pointers to the different particle quantities
            ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            const ParticleReal* const AMREX_RESTRICT Expp = Exp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Eypp = Eyp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Ezpp = Ezp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Bxpp = Bxp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Bypp = Byp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Bzpp = Bzp.dataPtr();

            BreitWheelerEvolveOpticalDepth evolve_opt =
                shr_ptr_bw_engine->build_evolve_functor();

            amrex::Real* AMREX_RESTRICT p_tau =
                pti.GetAttribs(particle_comps["tau"]).dataPtr();

            const auto me = PhysConst::m_e;

            // Loop over the particles and update their optical_depth
            amrex::ParallelFor(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    evolve_opt(
                        me*ux[i], me*uy[i], me*uz[i],
                        Expp[i], Eypp[i], Ezpp[i],
                        Bxpp[i], Bypp[i], Bzpp[i],
                        dt, p_tau[i]);
                }
            );
        }
    }
}

void
PhotonParticleContainer::DoBreitWheelerPti(WarpXParIter& pti,
                  amrex::Real dt)
{
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
        shr_ptr_bw_engine->build_evolve_functor();

    amrex::Real* AMREX_RESTRICT p_tau =
        pti.GetAttribs(particle_comps["tau"]).dataPtr();

    const auto me = PhysConst::m_e;

    amrex::ParallelFor(
        pti.numParticles(),
        [=] AMREX_GPU_DEVICE (long i) {
            evolve_opt(
                me*ux[i], me*uy[i], me*uz[i],
                Ex[i], Ey[i], Ez[i],
                Bx[i], By[i], Bz[i],
                dt, p_tau[i]);
        }
    );
}
#endif
