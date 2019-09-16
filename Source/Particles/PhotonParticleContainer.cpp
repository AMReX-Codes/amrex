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
        bool test_quantum_sync;
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

#ifdef WARPX_QED
    if(do_qed_breit_wheeler)
        InitTauBreitWheeler();
#endif

    if (maxLevel() > 0) {
        Redistribute();  // We then redistribute
    }
}

void
PhotonParticleContainer::PushPX(WarpXParIter& pti,
                                Cuda::ManagedDeviceVector<Real>& xp,
                                Cuda::ManagedDeviceVector<Real>& yp,
                                Cuda::ManagedDeviceVector<Real>& zp,
                                Real dt)
{

    // This wraps the momentum and position advance so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();
    // Extract pointers to the different particle quantities
    Real* const AMREX_RESTRICT x = xp.dataPtr();
    Real* const AMREX_RESTRICT y = yp.dataPtr();
    Real* const AMREX_RESTRICT z = zp.dataPtr();
    Real* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    Real* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    Real* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
    const Real* const AMREX_RESTRICT Ex = attribs[PIdx::Ex].dataPtr();
    const Real* const AMREX_RESTRICT Ey = attribs[PIdx::Ey].dataPtr();
    const Real* const AMREX_RESTRICT Ez = attribs[PIdx::Ez].dataPtr();
    const Real* const AMREX_RESTRICT Bx = attribs[PIdx::Bx].dataPtr();
    const Real* const AMREX_RESTRICT By = attribs[PIdx::By].dataPtr();
    const Real* const AMREX_RESTRICT Bz = attribs[PIdx::Bz].dataPtr();

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        copy_attribs(pti, x, y, z);
    }

    //No need to update momentum for photons (for now)

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
                                        Real t, Real dt)
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
// A function to initialize the Tau component according to the BW engine
void PhotonParticleContainer::InitTauBreitWheeler()
{
    BL_PROFILE("PhotonParticleContainer::InitOpticalDepth");
    //Looping over all the particles
    int num_levels = finestLevel() + 1;
    for (int lev=0; lev < num_levels; ++lev)
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti){
            auto taus = pti.GetAttribs(particle_comps["tau"]).dataPtr();
            amrex::ParallelFor(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    taus[i] = warpx_breit_wheeler_engine::get_optical_depth();
                }
                );
    }
}
#endif
