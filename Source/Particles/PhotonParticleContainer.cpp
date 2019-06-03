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

using namespace amrex;

PhotonParticleContainer::PhotonParticleContainer (AmrCore* amr_core, int ispecies,
                                                  const std::string& name)
    : PhysicalParticleContainer(amr_core, ispecies, name)
{

    // This will read <species>.[...] from the inputs file
    // where <species> is the name of your species
    ParmParse pp(species_name);

    // read <species>.size_in_inches in the input file, and
    // store it into member data.
    pp.query("size_in_inches", size_in_inches);

}

void PhotonParticleContainer::InitData()
{
    AddParticles(0); // Note - add on level 0
    if (maxLevel() > 0) {
        Redistribute();  // We then redistribute
    }
}

void
PhotonParticleContainer::PushPX(WarpXParIter& pti,
                                Cuda::ManagedDeviceVector<Real>& xp,
                                Cuda::ManagedDeviceVector<Real>& yp,
                                Cuda::ManagedDeviceVector<Real>& zp,
                                Cuda::ManagedDeviceVector<Real>& giv,
                                Real dt)
{
    // This wraps the call to warpx_particle_pusher so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();
    auto& uxp = attribs[PIdx::ux];
    auto& uyp = attribs[PIdx::uy];
    auto& uzp = attribs[PIdx::uz];
    auto& Exp = attribs[PIdx::Ex];
    auto& Eyp = attribs[PIdx::Ey];
    auto& Ezp = attribs[PIdx::Ez];
    auto& Bxp = attribs[PIdx::Bx];
    auto& Byp = attribs[PIdx::By];
    auto& Bzp = attribs[PIdx::Bz];
    const long np  = pti.numParticles();

    // Using new pusher for positions    
    const amrex_real zero_mass = 0.0;
    const amrex_real zero_charge = 0.0;
    warpx_particle_pusher_positions(&np,
                      xp.dataPtr(),
                      yp.dataPtr(),
                      zp.dataPtr(),
                      uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
                      giv.dataPtr(),
                      Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                      Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                      &zero_charge, &zero_mass, &dt,
                      &WarpX::particle_pusher_algo);
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
    // Currently, I guess photons do gather fields from the mesh.
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
