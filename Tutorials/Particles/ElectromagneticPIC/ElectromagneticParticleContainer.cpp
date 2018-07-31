#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_RedistributeStrategy.H>

#include "ElectromagneticParticleContainer.H"
#include "Constants.H"

#include "em_pic_F.H"

using namespace amrex;

namespace {
    
    void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
        int ny = nppc[1];
        int nz = nppc[2];
        
        int ix_part = i_part/(ny * nz);
        int iy_part = (i_part % (ny * nz)) % ny;
        int iz_part = (i_part % (ny * nz)) / ny;
        
        r[0] = (0.5+ix_part)/nx;
        r[1] = (0.5+iy_part)/ny;
        r[2] = (0.5+iz_part)/nz;    
    }
    
    void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
        Real ux_th = amrex::RandomNormal(0.0, u_std);
        Real uy_th = amrex::RandomNormal(0.0, u_std);
        Real uz_th = amrex::RandomNormal(0.0, u_std);
        
        u[0] = u_mean + ux_th;
        u[1] = u_mean + uy_th;
        u[2] = u_mean + uz_th;
    }
}

ElectromagneticParticleContainer::
ElectromagneticParticleContainer(const Geometry            & a_geom,
                                 const DistributionMapping & a_dmap,
                                 const BoxArray            & a_ba,
                                 const int                   a_species_id,
                                 const Real                  a_charge,
                                 const Real                  a_mass)
    : m_ba(a_ba), m_geom(a_geom), m_dmap(a_dmap),
      m_species_id(a_species_id), m_charge(a_charge), m_mass(a_mass)
{
    BL_PROFILE("ElectromagneticParticleContainer::ElectromagneticParticleContainer");
    
    const int ng = 1;
    m_mask_ptr.reset(new MultiFab(m_ba, m_dmap, 1, ng));
    m_mask_ptr->setVal(-1);
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const int grid_id = mfi.index();
        m_mask_ptr->setVal(grid_id, box, 0, 1);
    }
    m_mask_ptr->FillBoundary(m_geom.periodicity());
    
    m_redistribute_strategy.reset(new RedistributeStrategyCPU() );
}

void
ElectromagneticParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std, 
              const Real     a_thermal_momentum_mean,
              const Real     a_density, 
              const RealBox& a_bounds, 
              const int      a_problem)
{
    BL_PROFILE("ElectromagneticParticleContainer::InitParticles");
    
    const Real* dx = m_geom.CellSize();

    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0], 
                                     *a_num_particles_per_cell[1], 
                                     *a_num_particles_per_cell[2]);
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
    
    for(MFIter mfi(*m_mask_ptr); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const Real* plo = m_geom.ProbLo();
        const int grid_id = mfi.index();
        auto& particles = m_particles[grid_id];
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[3];
                Real u[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                if (a_problem == 0) {
                    get_gaussian_random_momentum(u, a_thermal_momentum_mean, a_thermal_momentum_std);
                }
                else if (a_problem == 1 ) {
                    u[0] = 0.01;
                    u[1] = 0.0;
                    u[2] = 0.0;
                } else {
                    amrex::Abort("problem type not valid");
                }

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;

                particles.x().push_back(x);
                particles.y().push_back(y);
                particles.z().push_back(z);
                
                particles.ux().push_back(u[0] * PhysConst::c);
                particles.uy().push_back(u[1] * PhysConst::c);
                particles.uz().push_back(u[2] * PhysConst::c);
                
                particles.w().push_back(a_density * scale_fac);
                
                particles.ex().push_back(0);
                particles.ey().push_back(0);
                particles.ez().push_back(0);
                particles.bx().push_back(0);
                particles.by().push_back(0);
                particles.bz().push_back(0);
                particles.ginv().push_back(0);
            }
        }
    }
}

void ElectromagneticParticleContainer::PushAndDeposeParticles(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                                              const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                                              MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("ElectromagneticParticleContainer::PushAndDeposeParticles");
    
    const Real* dx  = m_geom.CellSize();
    const Real* plo = m_geom.ProbLo();
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.numParticles();

        if (np == 0) continue;

#ifdef AMREX_USE_ACC
	FTOC(gather_magnetic_field)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np,
                                    gather_magnetic_field,
#endif
                                    np, 
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.bx().data(), particles.by().data(), particles.bz().data(),
                                    BL_TO_FORTRAN_3D(Bx[mfi]),
                                    BL_TO_FORTRAN_3D(By[mfi]),
                                    BL_TO_FORTRAN_3D(Bz[mfi]),
                                    plo, dx);
        
#ifdef AMREX_USE_ACC
	FTOC(gather_electric_field)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np,
                                    gather_electric_field,
#endif
                                    np, 
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.ex().data(), particles.ey().data(), particles.ez().data(),
                                    BL_TO_FORTRAN_3D(Ex[mfi]),
                                    BL_TO_FORTRAN_3D(Ey[mfi]),
                                    BL_TO_FORTRAN_3D(Ez[mfi]),
                                    plo, dx);
        
#ifdef AMREX_USE_ACC
        FTOC(push_momentum_boris)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np, 
                                    push_momentum_boris,
#endif
                                    np,
                                    particles.ux().data(), particles.uy().data(), particles.uz().data(),
                                    particles.ginv().data(),
                                    particles.ex().data(), particles.ey().data(), particles.ez().data(),
                                    particles.bx().data(), particles.by().data(), particles.bz().data(),
                                    m_charge, m_mass, dt);
        
#ifdef AMREX_USE_ACC
        FTOC(push_position_boris)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np,
                                    push_position_boris,
#endif
                                    np,
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.ux().data(), particles.uy().data(), particles.uz().data(),
                                    particles.ginv().data(), dt);
        
#ifdef AMREX_USE_ACC
        FTOC(deposit_current)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np, 
                                    deposit_current,
#endif
                                    BL_TO_FORTRAN_3D(jx[mfi]),
                                    BL_TO_FORTRAN_3D(jy[mfi]),
                                    BL_TO_FORTRAN_3D(jz[mfi]),
                                    np,
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.ux().data(), particles.uy().data(), particles.uz().data(),
                                    particles.ginv().data(), particles.w().data(),
                                    m_charge, plo, dt, dx);
        }
}

void ElectromagneticParticleContainer::PushParticleMomenta(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                                           const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz, Real dt)
{   
    BL_PROFILE("ElectromagneticParticleContainer::PushParticleMomenta");
    
    const Real* dx  = m_geom.CellSize();
    const Real* plo = m_geom.ProbLo();
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.numParticles();

        if (np == 0) continue;

#ifdef AMREX_USE_ACC
	FTOC(gather_magnetic_field)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np,
                                    gather_magnetic_field,
#endif
                                    np, 
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.bx().data(), particles.by().data(), particles.bz().data(),
                                    BL_TO_FORTRAN_3D(Bx[mfi]),
                                    BL_TO_FORTRAN_3D(By[mfi]),
                                    BL_TO_FORTRAN_3D(Bz[mfi]),
                                    plo, dx);
        
#ifdef AMREX_USE_ACC
	FTOC(gather_electric_field)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np,
                                    gather_electric_field,
#endif
                                    np, 
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.ex().data(), particles.ey().data(), particles.ez().data(),
                                    BL_TO_FORTRAN_3D(Ex[mfi]),
                                    BL_TO_FORTRAN_3D(Ey[mfi]),
                                    BL_TO_FORTRAN_3D(Ez[mfi]),
                                    plo, dx);
        
#ifdef AMREX_USE_ACC
        FTOC(push_momentum_boris)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np, 
                                    push_momentum_boris,
#endif
                                    np,
                                    particles.ux().data(), particles.uy().data(), particles.uz().data(),
                                    particles.ginv().data(),
                                    particles.ex().data(), particles.ey().data(), particles.ez().data(),
                                    particles.bx().data(), particles.by().data(), particles.bz().data(),
                                    m_charge, m_mass, dt);
    }
}

void ElectromagneticParticleContainer::PushParticlePositions(Real dt)
{
    BL_PROFILE("ElectromagneticParticleContainer::PushParticlePositions");
    
    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.numParticles();
        
        if (np == 0) continue;

        AMREX_FORT_LAUNCH_PARTICLES(np, set_gamma,
                                    np, 
                                    particles.ux().data(), particles.uy().data(), particles.uz().data(),
                                    particles.ginv().data());
        
#ifdef AMREX_USE_ACC
        FTOC(push_position_boris)(
#else
        AMREX_FORT_LAUNCH_PARTICLES(np, 
                                    push_position_boris,
#endif
                                    np,
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    particles.ux().data(), particles.uy().data(), particles.uz().data(),
                                    particles.ginv().data(), dt);
    }
}

void ElectromagneticParticleContainer::EnforcePeriodicBCs()
{
    BL_PROFILE("ElectromagneticParticleContainer::EnforcePeriodicBCs");

    for (MFIter mfi(*m_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        auto& particles = m_particles[mfi.index()];
        const int np    = particles.numParticles();

        if (np == 0) continue;

        const Real* plo = m_geom.ProbLo();
        const Real* phi = m_geom.ProbHi();
        
        AMREX_FORT_LAUNCH_PARTICLES(np, enforce_periodic,
                                    np,
                                    particles.x().data(),  particles.y().data(),  particles.z().data(),
                                    plo, phi);
    }
}

void ElectromagneticParticleContainer::OK ()
{
    m_redistribute_strategy->OK(m_particles, m_ba, m_dmap, m_geom, m_mask_ptr.get());
}

void ElectromagneticParticleContainer::Redistribute()
{
    m_redistribute_strategy->Redistribute(m_particles, m_ba, m_dmap, m_geom, m_mask_ptr.get());
} 
