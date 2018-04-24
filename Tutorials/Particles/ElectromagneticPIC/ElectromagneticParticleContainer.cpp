#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

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
    : ParticleContainer<0, 0, PIdx::nattribs> (a_geom, a_dmap, a_ba),
      m_species_id(a_species_id), m_charge(a_charge), m_mass(a_mass)
{
}

void
ElectromagneticParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std, 
              const Real     a_thermal_momentum_mean,
              const Real     a_density)
{
    BL_PROFILE("ElectromagneticParticleContainer::InitParticles");
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();

    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0], 
                                     *a_num_particles_per_cell[1], 
                                     *a_num_particles_per_cell[2]);
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
    
    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const Real* tile_lo = tile_realbox.lo();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[3];
                Real u[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                get_gaussian_random_momentum(u, a_thermal_momentum_mean, a_thermal_momentum_std);
                
                Real x = tile_lo[0] + (iv[0] + r[0])*dx[0];
                Real y = tile_lo[1] + (iv[1] + r[1])*dx[1];
                Real z = tile_lo[2] + (iv[2] + r[2])*dx[2];
                
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;
    
                attribs[PIdx::w ] = a_density * scale_fac;
                attribs[PIdx::ux] = u[0] * PhysConst::c;
                attribs[PIdx::uy] = u[1] * PhysConst::c;
                attribs[PIdx::uz] = u[2] * PhysConst::c;
                
                particle_tile.push_back(p);
                particle_tile.push_back_real(attribs);                
            }
        }
    }
}

void
ElectromagneticParticleContainer::
PushAndDeposeParticles(const amrex::MultiFab& Ex,
                       const amrex::MultiFab& Ey,
                       const amrex::MultiFab& Ez,
                       const amrex::MultiFab& Bx,
                       const amrex::MultiFab& By,
                       const amrex::MultiFab& Bz,
                             amrex::MultiFab& jx, 
                             amrex::MultiFab& jy, 
                             amrex::MultiFab& jz,
                             amrex::Real      dt)
{   
    BL_PROFILE("ElectromagneticParticleContainer::PushAndDeposeParticles");
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx  = geom.CellSize();
    const Real* plo = geom.ProbLo();
    
    for (ElectromagneticParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        const int np    = pti.numParticles();
        auto& attribs   = pti.GetAttribs();

        auto& uxp  = attribs[PIdx::ux  ];
        auto& uyp  = attribs[PIdx::uy  ];
        auto& uzp  = attribs[PIdx::uz  ];
        auto& Exp  = attribs[PIdx::Ex  ];
        auto& Eyp  = attribs[PIdx::Ey  ];
        auto& Ezp  = attribs[PIdx::Ez  ];
        auto& Bxp  = attribs[PIdx::Bx  ];
        auto& Byp  = attribs[PIdx::By  ];
        auto& Bzp  = attribs[PIdx::Bz  ];
        auto& wp   = attribs[PIdx::w   ];
        auto& ginv = attribs[PIdx::ginv];

        FORT_LAUNCH_PARTICLES(np,
                              gather_magnetic_field,
                              np, particles.data(),
                              Bxp.data(), Byp.data(), Bzp.data(),
                              BL_TO_FORTRAN_3D(Bx[pti]),
                              BL_TO_FORTRAN_3D(By[pti]),
                              BL_TO_FORTRAN_3D(Bz[pti]),
                              plo, dx);
        
        FORT_LAUNCH_PARTICLES(np,
                              gather_electric_field,
                              np, particles.data(),
                              Exp.data(), Eyp.data(), Ezp.data(), 
                              BL_TO_FORTRAN_3D(Ex[pti]),
                              BL_TO_FORTRAN_3D(Ey[pti]),
                              BL_TO_FORTRAN_3D(Ez[pti]),
                              plo, dx);
        
#ifdef AMREX_USE_CUDA           
        cudaDeviceSynchronize();
#endif
        FORT_LAUNCH_PARTICLES(np, 
                              push_momentum_boris,
                              np, uxp.data(), uyp.data(), uzp.data(), ginv.data(),
                              Exp.data(), Eyp.data(), Ezp.data(),
                              Bxp.data(), Byp.data(), Bzp.data(),
                              m_charge, m_mass, dt);
        
#ifdef AMREX_USE_CUDA                        
        cudaDeviceSynchronize();
#endif      
        
        FORT_LAUNCH_PARTICLES(np,
                              push_position_boris,
                              np, particles.data(),
                              uxp.data(), uyp.data(), uzp.data(),
                              ginv.data(), dt);
        
#ifdef AMREX_USE_CUDA            
        cudaDeviceSynchronize();
#endif
        
        FORT_LAUNCH_PARTICLES(np, deposit_current,
                              BL_TO_FORTRAN_3D(jx[pti]),
                              BL_TO_FORTRAN_3D(jy[pti]),
                              BL_TO_FORTRAN_3D(jz[pti]),
                              np, particles.data(),
                              uxp.data(), uyp.data(), uzp.data(),
                              ginv.data(), wp.data(),
                              m_charge, plo, dt, dx);            
    }
}

void 
ElectromagneticParticleContainer::
PushParticlesOnly(amrex::Real dt)
{
    BL_PROFILE("ElectromagneticParticleContainer::PushParticlesOnly");

    const int lev = 0;
    for (ElectromagneticParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        const int np    = pti.numParticles();
        auto& attribs   = pti.GetAttribs();

        auto& ux   = attribs[PIdx::ux  ];
        auto& uy   = attribs[PIdx::uy  ];
        auto& uz   = attribs[PIdx::uz  ];
        auto& ginv = attribs[PIdx::ginv];

        FORT_LAUNCH_PARTICLES(np, set_gamma,
                              np, ux.data(), uy.data(), uz.data(), ginv.data());
        
#ifdef AMREX_USE_CUDA
        cudaDeviceSynchronize();
#endif        
        
        FORT_LAUNCH_PARTICLES(np, push_position_boris,
                              np, particles.data(),
                              ux.data(), uy.data(), uz.data(), ginv.data(), dt);
    }
}

void
ElectromagneticParticleContainer::
EnforcePeriodicBCs()
{
    BL_PROFILE("ElectromagneticParticleContainer::EnforcePeriodicBCs");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    
    for (ElectromagneticParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        const int np    = pti.numParticles();

        const Real* plo = geom.ProbLo();
        const Real* phi = geom.ProbHi();
        
        FORT_LAUNCH_PARTICLES(np, enforce_periodic,
                              np, particles.data(), plo, phi);
    }
}
