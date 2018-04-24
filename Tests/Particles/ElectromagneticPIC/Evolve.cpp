#include "Evolve.H"
#include "NodalFlags.H"
#include "Constants.H"
#include "DeviceManager.H"
#include "em_pic_F.H"

using namespace amrex;

Real compute_dt(const Geometry& geom) 
{
    const static Real cfl = 1.0;
    const Real* dx = geom.CellSize();
    const Real dt  = cfl * 1./( std::sqrt(D_TERM(  1./(dx[0]*dx[0]),
                                                 + 1./(dx[1]*dx[1]),
                                                 + 1./(dx[2]*dx[2]))) * PhysConst::c );
    return dt;
}

void evolve_electric_field(      MultiFab& Ex,       MultiFab& Ey,       MultiFab& Ez,
                           const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                           const MultiFab& jx, const MultiFab& jy, const MultiFab& jz,
                           const Geometry& geom, Real dt)
{
    BL_PROFILE("evolve_electric_field");
    
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
    const Real foo = (PhysConst::c*PhysConst::c) * dt;
    
    const Real* dx = geom.CellSize();    
    const std::array<Real,3> dtsdx_c2 {foo/dx[0], foo/dx[1], foo/dx[2]};
    
    for (MFIter mfi(Ex); mfi.isValid(); ++mfi)
    {
        const Box& tbx  = mfi.tilebox(YeeGrid::Ex_nodal_flag);
        const Box& tby  = mfi.tilebox(YeeGrid::Ey_nodal_flag);
        const Box& tbz  = mfi.tilebox(YeeGrid::Ez_nodal_flag);

        FORT_LAUNCH(tbx, push_electric_field_x,
                    BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_3D(Ex[mfi]),
                    BL_TO_FORTRAN_3D(By[mfi]),
                    BL_TO_FORTRAN_3D(Bz[mfi]),
                    BL_TO_FORTRAN_3D(jx[mfi]),
                    mu_c2_dt, dtsdx_c2[1], dtsdx_c2[2]);

        FORT_LAUNCH(tby, push_electric_field_y,
                    BL_TO_FORTRAN_BOX(tby),
                    BL_TO_FORTRAN_3D(Ey[mfi]),
                    BL_TO_FORTRAN_3D(Bx[mfi]),
                    BL_TO_FORTRAN_3D(Bz[mfi]),
                    BL_TO_FORTRAN_3D(jy[mfi]),
                    mu_c2_dt, dtsdx_c2[0], dtsdx_c2[2]);

        FORT_LAUNCH(tbz, push_electric_field_z,
                    BL_TO_FORTRAN_BOX(tbz),
                    BL_TO_FORTRAN_3D(Ez[mfi]),
                    BL_TO_FORTRAN_3D(Bx[mfi]),
                    BL_TO_FORTRAN_3D(By[mfi]),
                    BL_TO_FORTRAN_3D(jz[mfi]),
                    mu_c2_dt, dtsdx_c2[0], dtsdx_c2[1]);
    }
}

void evolve_magnetic_field(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                 MultiFab& Bx,       MultiFab& By,       MultiFab& Bz,
                           const Geometry& geom, Real dt)
{
    BL_PROFILE("evolve_magnetic_field");
    
    const Real* dx = geom.CellSize();
    const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    for (MFIter mfi(Bx); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox(YeeGrid::Bx_nodal_flag);
        const Box& tby = mfi.tilebox(YeeGrid::By_nodal_flag);
        const Box& tbz = mfi.tilebox(YeeGrid::Bz_nodal_flag);

        FORT_LAUNCH(tbx, push_magnetic_field_x,
                    BL_TO_FORTRAN_BOX(tbx),
                    BL_TO_FORTRAN_3D(Bx[mfi]),
                    BL_TO_FORTRAN_3D(Ey[mfi]),
                    BL_TO_FORTRAN_3D(Ez[mfi]),
                    dtsdx[1], dtsdx[2]);

        FORT_LAUNCH(tby, push_magnetic_field_y,
                    BL_TO_FORTRAN_BOX(tby),
                    BL_TO_FORTRAN_3D(By[mfi]),
                    BL_TO_FORTRAN_3D(Ex[mfi]),
                    BL_TO_FORTRAN_3D(Ez[mfi]),
                    dtsdx[0], dtsdx[2]);

        FORT_LAUNCH(tbz, push_magnetic_field_z,
                    BL_TO_FORTRAN_BOX(tbz),
                    BL_TO_FORTRAN_3D(Bz[mfi]),
                    BL_TO_FORTRAN_3D(Ex[mfi]),
                    BL_TO_FORTRAN_3D(Ey[mfi]),
                    dtsdx[0], dtsdx[1]);
    }    
}

void push_particles_only(Vector<Particles*>& particles, Real dt)
{
    BL_PROFILE("particle_push");
    
    for (int ispec = 0; ispec < 2; ++ispec) {
        
        Particles& species = *particles[ispec];
        
        FORT_LAUNCH_PARTICLES(species.np, set_gamma,
                              species.np_d, 
                              species.ux(), species.uy(), species.uz(),
                              species.ginv());
        
#ifdef AMREX_USE_CUDA
        cudaDeviceSynchronize();
#endif        

        FORT_LAUNCH_PARTICLES(species.np, push_position_boris,
                              species.np_d,
                              species.x(),  species.y(),  species.z(),
                              species.ux(), species.uy(), species.uz(),
                              species.ginv(), dt);
    }
}

void enforce_periodic_bcs(Vector<Particles*>& particles, const Geometry& geom)
{
    BL_PROFILE("enforce_periodic");
    
    for (int ispec = 0; ispec < 2; ++ispec) {

        Particles& species = *particles[ispec];

        const Real* plo = geom.ProbLo();
        const Real* phi = geom.ProbHi();
        
        FORT_LAUNCH_PARTICLES(species.np, enforce_periodic,
                              species.np_d,
                              species.x(),  species.y(),  species.z(),
                              plo, phi);
    }
}
