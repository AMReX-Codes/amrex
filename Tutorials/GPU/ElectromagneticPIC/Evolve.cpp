#include "Evolve.H"
#include "NodalFlags.H"
#include "Constants.H"
#include "em_pic_F.H"

#include <AMReX_Managed.H>
#include <AMReX_Device.H>
#include <AMReX_CUDA_Utility.H>

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
    const CudaArray<Real,3> dtsdx_c2 {foo/dx[0], foo/dx[1], foo/dx[2]};
    
    for (MFIter mfi(Ex); mfi.isValid(); ++mfi)
    {
        const Box& tbx   = mfi.tilebox(YeeGrid::Ex_nodal_flag);
        const Box& tby   = mfi.tilebox(YeeGrid::Ey_nodal_flag);
        const Box& tbz   = mfi.tilebox(YeeGrid::Ez_nodal_flag);
        const Box & cc = convert(mfi.validbox(), IntVect(0,0,0));

        BaseFab<Real>* elecX = &(Ex[mfi]);
        BaseFab<Real>* elecY = &(Ey[mfi]);
        BaseFab<Real>* elecZ = &(Ez[mfi]);
        const BaseFab<Real>* magX  = &(Bx[mfi]);
        const BaseFab<Real>* magY  = &(By[mfi]);
        const BaseFab<Real>* magZ  = &(Bz[mfi]);
        const BaseFab<Real>* currDenX = &(jx[mfi]);
        const BaseFab<Real>* currDenY = &(jy[mfi]);
        const BaseFab<Real>* currDenZ = &(jz[mfi]);

        auto electric = [=] AMREX_CUDA_DEVICE ()
        {
            Box threadBox = getThreadBox(cc, tbx.type()); 

            if (threadBox.ok())
            {
                push_electric_field_x(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*elecX),
                                      BL_TO_FORTRAN_3D(*magY),
                                      BL_TO_FORTRAN_3D(*magZ),
                                      BL_TO_FORTRAN_3D(*currDenX),
                                      mu_c2_dt, dtsdx_c2[1], dtsdx_c2[2]);
            }

            threadBox = getThreadBox(cc, tby.type()); 
            if (threadBox.ok())
            {
                push_electric_field_y(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*elecY),
                                      BL_TO_FORTRAN_3D(*magX),
                                      BL_TO_FORTRAN_3D(*magZ),
                                      BL_TO_FORTRAN_3D(*currDenY),
                                      mu_c2_dt, dtsdx_c2[0], dtsdx_c2[2]);
            }

            threadBox = getThreadBox(cc, tbz.type()); 
            if (threadBox.ok())
            {
                push_electric_field_z(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*elecZ),
                                      BL_TO_FORTRAN_3D(*magX),
                                      BL_TO_FORTRAN_3D(*magY),
                                      BL_TO_FORTRAN_3D(*currDenZ),
                                      mu_c2_dt, dtsdx_c2[0], dtsdx_c2[1]);
            }
        };

        AMREX_BOX_L_LAUNCH(cc, electric);

    }
}

void evolve_magnetic_field(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                 MultiFab& Bx,       MultiFab& By,       MultiFab& Bz,
                           const Geometry& geom, Real dt)
{
    BL_PROFILE("evolve_magnetic_field");
    
    const Real* dx = geom.CellSize();
    const CudaArray<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    for (MFIter mfi(Bx); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox(YeeGrid::Bx_nodal_flag);
        const Box& tby = mfi.tilebox(YeeGrid::By_nodal_flag);
        const Box& tbz = mfi.tilebox(YeeGrid::Bz_nodal_flag);
        const Box & cc = convert(mfi.validbox(), IntVect(0,0,0));

        const BaseFab<Real>* elecX = &(Ex[mfi]);
        const BaseFab<Real>* elecY = &(Ey[mfi]);
        const BaseFab<Real>* elecZ = &(Ez[mfi]);
        BaseFab<Real>* magX  = &(Bx[mfi]);
        BaseFab<Real>* magY  = &(By[mfi]);
        BaseFab<Real>* magZ  = &(Bz[mfi]);

        auto magnetic = [=] AMREX_CUDA_DEVICE ()
        {
            Box threadBox = getThreadBox(cc, tbx.type()); 

            if (threadBox.ok())
            {
                push_magnetic_field_x(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*magX),
                                      BL_TO_FORTRAN_3D(*elecY),
                                      BL_TO_FORTRAN_3D(*elecZ),
                                      dtsdx[1], dtsdx[2]);
            }

            threadBox = getThreadBox(cc, tby.type()); 
            if (threadBox.ok())
            {
                push_magnetic_field_y(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*magY),
                                      BL_TO_FORTRAN_3D(*elecX),
                                      BL_TO_FORTRAN_3D(*elecZ),
                                      dtsdx[0], dtsdx[2]);
            }

            threadBox = getThreadBox(cc, tbz.type()); 
            if (threadBox.ok())
            {
                push_magnetic_field_z(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*magZ),
                                      BL_TO_FORTRAN_3D(*elecX),
                                      BL_TO_FORTRAN_3D(*elecY),
                                      dtsdx[0], dtsdx[1]);
            }
        };

        AMREX_BOX_L_LAUNCH(cc, magnetic);
    } 
}

void fill_boundary_electric_field(MultiFab& Ex, MultiFab& Ey, MultiFab& Ez, const Geometry& geom)
{    
    const auto& period = geom.periodicity();
    Ex.FillBoundary(period);
    Ey.FillBoundary(period);
    Ez.FillBoundary(period);    
}

void fill_boundary_magnetic_field(MultiFab& Bx, MultiFab& By, MultiFab& Bz, const Geometry& geom)
{    
    const auto& period = geom.periodicity();
    Bx.FillBoundary(period);
    By.FillBoundary(period);
    Bz.FillBoundary(period);    
}
