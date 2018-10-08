#include "Evolve.H"
#include "NodalFlags.H"
#include "Constants.H"
#include "em_pic_F.H"

#include <AMReX_CudaManaged.H>
#include <AMReX_Device.H>
#include <AMReX_CudaUtility.H>

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
        const Box& tbx = mfi.tilebox(YeeGrid::Ex_nodal_flag);
        const Box& tby = mfi.tilebox(YeeGrid::Ey_nodal_flag);
        const Box& tbz = mfi.tilebox(YeeGrid::Ez_nodal_flag);
        const Box& vbx = mfi.validbox();

        FArrayBox* elecX = &(Ex[mfi]);
        FArrayBox* elecY = &(Ey[mfi]);
        FArrayBox* elecZ = &(Ez[mfi]);
        const FArrayBox* magX  = &(Bx[mfi]);
        const FArrayBox* magY  = &(By[mfi]);
        const FArrayBox* magZ  = &(Bz[mfi]);
        const FArrayBox* currDenX = &(jx[mfi]);
        const FArrayBox* currDenY = &(jy[mfi]);
        const FArrayBox* currDenZ = &(jz[mfi]);

        auto electric = [=] AMREX_CUDA_DEVICE ()
        {
            Box threadBox = getThreadBox(vbx, tbx.type()); 

            if (threadBox.ok())
            {
                push_electric_field_x(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*elecX),
                                      BL_TO_FORTRAN_3D(*magY),
                                      BL_TO_FORTRAN_3D(*magZ),
                                      BL_TO_FORTRAN_3D(*currDenX),
                                      mu_c2_dt, dtsdx_c2[1], dtsdx_c2[2]);
            }

            threadBox = getThreadBox(vbx, tby.type()); 
            if (threadBox.ok())
            {
                push_electric_field_y(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*elecY),
                                      BL_TO_FORTRAN_3D(*magX),
                                      BL_TO_FORTRAN_3D(*magZ),
                                      BL_TO_FORTRAN_3D(*currDenY),
                                      mu_c2_dt, dtsdx_c2[0], dtsdx_c2[2]);
            }

            threadBox = getThreadBox(vbx, tbz.type()); 
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

        AMREX_BOX_L_LAUNCH(vbx, electric);

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
        const Box& vbx = mfi.validbox();

        const FArrayBox* elecX = &(Ex[mfi]);
        const FArrayBox* elecY = &(Ey[mfi]);
        const FArrayBox* elecZ = &(Ez[mfi]);
        FArrayBox* magX  = &(Bx[mfi]);
        FArrayBox* magY  = &(By[mfi]);
        FArrayBox* magZ  = &(Bz[mfi]);

        auto magnetic = [=] AMREX_CUDA_DEVICE ()
        {
            Box threadBox = getThreadBox(vbx, tbx.type()); 

            if (threadBox.ok())
            {
                push_magnetic_field_x(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*magX),
                                      BL_TO_FORTRAN_3D(*elecY),
                                      BL_TO_FORTRAN_3D(*elecZ),
                                      dtsdx[1], dtsdx[2]);
            }

            threadBox = getThreadBox(vbx, tby.type()); 
            if (threadBox.ok())
            {
                push_magnetic_field_y(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*magY),
                                      BL_TO_FORTRAN_3D(*elecX),
                                      BL_TO_FORTRAN_3D(*elecZ),
                                      dtsdx[0], dtsdx[2]);
            }

            threadBox = getThreadBox(vbx, tbz.type()); 
            if (threadBox.ok())
            {
                push_magnetic_field_z(BL_TO_FORTRAN_BOX(threadBox),
                                      BL_TO_FORTRAN_3D(*magZ),
                                      BL_TO_FORTRAN_3D(*elecX),
                                      BL_TO_FORTRAN_3D(*elecY),
                                      dtsdx[0], dtsdx[1]);
            }
        };

        AMREX_BOX_L_LAUNCH(vbx, magnetic);
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
