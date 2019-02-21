#include "Evolve.H"
#include "NodalFlags.H"
#include "Constants.H"
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

        push_electric_field_x(BL_TO_FORTRAN_BOX(tbx),
                                    BL_TO_FORTRAN_3D(Ex[mfi]),
                                    BL_TO_FORTRAN_3D(By[mfi]),
                                    BL_TO_FORTRAN_3D(Bz[mfi]),
                                    BL_TO_FORTRAN_3D(jx[mfi]),
                                    mu_c2_dt, dtsdx_c2[1], dtsdx_c2[2]);

        push_electric_field_y(BL_TO_FORTRAN_BOX(tby),
                                    BL_TO_FORTRAN_3D(Ey[mfi]),
                                    BL_TO_FORTRAN_3D(Bx[mfi]),
                                    BL_TO_FORTRAN_3D(Bz[mfi]),
                                    BL_TO_FORTRAN_3D(jy[mfi]),
                                    mu_c2_dt, dtsdx_c2[0], dtsdx_c2[2]);

        push_electric_field_z(BL_TO_FORTRAN_BOX(tbz),
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

        push_magnetic_field_x(BL_TO_FORTRAN_BOX(tbx),
                                    BL_TO_FORTRAN_3D(Bx[mfi]),
                                    BL_TO_FORTRAN_3D(Ey[mfi]),
                                    BL_TO_FORTRAN_3D(Ez[mfi]),
                                    dtsdx[1], dtsdx[2]);

        push_magnetic_field_y(BL_TO_FORTRAN_BOX(tby),
                                    BL_TO_FORTRAN_3D(By[mfi]),
                                    BL_TO_FORTRAN_3D(Ex[mfi]),
                                    BL_TO_FORTRAN_3D(Ez[mfi]),
                                    dtsdx[0], dtsdx[2]);

        push_magnetic_field_z(BL_TO_FORTRAN_BOX(tbz),
                                    BL_TO_FORTRAN_3D(Bz[mfi]),
                                    BL_TO_FORTRAN_3D(Ex[mfi]),
                                    BL_TO_FORTRAN_3D(Ey[mfi]),
                                    dtsdx[0], dtsdx[1]);
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

void check_solution(const MultiFab& jx, const Geometry& geom, Real time)
{
    BL_PROFILE("check_solution");

    const Real* dx = geom.CellSize();

    Box test_box = geom.Domain();
    test_box.setSmall(IntVect(AMREX_D_DECL(2, 2, 2)));
    test_box.setBig(IntVect(AMREX_D_DECL(30, 30, 30)));

    Real max_error = 0.0;
    for(MFIter mfi(jx); mfi.isValid(); ++mfi)
    {
        Real fab_max_error;
        const Box& tbx  = mfi.tilebox();
        check_langmuir_solution(BL_TO_FORTRAN_BOX(tbx),
                                BL_TO_FORTRAN_BOX(test_box),
                                BL_TO_FORTRAN_3D(jx[mfi]), time, &fab_max_error);
        max_error += fab_max_error;
    }

    ParallelDescriptor::ReduceRealMax(max_error);

    amrex::Print() << "Max error is: " << max_error << std::endl;
}
