#include "MyTest.H"
#include "MyTest_K.H"

#include <AMReX_MLTensorOp.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initData();
}

void
MyTest::solve ()
{
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        mlmg_lobc[idim] = LinOpBCType::Dirichlet;
        mlmg_hibc[idim] = LinOpBCType::Dirichlet;
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    std::unique_ptr<MLTensorOp> mltensor;
    if (do_overset) {
        mltensor.reset(new MLTensorOp({geom}, {grids}, {dmap}, {&oversetmask}, info));
    } else {
        mltensor.reset(new MLTensorOp({geom}, {grids}, {dmap}, info));
    }

    mltensor->setDomainBC(mlmg_lobc, mlmg_hibc);
    mltensor->setLevelBC(0, &exact);

    const Real a = 0.0; // 1.0e6;
    {
        MultiFab tmp(grids, dmap, 1, 0);
        tmp.setVal(a);

        mltensor->setACoeffs(0, tmp);

        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(grids, IntVect::TheDimensionVector(idim));
            face_bcoef[idim].define(ba, dmap, 1, 0);
        }
        amrex::average_cellcenter_to_face(amrex::GetArrOfPtrs(face_bcoef), eta, geom);
        mltensor->setShearViscosity(0, amrex::GetArrOfConstPtrs(face_bcoef));
    }

    MultiFab::Saxpy(rhs, a, exact, 0, 0, AMREX_SPACEDIM, 0);

    MLMG mlmg(*mltensor);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    // In region with overset mask = 0, phi has valid solution and rhs is zero.
    Real mlmg_err = mlmg.solve({&solution}, {&rhs}, 1.e-11, 0.0);
}

void
MyTest::writePlotfile ()
{
    Vector<std::string> varname = {"u", "v", "w", "uexact", "vexact", "wexact",
                                   "xerror", "yerror", "zerror", "xrhs", "yrhs", "zrhs", "eta"};
    MultiFab plotmf(grids, dmap, varname.size(),  0);
    MultiFab::Copy(plotmf, solution, 0,  0, 3, 0);
    MultiFab::Copy(plotmf, exact   , 0,  3, 3, 0);
    MultiFab::Copy(plotmf, solution, 0,  6, 3, 0);
    MultiFab::Copy(plotmf, rhs     , 0,  9, 3, 0);
    MultiFab::Copy(plotmf, eta     , 0, 12, 1, 0);
    MultiFab::Subtract(plotmf, exact, 0, 6, 3, 0);
    WriteMultiLevelPlotfile("plot", 1, {&plotmf},
                            varname, {geom}, 0.0, {0}, {IntVect(2)});
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        amrex::Print() << "\n";
        amrex::Print() << "  max-norm error = " << plotmf.norm0(6+idim) << std::endl;
        const auto dx = geom.CellSize();
        amrex::Print() << "    1-norm error = " << plotmf.norm1(6+idim) * (dx[0]*dx[1]*dx[2])
                       << std::endl;
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("plot_file", plot_file_name);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("do_overset", do_overset);
}

void
MyTest::initGrids ()
{
    RealBox rb({AMREX_D_DECL(-1.,-1.,-1.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, isperiodic.data());
    Box domain(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    geom.define(domain, rb, CoordSys::cartesian, isperiodic);

    grids.define(domain);
    grids.maxSize(max_grid_size);
}

void
MyTest::initData ()
{
    dmap.define(grids);

    solution.define(grids, dmap, AMREX_SPACEDIM, 1);
    exact.define(grids, dmap, AMREX_SPACEDIM, 1);
    rhs.define(grids, dmap, AMREX_SPACEDIM, 1);
    eta.define(grids, dmap, 1, 1);
    oversetmask.define(grids, dmap, 1, 0);

    Box overset_box = amrex::grow(geom.Domain(), -n_cell/4); // middle of the domain
    // Box overset_box = amrex::shift(geom.Domain(), 0, n_cell/2); // right half

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();
    auto loverset = do_overset;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = mfi.growntilebox(1);
        const Array4<Real> solnfab = solution.array(mfi);
        const Array4<Real> exactfab = exact.array(mfi);
        const Array4<Real> rhsfab = rhs.array(mfi);
        const Array4<Real> etafab = eta.array(mfi);
        const Array4<int> mask = oversetmask.array(mfi);

        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5)*dx[0] + problo[0];
            Real y = (j+0.5)*dx[1] + problo[1];
            Real z = (k+0.5)*dx[2] + problo[2];
            x = std::max(-1.0,std::min(1.0,x));
            y = std::max(-1.0,std::min(1.0,y));
            z = std::max(-1.0,std::min(1.0,z));
            Real u,v,w,urhs,vrhs,wrhs,seta;
            init(x,y,z,1.0,u,v,w,urhs,vrhs,wrhs,seta);
            exactfab(i,j,k,0) = u;
            exactfab(i,j,k,1) = v;
            exactfab(i,j,k,2) = w;
            rhsfab(i,j,k,0) = urhs;
            rhsfab(i,j,k,1) = vrhs;
            rhsfab(i,j,k,2) = wrhs;
            etafab(i,j,k) = seta;
            if (vbx.contains(IntVect(i,j,k)) and overset_box.contains(IntVect(i,j,k))) {
                solnfab(i,j,k,0) = u;
                solnfab(i,j,k,1) = v;
                solnfab(i,j,k,2) = w;
            } else {
                solnfab(i,j,k,0) = 0;
                solnfab(i,j,k,1) = 0;
                solnfab(i,j,k,2) = 0;
            }
            if (vbx.contains(IntVect(i,j,k))) {
                if (overset_box.contains(IntVect(i,j,k))) {
                    mask(i,j,k) = 0;
                } else {
                    mask(i,j,k) = 1;
                }
            }
        });
    }
}
