#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
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

//
// Solve L(phi) = rhs
//
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

    std::unique_ptr<MLABecLaplacian> mlabec;
    if (do_overset) {
        mlabec.reset(new MLABecLaplacian({geom}, {grids}, {dmap}, {&oversetmask}, info));
    } else {
        mlabec.reset(new MLABecLaplacian({geom}, {grids}, {dmap}, info));
    }

    mlabec->setDomainBC(mlmg_lobc, mlmg_hibc);
    mlabec->setLevelBC(0, &exact_phi);

    mlabec->setScalars(ascalar, bscalar);
    mlabec->setACoeffs(0, acoef);

    Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        const BoxArray& ba = amrex::convert(bcoef.boxArray(),
                                            IntVect::TheDimensionVector(idim));
        face_bcoef[idim].define(ba, bcoef.DistributionMap(), 1, 0);
    }
    amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                      bcoef, geom);
    mlabec->setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

    MLMG mlmg(*mlabec);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    // In region with overset mask = 0, phi has valid solution and rhs is zero.
    Real mlmg_err = mlmg.solve({&phi}, {&rhs}, 1.e-11, 0.0);
}

void
MyTest::writePlotfile ()
{
    Vector<std::string> varname = {"solution", "rhs", "exact_solution", "error", "acoef", "bcoef"};
    MultiFab plotmf(grids, dmap, varname.size(), 0);
    MultiFab::Copy(plotmf, phi       , 0, 0, 1, 0);
    MultiFab::Copy(plotmf, rhs       , 0, 1, 1, 0);
    MultiFab::Copy(plotmf, exact_phi , 0, 2, 1, 0);
    MultiFab::Copy(plotmf, phi       , 0, 3, 1, 0);
    MultiFab::Subtract(plotmf, plotmf, 2, 3, 1, 0); // error = soln - exact
    MultiFab::Copy(plotmf, acoef     , 0, 4, 1, 0);
    MultiFab::Copy(plotmf, bcoef     , 0, 5, 1, 0);
    auto dx = geom.CellSize();
    Real dvol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
    amrex::Print() << " max-norm error: " << plotmf.norminf(3)
                   << " 1-norm error: " << plotmf.norm1(3)*dvol << std::endl;
    WriteSingleLevelPlotfile("plot", plotmf, varname, geom, 0.0, 0);
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
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
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

    phi.define(grids, dmap, 1, 1);
    rhs.define(grids, dmap, 1, 0);
    exact_phi.define(grids, dmap, 1, 1);
    acoef.define(grids, dmap, 1, 0);
    bcoef.define(grids, dmap, 1, 1);
    oversetmask.define(grids, dmap, 1, 0);

    Box overset_box = amrex::grow(geom.Domain(), -n_cell/4); // middle of the domain
    // Box overset_box = amrex::shift(geom.Domain(), 0, n_cell/2); // right half

    const auto prob_lo = geom.ProbLoArray();
    const auto prob_hi = geom.ProbHiArray();
    const auto dx      = geom.CellSizeArray();
    auto a = ascalar;
    auto b = bscalar;
    auto loverset = do_overset;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox(1);

        auto phifab = phi.array(mfi);
        auto rhsfab = rhs.array(mfi);
        auto exact = exact_phi.array(mfi);
        auto alpha = acoef.array(mfi);
        auto beta = bcoef.array(mfi);
        auto mask = oversetmask.array(mfi);

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr amrex::Real w = 0.05;
            constexpr amrex::Real sigma = 10.;
            const amrex::Real theta = 0.5*std::log(3.) / (w + 1.e-50);

            constexpr amrex::Real pi = 3.1415926535897932;
            constexpr amrex::Real tpi =  2.*pi;
            constexpr amrex::Real fpi =  4.*pi;
            constexpr amrex::Real fac = static_cast<amrex::Real>(AMREX_SPACEDIM)*4.*pi*pi;

            amrex::Real xc = (prob_hi[0] + prob_lo[0])*0.5;
            amrex::Real yc = (prob_hi[1] + prob_lo[1])*0.5;
#if (AMREX_SPACEDIM == 2)
            amrex::Real zc = 0.0;
#else
            amrex::Real zc = (prob_hi[2] + prob_lo[2])*0.5;
#endif

            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
#if (AMREX_SPACEDIM == 2)
            amrex::Real z = 0.0;
#else
            amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);
#endif

            amrex::Real r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
            amrex::Real tmp = std::cosh(theta*(r-0.25));
            amrex::Real dbdrfac = (sigma-1.)/2./(tmp*tmp) * theta/r;
            dbdrfac *= b;

            // for domain boundary
            x = amrex::min(prob_hi[0], amrex::max(prob_lo[0], x));
            y = amrex::min(prob_hi[1], amrex::max(prob_lo[1], y));
#if (AMREX_SPACEDIM == 3)
            z = amrex::min(prob_hi[2], amrex::max(prob_lo[2], z));
#endif

            beta(i,j,k) = (sigma-1.)/2.*std::tanh(theta*(r-0.25)) + (sigma+1.)/2.;
            exact(i,j,k) = std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                   + .25 * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z);
            phifab(i,j,k) = 0.0;

            if (vbx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                alpha(i,j,k) = 1.;
                rhsfab(i,j,k) = beta(i,j,k)*b*fac*(std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                                 + std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z))
                            + dbdrfac*((x-xc)*(tpi*std::sin(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                              + pi*std::sin(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z))
                                     + (y-yc)*(tpi*std::cos(tpi*x) * std::sin(tpi*y) * std::cos(tpi*z)
                                              + pi*std::cos(fpi*x) * std::sin(fpi*y) * std::cos(fpi*z))
                                     + (z-zc)*(tpi*std::cos(tpi*x) * std::cos(tpi*y) * std::sin(tpi*z)
                                              + pi*std::cos(fpi*x) * std::cos(fpi*y) * std::sin(fpi*z)))
                                            + a * (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                          + 0.25 * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));
                if (loverset and overset_box.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    mask(i,j,k) = 0;
                    phifab(i,j,k) = exact(i,j,k);
                } else {
                    mask(i,j,k) = 1;
                }
            }
        });
    }
}
