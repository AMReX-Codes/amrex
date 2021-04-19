#include "MyTest.H"

#include <AMReX_Config.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>


#include <AMReX_EB_slopes_K.H>
#include <AMReX_Slopes_K.H>

#include <cmath>

using namespace amrex;


MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

    initData();
}

MyTest::~MyTest ()
{
}

void
MyTest::compute_gradient ()
{
    int ilev = 0;

    int ncomp = phi[0].nComp();

    const Box& domain_box = geom[ilev].Domain();

    const int domlo_x = domain_box.smallEnd(0);
    const int domhi_x = domain_box.bigEnd(0);
    const int domlo_y = domain_box.smallEnd(1);
    const int domhi_y = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domlo_z = domain_box.smallEnd(2);
    const int domhi_z = domain_box.bigEnd(2);
#endif

    const bool on_x_face = !(geom[0].isPeriodic(0));
    const bool on_y_face = !(geom[0].isPeriodic(1));
#if (AMREX_SPACEDIM == 3)
    const bool on_z_face = !(geom[0].isPeriodic(2));
#endif

    MultiFab dummy(grids[ilev],dmap[ilev],1,0);

    for (MFIter mfi(dummy); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();

        Array4<const Real> const& phi_arr     = phi[ilev].array(mfi);
        Array4<      Real> const& grad_x_arr  = grad_x[ilev].array(mfi);
        Array4<      Real> const& grad_y_arr  = grad_y[ilev].array(mfi);
        Array4<      Real> const& grad_z_arr  = grad_z[ilev].array(mfi);
        Array4<      Real> const& ccentr_arr  = ccentr[ilev].array(mfi);

        Array4<Real const> const& fcx   = (factory[ilev]->getFaceCent())[0]->const_array(mfi);
        Array4<Real const> const& fcy   = (factory[ilev]->getFaceCent())[1]->const_array(mfi);
#if (AMREX_SPACEDIM == 3)
        Array4<Real const> const& fcz   = (factory[ilev]->getFaceCent())[2]->const_array(mfi);
#endif

        Array4<Real const> const& ccent = (factory[ilev]->getCentroid()).array(mfi);

        const FabArray<EBCellFlagFab>* flags = &(factory[ilev]->getMultiEBCellFlagFab());
        Array4<EBCellFlag const> const& flag = flags->const_array(mfi);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            ccentr_arr(i,j,k,n) = ccent(i,j,k,n);

            // There is no need to set these to zero other than it makes using
            // amrvis a lot more friendly.
            if ( flag(i,j,k).isCovered() )
            {
              grad_x_arr(i,j,k,n)  = 0.0;
              grad_y_arr(i,j,k,n)  = 0.0;
              grad_z_arr(i,j,k,n)  = 0.0;

            } else {

            // First get EB-aware slope that doesn't know about extdir
           bool edlo_x = 0;
           bool edhi_x = 0;
           bool edlo_y = 0;
           bool edhi_y = 0;
           bool edlo_z = 0;
           bool edhi_z = 0;

#if (AMREX_SPACEDIM == 2)

           bool needs_bdry_stencil = (on_x_face and on_y_face);

           if (needs_bdry_stencil)
           {
              edlo_x = 1;
              edhi_x = 1;
              edlo_y = 1;
              edhi_y = 1;
           }

           auto slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,
                   phi_arr,ccent,
                   fcx,fcy,flag,
                   edlo_x,edlo_y,
                   edhi_x,edhi_y,
                   domlo_x, domlo_y,
                   domhi_x,domhi_y);

           grad_x_arr(i,j,k,n) = slopes[0];
           grad_y_arr(i,j,k,n) = slopes[1];

#elif (AMREX_SPACEDIM == 3)

           if (on_x_face){
               edlo_x = 1;
               edhi_x = 1;
           }

           if (on_y_face){
               edlo_y = 1;
               edhi_y = 1;
           }

           if (on_z_face){
               edlo_z = 1;
               edhi_z = 1;
           }

           auto slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,
                 phi_arr,ccent,
                 fcx,fcy,fcz,flag,
                 edlo_x,edlo_y,edlo_z,
                 edhi_x,edhi_y,edhi_z,
                 domlo_x,domlo_y,domlo_z,
                 domhi_x,domhi_y,domhi_z);

           grad_x_arr(i,j,k,n) = slopes[0];
           grad_y_arr(i,j,k,n) = slopes[1];
           grad_z_arr(i,j,k,n) = slopes[2];

#endif
           }
        });
    }
}

void
MyTest::writePlotfile ()
{
    Vector<MultiFab> plotmf(max_level+1);
    for (int ilev = 0; ilev <= max_level; ++ilev) {

#if (AMREX_SPACEDIM == 2)

        plotmf[ilev].define(grids[ilev],dmap[ilev],3,0);
        MultiFab::Copy(plotmf[ilev], phi[ilev]    , 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], grad_x[ilev] , 0, 1, 1, 0);
        MultiFab::Copy(plotmf[ilev], grad_y[ilev] , 0, 2, 1, 0);
    }
    WriteMultiLevelPlotfile(plot_file_name, max_level+1,
                            amrex::GetVecOfConstPtrs(plotmf),
                            {"phi",
                            "dphidx",
                            "dphidy",
                             },
                            geom, 0.0, Vector<int>(max_level+1,0),
                            Vector<IntVect>(max_level,IntVect{2}));

    Vector<MultiFab> plotmf_analytic(max_level+1);
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        plotmf_analytic[ilev].define(grids[ilev],dmap[ilev],3,0);
        MultiFab::Copy(plotmf_analytic[ilev], phi[ilev]            , 0, 0, 1, 0);
        MultiFab::Copy(plotmf_analytic[ilev], grad_x_analytic[ilev], 0, 1, 1, 0);
        MultiFab::Copy(plotmf_analytic[ilev], grad_y_analytic[ilev], 0, 2, 1, 0);
    }
    WriteMultiLevelPlotfile(plot_file_name + "-analytic", max_level+1,
                            amrex::GetVecOfConstPtrs(plotmf_analytic),
                            {"phi",
                            "dphidx",
                            "dphidy",
                             },
                            geom, 0.0, Vector<int>(max_level+1,0),
                            Vector<IntVect>(max_level,IntVect{2}));
#else


    plotmf[ilev].define(grids[ilev],dmap[ilev],4,0);
    plotmf[ilev].setVal(0.);
        MultiFab::Copy(plotmf[ilev], phi[ilev]   , 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], grad_x[ilev], 0, 1, 1, 0);
        MultiFab::Copy(plotmf[ilev], grad_y[ilev], 0, 2, 1, 0);
        MultiFab::Copy(plotmf[ilev], grad_z[ilev], 0, 3, 1, 0);
    }
    WriteMultiLevelPlotfile(plot_file_name, max_level+1,
                            amrex::GetVecOfConstPtrs(plotmf),
                            {"phi",
                            "dphidx",
                            "dphidy",
                            "dphidz"
                             },
                            geom, 0.0, Vector<int>(max_level+1,0),
                            Vector<IntVect>(max_level,IntVect{2}));

    Vector<MultiFab> plotmf_analytic(max_level+1);
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        plotmf_analytic[ilev].define(grids[ilev],dmap[ilev],4,0);
        MultiFab::Copy(plotmf_analytic[ilev], phi[ilev]            , 0, 0, 1, 0);
        MultiFab::Copy(plotmf_analytic[ilev], grad_x_analytic[ilev], 0, 1, 1, 0);
        MultiFab::Copy(plotmf_analytic[ilev], grad_y_analytic[ilev], 0, 2, 1, 0);
        MultiFab::Copy(plotmf_analytic[ilev], grad_z_analytic[ilev], 0, 3, 1, 0);
    }
    WriteMultiLevelPlotfile(plot_file_name + "-analytic", max_level+1,
                            amrex::GetVecOfConstPtrs(plotmf_analytic),
                            {"phi",
                            "dphidx",
                            "dphidy",
                            "dphidz"
                             },
                            geom, 0.0, Vector<int>(max_level+1,0),
                            Vector<IntVect>(max_level,IntVect{2}));

#endif

}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.queryarr("is_periodic", is_periodic);

    pp.query("eb_is_dirichlet", eb_is_dirichlet);
    pp.query("eb_is_homog_dirichlet", eb_is_homog_dirichlet);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(eb_is_dirichlet + eb_is_homog_dirichlet < 2,
       "Cannot set both Homogeneous Dirichlet and Dirichlet EB!!!");

    pp.query("plot_file", plot_file_name);

    pp.queryarr("prob_lo", prob_lo);
    pp.queryarr("prob_hi", prob_hi);

    scalars.resize(2);
    if (is_periodic[0]) {
        scalars[0] = 0.0;
        scalars[1] = 1.0;
    }
    else if (is_periodic[1]) {
        scalars[0] = 1.0;
        scalars[1] = 0.0;
    }
    else {
        scalars[0] = 1.0;
        scalars[1] = 1.0;
    }
    pp.queryarr("scalars", scalars);

    pp.query("use_linear_1d", use_linear_1d);
    pp.query("linear_1d_askew", linear_1d_askew);
    pp.queryarr("linear_1d_pt_on_top_wall",linear_1d_pt_on_top_wall);
    pp.query("linear_1d_height",linear_1d_height);
    pp.query("linear_1d_rotation",linear_1d_rotation);
    pp.queryarr("linear_1d_askew_rotation",linear_1d_askew_rotation);
    pp.query("linear_1d_flow_dir", linear_1d_flow_dir);
    pp.query("linear_1d_height_dir", linear_1d_height_dir);
    pp.query("linear_1d_bottom", linear_1d_bottom);
    pp.query("linear_1d_no_flow_dir", linear_1d_no_flow_dir);
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])}, {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(is_periodic[0],is_periodic[1],is_periodic[2])};
    Geometry::Setup(&rb, 0, isperiodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio);
    }
}
