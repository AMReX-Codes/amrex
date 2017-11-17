
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <prob_par.H>

using namespace amrex;

namespace {
    static int max_iter = 100;
    static int max_fmg_iter = 20;
    static int verbose  = 2;
    static int cg_verbose = 0;
    static int linop_maxorder = 2;
}

void solve_with_mlmg (const Vector<Geometry>& geom,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      Vector<MultiFab>& rhs)
{
    BL_PROFILE("solve_with_mlmg");

    {
        ParmParse pp;
        pp.query("max_iter", max_iter);
        pp.query("max_fmg_iter", max_fmg_iter);
        pp.query("verbose", verbose);
        pp.query("cg_verbose", cg_verbose);
        pp.query("linop_maxorder", linop_maxorder);
    }

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    Vector<BoxArray> grids;
    Vector<DistributionMapping> dmap;

    Vector<MultiFab*> psoln;
    Vector<MultiFab const*> prhs;

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids.push_back(soln[ilev].boxArray());
        dmap.push_back(soln[ilev].DistributionMap());
        psoln.push_back(&(soln[ilev]));
        prhs.push_back(&(rhs[ilev]));
    }

    MLABecLaplacian mlabec(geom, grids, dmap);

    mlabec.setMaxOrder(linop_maxorder);

    // BC
    mlabec.setDomainBC({prob::bc_type,prob::bc_type,prob::bc_type},
                       {prob::bc_type,prob::bc_type,prob::bc_type});
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        mlabec.setLevelBC(ilev, psoln[ilev]);
    }

    mlabec.setScalars(prob::a, prob::b);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        mlabec.setACoeffs(ilev, alpha[ilev]);

        std::array<MultiFab,AMREX_SPACEDIM> bcoefs;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(beta[ilev].boxArray(), IntVect::TheDimensionVector(idim));
            bcoefs[idim].define(ba, beta[ilev].DistributionMap(), 1, 0);
        }
        amrex::average_cellcenter_to_face({AMREX_D_DECL(&bcoefs[0],
                                                        &bcoefs[1],
                                                        &bcoefs[2])},
                                          beta[ilev], geom[ilev]);
        mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcoefs));
    }
    
    MLMG mlmg(mlabec);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setCGVerbose(cg_verbose);

    mlmg.solve(psoln, prhs, tol_rel, tol_abs);
}

