
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

#include <prob_par.H>

using namespace amrex;

void solve_with_mlmg (const Vector<Geometry>& geom,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      const Vector<MultiFab>& rhs)
{
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

    MLABecLaplacian mlabc(geom, grids, dmap);
    // set bc ???
    mlabc.setScalars(prob::a, prob::b);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        mlabc.setACoeffs(ilev, alpha[ilev]);
//        mlabs.setBCoeffs();  // for each direction
    }
    
//    MLMG mlmg(mlabc);
//    mlmg.solve(psoln, prhs, tol_rel, tol_abs);

}

