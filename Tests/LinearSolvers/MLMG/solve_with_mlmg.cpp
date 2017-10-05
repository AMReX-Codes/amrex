
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

#include <prob_par.H>

using namespace amrex;

void solve_with_mlmg (const Array<Geometry>& geom, int rr,
                      Array<MultiFab>& soln,
                      const Array<MultiFab>& alpha, const Array<MultiFab>& beta,
                      const Array<MultiFab>& rhs)
{
    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    Array<int> ref_ratio(nlevels-1,rr);

    Array<BoxArray> grids;
    Array<DistributionMapping> dmap;

    Array<MultiFab*> psoln;
    Array<MultiFab const*> prhs;

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids.push_back(soln[ilev].boxArray());
        dmap.push_back(soln[ilev].DistributionMap());
        psoln.push_back(&(soln[ilev]));
        prhs.push_back(&(rhs[ilev]));
    }

//    MLABecLaplacian mlabc(geom, grids, dmap, ref_ratio);
    // set bc ???
//    mlabc.setScalars(a, b);
//    for  each level
//      mlabc.setACoeffs();
//      mlabs.setBCoeffs();  // for each direction
    
//    MLMG mlmg(mlabc);
//    mlmg.solve(psoln, prhs, tol_rel, tol_abs);

}

