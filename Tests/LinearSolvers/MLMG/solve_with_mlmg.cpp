
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>

#include <prob_par.H>

using namespace amrex;

namespace my {
    extern int max_iter;
    extern int verbose;
}

void solve_with_mlmg (const Vector<Geometry>& geom,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      Vector<MultiFab>& rhs)
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

    MLABecLaplacian mlabec(geom, grids, dmap);
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

        if (ilev == 0) {
            mlabec.setDirichletBC(0, *psoln[0]);
        } else {
            mlabec.setDirichletBC(ilev, *psoln[ilev], psoln[ilev-1]);            
        }
    }
    
    MLMG mlmg(mlabec);
    mlmg.setMaxIter(my::max_iter);
    mlmg.setVerbose(my::verbose);
    mlmg.solve(psoln, prhs, tol_rel, tol_abs);
}

