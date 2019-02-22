#include <iostream>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>

using namespace amrex;

void 
solve_with_mlmg(const Vector<MultiFab*>& rhs,
                const Vector<MultiFab*>& phi, 
                const Vector< Vector<MultiFab*> >& grad_phi_edge, 
                const Vector<Geometry>& geom, int base_level, int finest_level,
                Real tol, Real abs_tol)
{
    constexpr int max_coarsening_level = 30;
    constexpr bool agglomeration = true;
    constexpr bool consolidation = true;
    constexpr int linop_maxorder = 3;
    constexpr int max_iter = 100;
    constexpr int max_fmg_iter = 100;
    constexpr int verbose = 2;
    constexpr int bottom_verbose = 0;
    
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    int nlevels = finest_level - base_level + 1;
    
    // Have to do some packing because these arrays does not always start with base_level
    Vector<Geometry> geom_p(nlevels);
    Vector<const MultiFab*> rhs_p(nlevels);
    Vector<MultiFab*> phi_p(nlevels);
    Vector<BoxArray> grids(nlevels);
    Vector<DistributionMapping> dmap(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
	geom_p[ilev] = geom[ilev+base_level];
	rhs_p[ilev]  = rhs[ilev+base_level];
	phi_p[ilev]  = phi[ilev+base_level];
        grids[ilev]  = rhs[ilev+base_level]->boxArray();
        dmap[ilev]   = rhs[ilev+base_level]->DistributionMap();
    }
    
    MLPoisson mlpoisson(geom_p, grids, dmap, info);

    mlpoisson.setMaxOrder(linop_maxorder);

    // This is a 3d problem with Periodic BCs
    mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                        LinOpBCType::Periodic,
                                        LinOpBCType::Periodic)},
        {AMREX_D_DECL(LinOpBCType::Periodic,
                      LinOpBCType::Periodic,
                      LinOpBCType::Periodic)});
    
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        mlpoisson.setLevelBC(ilev, phi[ilev]);
    }

    MLMG mlmg(mlpoisson);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.solve(phi_p, rhs_p, tol, abs_tol);
}
