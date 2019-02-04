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

// void 
// solve_with_f90(const Vector<MultiFab*>& rhs,
// 	       const Vector<MultiFab*>& phi, 
//                const Vector< Vector<MultiFab*> >& grad_phi_edge, 
//                const Vector<Geometry>& geom, int base_level, int finest_level,
// 	       Real tol, Real abs_tol)
// {
//     int nlevs = finest_level - base_level + 1;

//     int mg_bc[2*BL_SPACEDIM];

//     // This tells the solver that we are using periodic bc's
//     if (Geometry::isAllPeriodic())
//     {
//         for (int dir = 0; dir < BL_SPACEDIM; ++dir)
//         {
//             mg_bc[2*dir + 0] = 0;
//             mg_bc[2*dir + 1] = 0;
//         }
//     } else {
// 	amrex::Abort("non periodic boundaries not supported here");
//     }

//     // Have to do some packing because these arrays does not always start with base_level
//     Vector<Geometry> geom_p(nlevs);
//     Vector<MultiFab*> rhs_p(nlevs);
//     Vector<MultiFab*> phi_p(nlevs);
//     for (int ilev = 0; ilev < nlevs; ++ilev) {
// 	geom_p[ilev] = geom[ilev+base_level];
// 	rhs_p[ilev]  = rhs[ilev+base_level];
// 	phi_p[ilev]  = phi[ilev+base_level];
//     }
    
//     // Refinement ratio is hardwired to 2 here.
//     IntVect crse_ratio = (base_level == 0) ? 
// 	IntVect::TheZeroVector() : IntVect::TheUnitVector() * 2;

//     FMultiGrid fmg(geom_p, base_level, crse_ratio);

//     if (base_level == 0) {
// 	fmg.set_bc(mg_bc, *phi[base_level]);
//     } else {
// 	fmg.set_bc(mg_bc, *phi[base_level-1], *phi[base_level]);
//     }

//     fmg.set_const_gravity_coeffs();

//     int always_use_bnorm = 0;
//     int need_grad_phi = 1;
//     fmg.set_verbose(0);
//     fmg.solve(phi_p, rhs_p, tol, abs_tol, always_use_bnorm, need_grad_phi);
   
//     for (int ilev = 0; ilev < nlevs; ++ilev)
//     {
// 	int amr_level = ilev + base_level;
// 	fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
//     }
// }
