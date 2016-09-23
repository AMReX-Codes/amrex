#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MGT_Solver.H>
#include <stencil_types.H>

void 
solve_with_f90(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi, 
               const Array<Geometry>& geom, int base_level, int finest_level, Real tol, Real abs_tol)
{
    int nlevs = finest_level - base_level + 1;

    int mg_bc[2*BL_SPACEDIM];

    // This tells the solver that we are using periodic bc's
    if (Geometry::isAllPeriodic())
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir)
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
    } else {
	BoxLib::Abort("non periodic boundaries not supported here");
    }

    // Have to do some packing because these arrays does not always start with base_level
    PArray<Geometry> geom_p(nlevs);
    PArray<MultiFab> rhs_p(nlevs);
    PArray<MultiFab> phi_p(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
	geom_p.set(ilev, &geom[ilev+base_level]);
	rhs_p.set (ilev,  &rhs[ilev+base_level]);
	phi_p.set (ilev,  &phi[ilev+base_level]);
    }
    
    // Refinement ratio is hardwired to 2 here.
    IntVect crse_ratio = (base_level == 0) ? 
	IntVect::TheZeroVector() : IntVect::TheUnitVector() * 2;

    // This is the cell-centered BoxArray
    std::vector<BoxArray> bav(1);
    bav[0] = rhs[0].boxArray();

    std::vector<DistributionMapping> dmv(1);
    dmv[0] = rhs[0].DistributionMap();

    bool nodal     = true;
    bool have_rhcc = false;
    int stencil = ND_CROSS_STENCIL;

    int verbose = 1;

    MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal, stencil, have_rhcc, 0, 1, verbose);

    mgt_solver.set_nodal_const_coefficients(1);

#if 0
//  mgt_solver.nodal_project(&phi[c_lev], &vel[c_lev], &rhs_cc[c_lev], rhnd,
//                           rel_tol, abs_tol, &lo_inflow[0], &hi_inflow[0]);


    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    fmg.solve(phi_p, rhs_p, tol, abs_tol, always_use_bnorm, need_grad_phi);
   
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
	int amr_level = ilev + base_level;
	fmg.get_fluxes(grad_phi[amr_level], ilev);
    }
#endif
}
