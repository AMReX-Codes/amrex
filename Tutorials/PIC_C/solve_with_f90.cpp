#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>

void 
solve_with_f90(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, Array< PArray<MultiFab> >& grad_phi_edge, 
               const PArray<Geometry>& geom, int base_level, Real tol, Real abs_tol)
{
    int nlevs = rhs.size() - base_level;

    BCRec phys_bc;
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
	BoxLib::Abort("non periodic boundraies not supported here");
    }

    for (int n=0; n<BL_SPACEDIM; n++) {
	phys_bc.setLo(n, Interior);  // for periodic boundary
	phys_bc.setHi(n, Interior);
    }

    MacBndry bndry(rhs[base_level].boxArray(), 1, geom[base_level]);
    //
    // Set Dirichlet boundary condition for phi in phi grow cells, use to
    // initialize bndry.
    //
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

    MultiFab* phi_p[nlevs];
    MultiFab* rhs_p[nlevs];
    for (int lev = base_level; lev < rhs.size(); lev++)
    {
        phi_p[lev-base_level] = &phi[lev];
        rhs_p[lev-base_level] = &rhs[lev];
    }

    // Need to set the boundary values here so they can get copied into "bndry"  
    if (base_level == 0)
    {
        bndry.setBndryValues(phi[base_level], src_comp, dest_comp, num_comp, phys_bc);
    }
#if 0
    else
    {
        MultiFab CPhi;
        Real cur_time = LevelData[level].get_state_data(State_Type).curTime();
        GetCrsePhi(level,CPhi,cur_time);
        BoxArray crse_boxes = BoxArray(grids[level]).coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes,in_rad,out_rad,extent_rad,num_comp);
        crse_br.copyFrom(CPhi,CPhi.nGrow(),src_comp,dest_comp,num_comp);

        bndry.setBndryValues(crse_br,src_comp,*(phi_p[base_level]),src_comp,
                             dest_comp,num_comp,crse_ratio,phys_bc);
    }
#endif

    std::vector<BoxArray>            bav  (nlevs);
    std::vector<DistributionMapping> dmv  (nlevs);
    std::vector<Geometry>            fgeom(nlevs);

    for (int lev = base_level; lev < rhs.size(); lev++)
    {
        bav[lev-base_level] = phi[lev].boxArray();
        dmv[lev-base_level] = phi[lev].DistributionMap();
        fgeom[lev-base_level] = geom[lev];
    }

    Array< Array<Real> > xa(nlevs);
    Array< Array<Real> > xb(nlevs);

    for (int lev = 0; lev < nlevs; ++lev)
    {
	xa[lev].resize(BL_SPACEDIM);
	xb[lev].resize(BL_SPACEDIM);

	if (lev + base_level == 0)
	{
	    for (int i = 0; i < BL_SPACEDIM; ++i)
	    {
		xa[lev][i] = 0;
		xb[lev][i] = 0;
	    }
	}
	else
	{
	    const Real* dx_crse = geom[lev+base_level-1].CellSize();
	    for (int i = 0; i < BL_SPACEDIM; ++i)
	    {
		xa[lev][i] = 0.5 * dx_crse[i];
		xb[lev][i] = 0.5 * dx_crse[i];
	    }
	}
    }

    {
        int stencil_type = CC_CROSS_STENCIL;
        int always_use_bnorm = 0;
        int need_grad_phi = 1;
        Real final_resnorm;
        MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type, false, 0, 1, 2);
        mgt_solver.set_const_gravity_coeffs(xa, xb);
        mgt_solver.solve(phi_p, rhs_p, bndry, tol, abs_tol, always_use_bnorm, 
			 final_resnorm,need_grad_phi); 

        for (int lev = base_level; lev < rhs.size(); lev++)
        {
            const Real* dx = geom[lev].CellSize();
            mgt_solver.get_fluxes(lev-base_level,grad_phi_edge[lev],dx);
        }
    }
}
