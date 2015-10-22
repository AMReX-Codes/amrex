// This solves (a - beta * (del dot b grad)) phi = rhs using MGT_Sover, which
// in turn uses Fortran multigrid solver F_MG.  Here a, phi and rhs are 
// cell-centered arrays (i.e., MultiFab), beta is a scalar, and b is face-centered
// array represented by BL_SPACEDIM face-centered MultiFabs.
//
// There are three special cases:
//   (1) (alpha - beta * (del dot b grad)) phi = rhs, where alpha is a scalar;
//   (2) - beta * (del dot b grad) phi = rhs;
//   (3) - del dot grad phi = rhs.
//
// Steps to use this solver:
//
// Step 1: 
//   Build FMultiGrid with Geometry, DistributionMapping, and BoxArray.  The
//   baselevel and crse_ratio arguments are optional.  They are needed when
//   the problem to solve has a base level greater than AMR level 0.  That is
//   the lowest local level in the solver is at level baselevel in the AMR 
//   hierachy.  Also note that geom[0] , dmap[0], and ba[0] are on the local
//   level 0.  In the case of baselevel != 0, crse_ratio must be set to the
//   refinement ratio between AMR levels baselevel-1 and baselevel.
//
// Step 2: OPTIONAL
//   Set the number of component, stencil type and verbosity.
//
// Step 3:
//   Use one of the set_bc functions to pass in multigrid boundary conditions
//   and set boundary values.  Here the valuse in mg_bc should, for periodic,
//   be 0 and, for others, match the ones in F_MG/mg_cpp_f.h. 
//
// Step 4:
//   Use one of the set_*_coeffs functions to set various coeffs.
//
// Step 5:
//   Call solve.
//
// Step 6: OPTIONAL
//   Call get_grad_phi to get grad_phi.  The need_grad_phi flag in last step
//   must be 1.  Note that grad_phi is PArray with BL_SPACEDIM face-cented
//   MultiFabs.
//
// If one only wants to test a solution by compute residual without actually
// solving the problem.  Replace Step 5 with a call to compute_residual.

#ifndef _FMultiGrid_H_
#define _FMultiGrid_H_

#include <MGT_Solver.H>
#include <stencil_types.H>
#include <MacBndry.H>

class FMultiGrid
{
public:
    
    FMultiGrid (const Geometry            & geom,
		const DistributionMapping & dmap,
		const BoxArray            & ba,
		int                         baselevel = 0,
		int                         crse_ratio = 0);

    FMultiGrid (const std::vector<Geometry>           & geom,
		const std::vector<DistributionMapping>& dmap,
		const std::vector<BoxArray>           & ba,
		int                                     baselevel = 0,
		int                                     crse_ratio = 0);

    FMultiGrid (const PArray<Geometry>           & geom,
		const PArray<DistributionMapping>& dmap,
		const PArray<BoxArray>           & ba,
		int                                baselevel = 0,
		int                                crse_ratio = 0);

    ~FMultiGrid ();

    void set_baselevel (int baselevel, int crse_ratio) {
	m_baselevel = baselevel;
	m_crse_ratio = crse_ratio;
	if (m_baselevel > 0 && m_crse_ratio < 2) 
	    BoxLib::Abort("FMultiGrid::set_baselevel must set crse_ratio if baselevel > 0");
    }
    void set_ncomp   (int ncomp)   {m_ncomp   = ncomp;}
    void set_stencil (int stencil) {m_stencil = stencil;}
    void set_verbose (int verbose) {m_verbose = verbose;}

    // Set homogeneous Dirichlet boundary
    void set_bc (int * mg_bc);

    // Set boundry values using phi
    void set_bc (int           * mg_bc,
		 const MultiFab& phi);

    // Set boundary values with crse level data, crse_phi
    void set_bc (int           * mg_bc,
		 const MultiFab& crse_phi,
		 const MultiFab& phi);

    void set_abeclap_coeffs (Real alpha, // A special case of a is const
			     Real beta,
			     const Array<PArray<MultiFab> >& b);

    void set_abeclap_coeffs (const PArray<MultiFab>& a,
			     Real beta,
			     const Array<PArray<MultiFab> >& b);

    void set_gravity_coeffs (const Array< PArray<MultiFab> >& b);

    void set_const_gravity_coeffs ();

    Real solve (MultiFab& phi,
		MultiFab& rhs,
		Real rel_tol, Real abs_tol,
		int always_use_bnorm = 0,
		int need_grad_phi = 0,
		int verbose = 0);

    Real solve (PArray<MultiFab>& phi,
		PArray<MultiFab>& rhs,
		Real rel_tol, Real abs_tol,
		int always_use_bnorm = 0,
		int need_grad_phi = 0,
		int verbose = 0);

    void get_grad_phi (int lev, PArray<MultiFab>& grad_phi);

    void compute_residual (MultiFab & phi,
			   MultiFab & rhs,
			   MultiFab & res);

    void compute_residual (PArray<MultiFab> & phi,
			   PArray<MultiFab> & rhs,
			   PArray<MultiFab> & res);

private:

    void build_mgt_solver ();
    void make_xa_xb (Array < Array<Real> > & xa,
		     Array < Array<Real> > & xb);

    int m_nlevels;
    int m_baselevel;
    int m_crse_ratio;
    int m_ncomp;
    int m_stencil;
    int m_verbose;

    PArray<Geometry>            m_geom;
    PArray<DistributionMapping> m_dmap;
    PArray<BoxArray>            m_ba;

    int m_coeffs_set;
    int m_bc_set;

    MacBndry* m_bndry;
    int m_mg_bc[2*BL_SPACEDIM];

    MGT_Solver* m_mgt_solver;
};

#endif
