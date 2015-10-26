// This solves (alpha * a - beta * (del dot b grad)) phi = rhs using MGT_Sover,
// which in turn uses Fortran multigrid solver F_MG.  Here a, phi and rhs are 
// cell-centered arrays (i.e., MultiFab), alpha and beta is a scalar, and b is
// face-centered array represented by BL_SPACEDIM face-centered MultiFabs.
//
// Steps to use this solver:
//
// Step 1: 
//   Build FMultiGrid with Geometry.  The baselevel and crse_ratio arguments are 
//   optional.  They are needed when the problem to solve has a base level greater
//   than AMR level 0.  That is the lowest local level in the solver is at level
//   baselevel in the AMR hierachy.  Also note that geom[0] is on the local level 0.
//   In the case of baselevel != 0, crse_ratio must be set to the refinement ratio
//   between AMR levels baselevel-1 and baselevel.
//
// Step 2: OPTIONAL
//   Set stencil type and verbosity.
//
// Step 3:
//   Use one of the set_bc functions to pass in multigrid boundary conditions
//   and set boundary values.  Here the valuse in mg_bc should match the ones in
//   F_MG/mg_cpp_f.h.  For periodic boundaries, it can set to either BC_PER 
//   (i.e., -1) or BC_INT (i.e., 0).
//
// Step 4:
//   * For a constant gravity problem, call set_const_gravity_coeffs.
//   * For a curvilinear gravity problem, call set_gravity_coeffs.
//   * For a general problem, call set_alpha, set_acoef, set_beta, and set_bcoef.
//
// Step 5:
//   Call solve.  This will store the solution in phi.
//
// Step 6: OPTIONAL
//   Call get_grad_phi to get grad_phi.  The need_grad_phi flag in last step
//   must be 1.  Note that grad_phi is PArray with BL_SPACEDIM face-cented
//   MultiFabs or Array<PArray> in case of multiple levels.
//
// If one only wants to test a solution by compute residual without actually
// solving the problem.  Replace Step 5 with a call to compute_residual.
//
// Various MultiFabs are passed in, and the pointers to these MultiFabs,
// not copies, are stored.  Therefore, these MultiFabs must still be alive
// when solve is called.

#ifndef _FMultiGrid_H_
#define _FMultiGrid_H_

#include <MGT_Solver.H>
#include <stencil_types.H>
#include <MacBndry.H>
#include <mg_cpp_f.h>

class FMultiGrid
{
public:

    // single level
    FMultiGrid (const             Geometry  & geom,
		int                           baselevel = 0,
		IntVect                       crse_ratio = IntVect::TheZeroVector());

    // multi level
    FMultiGrid (const std::vector<Geometry> & geom,
		int                           baselevel = 0,
		IntVect                       crse_ratio = IntVect::TheZeroVector());

    // multi level
    FMultiGrid (const      PArray<Geometry> & geom,
		int                           baselevel = 0,
		IntVect                       crse_ratio = IntVect::TheZeroVector());

    ~FMultiGrid () { delete m_bndry;   delete m_mgt_solver; }

    void set_baselevel (int baselevel, IntVect& crse_ratio) {
	m_baselevel = baselevel;
	m_crse_ratio = crse_ratio;
    }
    void set_stencil (int stencil) {m_stencil = stencil;}
    void set_verbose (int verbose) {m_verbose = verbose;}

    // Set homogeneous Dirichlet boundary
    void set_bc (int * mg_bc);

    // Set boundry values using phi
    void set_bc (int     * mg_bc,
		 MultiFab& phi);

    // Set boundary values with crse level data, crse_phi
    void set_bc (int     * mg_bc,
		 MultiFab& crse_phi,
		 MultiFab& phi);

    // For const_gravity_eq
    void set_const_gravity_coeffs ();

    // For gravity_eq
    void set_gravity_coeffs (       PArray<MultiFab>  & b); // single level
    void set_gravity_coeffs (Array< PArray<MultiFab> >& b); // multi level
    
    // For general_eq
    void set_alpha (Real alpha);
    //
    void set_acoef (       MultiFab & a);
    void set_acoef (PArray<MultiFab>& a);
    //
    void set_beta (Real beta);    
    //
    void set_bcoef (      PArray<MultiFab>   & b);
    void set_bcoef (Array<PArray<MultiFab> > & b);

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

    void get_grad_phi (       PArray<MultiFab>  & grad_phi);
    void get_grad_phi ( Array<PArray<MultiFab> >& grad_phi);
    void get_grad_phi (PArray<PArray<MultiFab> >& grad_phi);

    void compute_residual (MultiFab & phi,
			   MultiFab & rhs,
			   MultiFab & res);
    void compute_residual (PArray<MultiFab> & phi,
			   PArray<MultiFab> & rhs,
			   PArray<MultiFab> & res);

private:

    int     m_nlevels;
    int     m_baselevel;
    IntVect m_crse_ratio;
    int     m_stencil;
    int     m_verbose;

    std::vector<Geometry> m_geom;

    MacBndry   * m_bndry;
    MGT_Solver * m_mgt_solver;

    struct Boundary
    {
	bool initilized;
	Array<int> mg_bc;
	MultiFab* crse_phi;
	MultiFab* phi;
	Boundary () : initilized(false) {}
	Boundary (const int * bc, MultiFab* cmf = 0, MultiFab * mf = 0) 
	    : initilized(true),
	      mg_bc(bc, 2*BL_SPACEDIM),
	      crse_phi(cmf),
	      phi(mf)
	    { }
	int * get_mg_bc () { return mg_bc.dataPtr(); }
	void set_bndry_values (MacBndry& bndry, IntVect crse_ratio);
    };
    Boundary m_bc;

    enum Eq_t { invalid_eq,
		const_gravity_eq, // del dot grad phi = rhs  
		gravity_eq,       // (del dot b grad) phi = rhs  
		general_eq};

    struct ABecCoeff
    {
	Eq_t eq_type;
	bool alpha_set;
	bool a_set;
	bool beta_set;
	bool b_set;
	Real alpha;
	PArray<MultiFab> a;
	Real beta;
	Array<PArray<MultiFab> > b;
	ABecCoeff () : eq_type(invalid_eq), alpha_set(false), a_set(false), 
		       beta_set(false), b_set(false), alpha(0), beta(0) {}
	void set_coeffs (MGT_Solver & mgt_solver, FMultiGrid& fmg);
    };
    ABecCoeff m_coeff;

    static void Copy (PArray<MultiFab>& dst, MultiFab& src);
    static void Copy (PArray<MultiFab>& dst, PArray<MultiFab>& src);
    static void Copy (Array<PArray<MultiFab> >& dst, PArray<MultiFab>& src);
    static void Copy (Array<PArray<MultiFab> >& dst, Array<PArray<MultiFab> >& src);

    void init_mgt_solver (PArray<MultiFab>& phi);
};

#endif
