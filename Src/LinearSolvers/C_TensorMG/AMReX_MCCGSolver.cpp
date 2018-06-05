#include <algorithm>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_MCCGSolver.H>

namespace amrex {

namespace
{
    bool initialized = false;
}
//
// Set default values for these in Initialize()!!!
//
int    MCCGSolver::def_maxiter;
int    MCCGSolver::def_verbose;
int    MCCGSolver::def_isExpert;
double MCCGSolver::def_unstable_criterion;

void
MCCGSolver::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    MCCGSolver::def_maxiter            = 40;
    MCCGSolver::def_verbose            = 0;
    MCCGSolver::def_isExpert           = 0;
    MCCGSolver::def_unstable_criterion = 10;

    ParmParse pp("cg");

    pp.query("maxiter",  def_maxiter);
    pp.query("v",        def_verbose);
    pp.query("isExpert", def_isExpert);
    pp.query("unstable_criterion", def_unstable_criterion);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        amrex::OutStream() << "def_maxiter            = " << def_maxiter            << '\n'
                           << "def_unstable_criterion = " << def_unstable_criterion << '\n'
                           << "def_isExpert           = " << def_isExpert           << '\n';
    }

    amrex::ExecOnFinalize(MCCGSolver::Finalize);

    initialized = true;
}

void
MCCGSolver::Finalize ()
{
    initialized = false;
}

MCCGSolver::MCCGSolver (MCLinOp& _lp,
			bool     _use_mg_precond,
			int      _lev)
    :
    mg_precond(NULL),
    use_mg_precond(_use_mg_precond),
    isExpert((int)def_isExpert),
    Lp(_lp),
    lev(_lev)
{
    Initialize();
    maxiter = def_maxiter;
    verbose = def_verbose;
    set_mg_precond();
}

void
MCCGSolver::set_mg_precond ()
{
    delete mg_precond;
    if (use_mg_precond)
    {
	mg_precond = new MCMultiGrid(Lp);
	mg_precond->setNumIter(1);
    }
}

MCCGSolver::~MCCGSolver ()
{
    delete mg_precond;
}

Real
MCCGSolver::norm (const MultiFab& res)
{
    //
    // Compute max-norm.
    //
    const int p = 0, ncomp = res.nComp();

    Real restot = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:restot)
#endif
    for (MFIter mfi(res,true); mfi.isValid(); ++mfi)
    {
	restot = std::max(restot, res[mfi].norm(mfi.tilebox(), p, 0, ncomp));
    }
    ParallelDescriptor::ReduceRealMax(restot);
    return restot;
}

void
MCCGSolver::solve (MultiFab&       sol,
		   const MultiFab& rhs,
		   Real            eps_rel,
		   Real            eps_abs,
		   MCBC_Mode       bc_mode)
{
    //
    // algorithm:
    //   k=0;r=rhs-A*soln_0;
    //   while (||r_k||^2_2 > eps^2*||r_o||^2_2 && k < maxiter {
    //      k++
    //      solve Mz_k-1 = r_k-1 (if preconditioning, else z_k-1 = r_k-1)
    //      rho_k-1 = r_k-1^T z_k-1
    //      if (k=1) { p_1 = z_0 }
    //      else { beta = rho_k-1/rho_k-2; p = z + beta*p }
    //      w = Ap
    //      alpha = rho_k-1/p^tw
    //      x += alpha p
    //      r -= alpha w
    //   }
    //
    BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
    BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));

    int nghost = 1, ncomp  = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();

    MultiFab s(ba, dm, ncomp, nghost);
    MultiFab r(ba, dm, ncomp, nghost);
    MultiFab z(ba, dm, ncomp, nghost);
    MultiFab w(ba, dm, ncomp, nghost);
    MultiFab p(ba, dm, ncomp, nghost);
    //
    // Copy initial guess into a temp multifab guaranteed to have ghost cells.
    //
    int srccomp=0;  int destcomp=0;  nghost=0;
    s.copy(sol,srccomp,destcomp,ncomp);

    /* Note:
	 This routine assumes the MCLinOp is linear, and that when bc_mode =
	 MCHomogeneous_BC, MCLinOp::apply() on a zero vector will return a zero
	 vector.  Given that, we define the problem we solve here from the
	 original equation:

	      Lp(sol) = rhs --> Lp(s) + Lp(sol,bc_mode=MCHomogeneous_BC) = rhs

	 where s is set to the incoming solution guess.  Rewriting,

	      Lp(sol,bc_mode=MCHomogeneous_BC) = r     [ = rhs - Lp(s) ].

	 CG needs the residual of this equation on our initial guess.  But
	 because we made the above assumption,

	      r - Lp(sol,bc_mode=MCHomogeneous_BC) = r = rhs - Lp(s)

	 Which is simply the residual of the original equation evaluated at
	 the initial guess.  Thus we get by with only one call to Lp.residual.
	 Without this assumption, we'd need two.
         */
    Lp.residual(r, rhs, s, lev, bc_mode);
    //
    // Set initial guess for correction to 0.
    //
    sol.setVal(0);
    //
    // Set bc_mode=homogeneous.
    //
    MCBC_Mode temp_bc_mode=MCHomogeneous_BC;
    Real rnorm  = norm(r);
    Real rnorm0 = rnorm;
    Real minrnorm = rnorm;
    int ret = 0; // will return this value if all goes well

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        for (int k = 0; k < lev; k++)
            amrex::OutStream() << "   ";
        amrex::OutStream() << "MCCGsolver: Initial error (error0) =  " << rnorm0 << '\n';
    }

    Real beta = 0, rho = 0, rhoold = 0;

    /* WARNING:
	 The MultiFab copies used below to update z and p require nghost=0
	 to avoid the possibility of filling valid regions with uninitialized
	 data in the invalid regions of neighboring grids.  The default
	 behavior in MultiFab copies will likely be changed in the future.
	 */

    //
    // Note: if eps_rel or eps_abs < 0: that test is effectively bypassed.
    //
    for (int nit = 0;
         (nit < maxiter) && (rnorm > eps_rel*rnorm0) && (rnorm > eps_abs);
         ++nit)
    {
	if (use_mg_precond)
	{
            //
	    // solve Mz_k-1 = r_k-1  and  rho_k-1 = r_k-1^T z_k-1
            //
	    z.setVal(0);
	    mg_precond->solve( z, r, eps_rel, eps_abs, temp_bc_mode );
	}
        else
        {
            //
	    // No preconditioner, z_k-1 = r_k-1  and  rho_k-1 = r_k-1^T r_k-1.
            //
	    srccomp=0;  destcomp=0;  
	    z.copy(r, srccomp, destcomp, ncomp);
	}

	rho = MultiFab::Dot(r, 0, z, 0, ncomp, 0);
	
	if (nit == 0)
	{
            //
	    // k=1, p_1 = z_0.
            //
	    srccomp=0;  destcomp=0;  nghost=0;
	    p.copy(z, srccomp, destcomp, ncomp);
	}
        else
        {
            //
	    // k>1, beta = rho_k-1/rho_k-2 and  p = z + beta*p
            //
	    beta = rho/rhoold;
	    advance( p, beta, z );
	}
        //
	// w = Ap, and compute Transpose(p).w
        //
	Real pw = axp( w, p, temp_bc_mode );
	//
	// alpha = rho_k-1/p^tw.
        //
	Real alpha = rho/pw;
	
	if (verbose > 2 && ParallelDescriptor::IOProcessor())
        {
            for (int k = 0; k < lev; k++)
                amrex::OutStream() << "   ";
            amrex::OutStream() << "MCCGSolver:"
                               << " nit " << nit
                               << " pw "  << pw 
                               << " rho " << rho
                               << " alpha " << alpha;
            if (nit == 0)
                amrex::OutStream() << " beta undefined ...";
            else
                amrex::OutStream() << " beta " << beta << " ...";
	}
        //
	// x += alpha p  and  r -= alpha w
        //
	rhoold = rho;
	update( sol, alpha, r, p, w );
	rnorm = norm(r);
        if (rnorm > def_unstable_criterion*minrnorm)
        {
            ret = 2;
            break;
        }
        else if (rnorm < minrnorm)
        {
            minrnorm = rnorm;
        }

	if (verbose > 1 ||
            (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
              (eps_abs > 0. && rnorm < eps_abs)) && verbose))
	{
	    if (ParallelDescriptor::IOProcessor())
	    {
		for (int k = 0; k < lev; k++) {
                    amrex::OutStream() << "   ";
                }
		amrex::OutStream() << "MCCGSolver: Iteration "
                                   << nit
                                   << " error/error0 "
                                   << rnorm/rnorm0 << '\n';
	    }
	}
    }

    if (ret != 0 && isExpert == false)
    {
        amrex::Error("MCCGSolver:: apparent accuracy problem; try expert setting or change unstable_criterion");
    }
    if (ret==0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        amrex::Error("MCCGSolver:: failed to converge!");
    }
    //
    // Omit ghost update since maybe not initialized in calling routine.
    //
    if (ret == 0)
    {
        srccomp=0; nghost=0;
        sol.plus(s,srccomp,ncomp,nghost);
    }
}

void
MCCGSolver::advance (MultiFab&       p,
		     Real            beta,
		     const MultiFab& z)
{
    //
    // Compute p = z  +  beta p
    //
    int ncomp = p.nComp();
    int nghost = 0;
    MultiFab::Xpay(p, beta, z, 0, 0, ncomp, nghost);
}

void
MCCGSolver::update (MultiFab&       sol,
		    Real            alpha,
		    MultiFab&       r,
		    const MultiFab& p,
		    const MultiFab& w)
{
    //
    // Compute sol =+ alpha p  and  r -= alpha w
    //
    int ncomp = r.nComp();
    MultiFab::Saxpy(sol,  alpha, p, 0, 0, ncomp, 0);
    MultiFab::Saxpy(r  , -alpha, w, 0, 0, ncomp, 0);
}

Real
MCCGSolver::axp (MultiFab& w,
		 MultiFab& p,
		 MCBC_Mode bc_mode)
{
    //
    // Compute w = A.p, and return Transpose(p).w
    //
    Lp.apply(w, p, lev, bc_mode);
    int nghost = 0;
    return MultiFab::Dot(w,0,p,0,p.nComp(),nghost);
}

}
