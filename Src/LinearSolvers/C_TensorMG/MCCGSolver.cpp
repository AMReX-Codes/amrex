
//
// $Id: MCCGSolver.cpp,v 1.14 2004-01-07 21:18:43 car Exp $
//
#include <winstd.H>

#include <algorithm>

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <CG_F.H>
#include <MCCGSolver.H>

int MCCGSolver::initialized = 0;
int MCCGSolver::def_maxiter = 40;
int MCCGSolver::def_verbose = 0;
int MCCGSolver::def_isExpert = 0;
double MCCGSolver::def_unstable_criterion = 10.;

void
MCCGSolver::initialize ()
{
    ParmParse pp("cg");
    pp.query("maxiter", def_maxiter);
    pp.query("v", def_verbose);
    pp.query("isExpert", def_isExpert);
    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
	std::cout << "def_maxiter=" << def_maxiter << '\n';
	std::cout << "def_unstable_criterion=" << def_unstable_criterion << '\n';
        std::cout << "def_isExpert=" << def_isExpert << '\n';
    }
    initialized = 1;
}

MCCGSolver::MCCGSolver (MCLinOp& _Lp,
			bool     _use_mg_precond,
			int      _lev)
    :
    isExpert((int)def_isExpert),
    Lp(_Lp),
    lev(_lev),
    use_mg_precond(_use_mg_precond),
    mg_precond(NULL)
{
    if (!initialized)
	initialize();
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
    const int p = 0;
    Real restot = 0.0;
    Real resk   = 0.0;
    const int ncomp = res.nComp();
    const BoxArray& gbox = res.boxArray();
    for (MFIter mfi(res); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(mfi.validbox() == gbox[mfi.index()]);
        resk = res[mfi].norm(mfi.validbox(), p, 0, ncomp);
        if (p == 0)
	{
	    restot = std::max(restot, resk);
        }
        else if (p == 2)
        {
	    restot += resk*resk;
        }
        else
	    BoxLib::Error("BOGUS P IN NORM");
    }
    if (p == 0)
    {
	ParallelDescriptor::ReduceRealMax(restot);
    }
    else if (p == 2)
    {
	ParallelDescriptor::ReduceRealSum(restot);
	restot = std::sqrt(restot);
    }
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

    int nghost = 1;
    int ncomp  = sol.nComp();
    MultiFab* s = new MultiFab(sol.boxArray(), ncomp, nghost);
    MultiFab* r = new MultiFab(sol.boxArray(), ncomp, nghost);
    MultiFab* z = new MultiFab(sol.boxArray(), ncomp, nghost);
    MultiFab* w = new MultiFab(sol.boxArray(), ncomp, nghost);
    MultiFab* p = new MultiFab(sol.boxArray(), ncomp, nghost);
    //
    // Copy initial guess into a temp multifab guaranteed to have ghost cells.
    //
    int srccomp=0;  int destcomp=0;  nghost=0;
    s->copy(sol,srccomp,destcomp,ncomp);

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
    Lp.residual((*r),  rhs, (*s), lev, bc_mode);
    //
    // Set initial guess for correction to 0.
    //
    sol.setVal(0.0);
    //
    // Set bc_mode=homogeneous.
    //
    MCBC_Mode temp_bc_mode=MCHomogeneous_BC;
    Real rnorm  = norm(*r);
    Real rnorm0 = rnorm;
    Real minrnorm = rnorm;
    int ret = 0; // will return this value if all goes well

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        for (int k = 0; k < lev; k++)
            std::cout << "   ";
        std::cout << "MCCGsolver: Initial error (error0) =  " << rnorm0 << '\n';
    }

    Real beta, rho, rhoold=0.0;

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
	    z->setVal(0.);
	    mg_precond->solve( *z, *r, eps_rel, eps_abs, temp_bc_mode );
	}
        else
        {
            //
	    // No preconditioner, z_k-1 = r_k-1  and  rho_k-1 = r_k-1^T r_k-1.
            //
	    srccomp=0;  destcomp=0;  
	    z->copy(*r, srccomp, destcomp, ncomp);
	}

	rho = 0.0;
	int ncomp = z->nComp();
	const BoxArray& gbox = r->boxArray();

        for (MFIter rmfi(*r); rmfi.isValid(); ++rmfi)
	{
            Real trho;
	    FORT_CGXDOTY(
		&trho,
		(*z)[rmfi].dataPtr(), 
                ARLIM((*z)[rmfi].loVect()), ARLIM((*z)[rmfi].hiVect()),
		(*r)[rmfi].dataPtr(), 
                ARLIM((*r)[rmfi].loVect()), ARLIM((*r)[rmfi].hiVect()),
		rmfi.validbox().loVect(), rmfi.validbox().hiVect(), &ncomp);
	    rho += trho;
	}
        ParallelDescriptor::ReduceRealSum(rho);
	
	if (nit == 0)
	{
            //
	    // k=1, p_1 = z_0.
            //
	    srccomp=0;  destcomp=0;  nghost=0;
	    p->copy(*z, srccomp, destcomp, ncomp);
	}
        else
        {
            //
	    // k>1, beta = rho_k-1/rho_k-2 and  p = z + beta*p
            //
	    beta = rho/rhoold;
	    advance( (*p), beta, (*z) );
	}
        //
	// w = Ap, and compute Transpose(p).w
        //
	Real pw = axp( (*w), (*p), temp_bc_mode );
	//
	// alpha = rho_k-1/p^tw.
        //
	Real alpha = rho/pw;
	
	if (verbose > 2 && ParallelDescriptor::IOProcessor())
        {
            for (int k = 0; k < lev; k++)
                std::cout << "   ";
            std::cout << "MCCGSolver:"
                      << " nit " << nit
                      << " pw "  << pw 
                      << " rho " << rho
                      << " alpha " << alpha;
            if (nit == 0)
                std::cout << " beta undefined ...";
            else
                std::cout << " beta " << beta << " ...";
	}
        //
	// x += alpha p  and  r -= alpha w
        //
	rhoold = rho;
	update( sol, alpha, (*r), (*p), (*w) );
	rnorm = norm(*r);
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
		for (int k = 0; k < lev; k++)
                    std::cout << "   ";
		std::cout << "MCCGSolver: Iteration "
                          << nit
                          << " error/error0 "
                          << rnorm/rnorm0 << '\n';
	    }
	}
    }

    if (ret != 0 && isExpert == false)
    {
        BoxLib::Error("MCCGSolver:: apparent accuracy problem; try expert setting or change unstable_criterion");
    }
    if (ret==0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        BoxLib::Error("MCCGSolver:: failed to converge!");
    }
    //
    // Omit ghost update since maybe not initialized in calling routine.
    // BoxLib_1.99 has no MultiFab::plus(MultiFab&) member, which would
    // operate only in valid regions; do explicitly.  Add to boundary
    // values stored in initialsolution.
    //
    if (ret == 0)
    {
        srccomp=0; nghost=0;
        sol.plus((*s),srccomp,ncomp,nghost);
    }

    delete s;
    delete r;
    delete w;
    delete p;
    delete z;
}

void
MCCGSolver::advance (MultiFab&       p,
		     Real            beta,
		     const MultiFab& z)
{
    //
    // Compute p = z  +  beta p
    //
    int             ncomp = p.nComp();
    const BoxArray& gbox  = Lp.boxArray(lev);
    const BoxArray& zbox  = z.boxArray();

    for (MFIter pmfi(p); pmfi.isValid(); ++pmfi)
    {
        const Box& bx = zbox[pmfi.index()];

        BL_ASSERT(bx == gbox[pmfi.index()]);

	FORT_CGADVCP(
	    p[pmfi].dataPtr(),
            ARLIM(p[pmfi].loVect()), ARLIM(p[pmfi].hiVect()),
	    z[pmfi].dataPtr(),
            ARLIM(z[pmfi].loVect()), ARLIM(z[pmfi].hiVect()),
	    &beta, bx.loVect(), bx.hiVect(),&ncomp);
    }
}

void
MCCGSolver::update (MultiFab&       sol,
		    Real            alpha,
		    MultiFab&       r,
		    const MultiFab& p,
		    const MultiFab& w)
{
    //
    // Compute x =+ alpha p  and  r -= alpha w
    //
    const BoxArray& gbox = Lp.boxArray(lev);
    int ncomp = r.nComp();
    for (MFIter solmfi(sol); solmfi.isValid(); ++solmfi)
    {
	FORT_CGUPDATE(
	    sol[solmfi].dataPtr(),
            ARLIM(sol[solmfi].loVect()),ARLIM(sol[solmfi].hiVect()),
	    r[solmfi].dataPtr(),
            ARLIM(r[solmfi].loVect()),ARLIM(r[solmfi].hiVect()),
	    &alpha,
	    w[solmfi].dataPtr(),
            ARLIM(w[solmfi].loVect()),ARLIM(w[solmfi].hiVect()),
	    p[solmfi].dataPtr(),
            ARLIM(p[solmfi].loVect()), ARLIM(p[solmfi].hiVect()),
	    solmfi.validbox().loVect(), solmfi.validbox().hiVect(),&ncomp);
    }
}

Real
MCCGSolver::axp (MultiFab& w,
		 MultiFab& p,
		 MCBC_Mode bc_mode)
{
    //
    // Compute w = A.p, and return Transpose(p).w
    //
    Real pw = 0.0;
    const BoxArray& gbox = Lp.boxArray(lev);
    Lp.apply(w, p, lev, bc_mode);
    int ncomp = p.nComp();
    for (MFIter pmfi(p); pmfi.isValid(); ++pmfi)
    {
	Real tpw;
	FORT_CGXDOTY(
	    &tpw,
	    p[pmfi].dataPtr(),ARLIM(p[pmfi].loVect()),ARLIM(p[pmfi].hiVect()),
	    w[pmfi].dataPtr(),ARLIM(w[pmfi].loVect()),ARLIM(w[pmfi].hiVect()),
	    pmfi.validbox().loVect(), pmfi.validbox().hiVect(),&ncomp);
	pw += tpw;
    }
    ParallelDescriptor::ReduceRealSum(pw);
    return pw;
}
