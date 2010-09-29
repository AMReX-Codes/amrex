#include <iostream>

#include "hg_multi.H"

#if defined( BL_FORT_USE_UNDERSCORE )
#define   FORT_HGCG1            hgcg1_
#define   FORT_HGCG2            hgcg2_
#define   FORT_HGIP             hgip_
#elif defined( BL_FORT_USE_UPPERCASE )
#define   FORT_HGCG1            HGCG1
#define   FORT_HGCG2            HGCG2
#define   FORT_HGIP             HGIP
#elif defined( BL_FORT_USE_LOWERCASE )
#define   FORT_HGCG1            hgcg1
#define   FORT_HGCG2            hgcg2
#define   FORT_HGIP             hgip
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C"
{
    void FORT_HGCG1         (Real*, const Real*, Real*, Real*,
			     const Real*, const Real*,
			     const Real*, intS, const Real&, Real&);
    void FORT_HGCG2         (Real*, const Real*, intS, const Real&);
    void FORT_HGIP          (const Real*, const Real*,
			     const Real*, intS, Real&);
}

void
holy_grail_amr_multigrid::cgsolve (int mglev)
{
    BL_ASSERT(mglev == 0);

    MultiFab& r = cgwork[0];
    MultiFab& p = cgwork[1];
    MultiFab& z = cgwork[2];
    MultiFab& x = cgwork[3];
    MultiFab& w = cgwork[4];
    MultiFab& c = cgwork[5];
    MultiFab& zero_array = cgwork[6];
    MultiFab& ipmask = cgwork[7];
    //
    // x (corr[0]) should be all 0.0 at this point
    //
    for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
    {
	r[r_mfi].copy(resid[mglev][r_mfi]);
	r[r_mfi].negate();
    }

    if (singular)
    {
        //
	// Singular systems are very sensitive to solvability
        //
	w.setVal(1.0);
	Real aa = inner_product(r, w) / mg_domain[mglev].volume();
	r.plus(-aa, 0);
    }

    Real rho = 0.0;
    for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
    {
	z[r_mfi].copy(r[r_mfi]);
	z[r_mfi].mult(c[r_mfi]);
	const Box& reg = p[r_mfi].box();
	FORT_HGIP(z[r_mfi].dataPtr(), r[r_mfi].dataPtr(),
                  ipmask[r_mfi].dataPtr(), DIMLIST(reg), rho);
	p[r_mfi].copy(z[r_mfi]);
    }
    ParallelDescriptor::ReduceRealSum(rho);
    if ( pcode >= 3 && ParallelDescriptor::IOProcessor() )
    {
	std::cout << "      HG: cgsolve rho = " << rho << std::endl;
    }

    const Real tol = HG::cgsolve_tolfact * rho;

    int i = 0;
    while (tol > 0.0)
    {
	if ( ++i > HG::cgsolve_maxiter )
	  {
	    if ( ParallelDescriptor::IOProcessor() )
	      {
		BoxLib::Warning( "cgsolve: Conjugate-gradient iteration failed to converge" );
	      }
	    break;
	}
	Real rho_old = rho;
        //
	// safe to set the clear flag to 0 here---bogus values make it
	// into r but are cleared from z by the mask in c
        //
	level_residual(w, zero_array, p, 0, false, 0);
	Real alpha = 0.0;
	for (MFIter p_mfi(p); p_mfi.isValid(); ++p_mfi)
	{
	    const Box& reg = p[p_mfi].box();
	    FORT_HGIP(p[p_mfi].dataPtr(), w[p_mfi].dataPtr(),
                      ipmask[p_mfi].dataPtr(), DIMLIST(reg), alpha);
	}
	ParallelDescriptor::ReduceRealSum(alpha);
	alpha = rho / alpha;
	rho = 0.0;
	for (MFIter r_mfi(r); r_mfi.isValid(); ++r_mfi)
	{
	    const Box& reg = p[r_mfi].box();
	    FORT_HGCG1(r[r_mfi].dataPtr(),
                       p[r_mfi].dataPtr(),
                       z[r_mfi].dataPtr(),
		       x[r_mfi].dataPtr(),
                       w[r_mfi].dataPtr(),
                       c[r_mfi].dataPtr(),
		       ipmask[r_mfi].dataPtr(),
		       DIMLIST(reg), alpha, rho);
	}
	ParallelDescriptor::ReduceRealSum(rho);
	if (pcode >= 3  && ParallelDescriptor::IOProcessor())
	{
	    std::cout << "      HG: cgsolve iter(" << i << ") rho=" << rho << std::endl;
	}
	if (rho <= tol)
	    break;
	alpha = rho / rho_old;
	for (MFIter p_mfi(p); p_mfi.isValid(); ++p_mfi)
	{
	    const Box& reg = p[p_mfi].box();
	    FORT_HGCG2(p[p_mfi].dataPtr(),z[p_mfi].dataPtr(),DIMLIST(reg),alpha);
	}
    }

    if (pcode >= 3  && ParallelDescriptor::IOProcessor())
    {
	std::cout << "      HG: "
                  << i << " iterations required for conjugate-gradient" << std::endl;
    }
}
