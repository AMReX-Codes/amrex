#include <SMC.H>
#include <SMC_F.H>

using namespace amrex;

void SMC::compute_dU0 (FArrayBox& UU_fab, FArrayBox& Uprime_fab, FArrayBox& Q_fab, FArrayBox& mu_fab, FArrayBox& xi_fab, FArrayBox& lam_fab, FArrayBox& Ddiag_fab, Box bx, Box tbx0, Box tbxg)
{
    ctoprim_3d(tbxg.loVect(), tbxg.hiVect(), bx.loVect(), bx.hiVect(),
       	    UU_fab.dataPtr(), Q_fab.dataPtr(), ngrow, ngrow);
    chemterm_3d(tbx0.loVect(), tbx0.hiVect(),
	    BL_TO_FORTRAN_3D(     Q_fab),
	    BL_TO_FORTRAN_3D(Uprime_fab));
    get_trans_prop_3d(tbxg.loVect(), tbxg.hiVect(), bx.loVect(), bx.hiVect(),
	    Q_fab.dataPtr(), mu_fab.dataPtr(), xi_fab.dataPtr(), 
	    lam_fab.dataPtr(), Ddiag_fab.dataPtr(), ngrow);
}

void SMC::compute_dU1 (FArrayBox& UU_fab, FArrayBox& Uprime_fab, FArrayBox& Q_fab, FArrayBox& mu_fab, FArrayBox& xi_fab, FArrayBox& lam_fab, FArrayBox& Ddiag_fab, Box tbx)
{
    hypterm_3d(tbx.loVect(), tbx.hiVect(), dx,
	    BL_TO_FORTRAN_3D(    UU_fab),
	    BL_TO_FORTRAN_3D(     Q_fab),
	    BL_TO_FORTRAN_3D(Uprime_fab));

    narrow_diffterm_3d(tbx.loVect(), tbx.hiVect(), dx,
	    BL_TO_FORTRAN_3D(     Q_fab),
	    BL_TO_FORTRAN_3D(Uprime_fab),
	    mu_fab.dataPtr(), xi_fab.dataPtr(), 
	    lam_fab.dataPtr(), Ddiag_fab.dataPtr());
}

// Ua = a*Ua + b*Ub + c*Up
    void
SMC::update_rk3 (Real a, FArrayBox& Ua_fab, Real b, const FArrayBox& Ub_fab, Real c, const FArrayBox& Uc_fab, int stage, Box bx)
{
    int ncomp = U.nComp();
    Ua_fab.mult(a, bx, 0, ncomp);
    Ua_fab.saxpy(b, Ub_fab, bx, bx, 0, 0, ncomp);
    Ua_fab.saxpy(c, Uc_fab, bx, bx, 0, 0, ncomp);
}

    void
SMC::reset_density (FArrayBox& UU_fab, Box bx, Box tbx, int ngrow)
{
    reset_rho_3d(tbx.loVect(), tbx.hiVect(), bx.loVect(), bx.hiVect(), UU_fab.dataPtr(), ngrow);
}

