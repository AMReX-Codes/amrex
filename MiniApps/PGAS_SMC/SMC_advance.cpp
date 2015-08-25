#include <SMC.H>
#include <SMC_F.H>

void
SMC::advance (int istep)
{
    compute_dUdt(U, istep);
    Utmp.setVal(0.0);
    update_rk3(0.0, Utmp, 1.0, U, dt, Uprime);
    reset_density(Utmp);

    compute_dUdt(Utmp);
    update_rk3(0.25, Utmp, 0.75, U, 0.25*dt, Uprime);
    reset_density(Utmp);

    compute_dUdt(Utmp);
    update_rk3(1./3., U, 2./3., Utmp, (2./3.)*dt, Uprime);
    reset_density(U);
}

void 
SMC::compute_dUdt (MultiFab& UU, int istep)
{
    Uprime.setVal(0.0);
}

// Ua = a*Ua + b*Ub + c*Up
void
SMC::update_rk3 (Real a, MultiFab& Ua, Real b, const MultiFab& Ub, Real c, const MultiFab& Uc)
{
    int ncomp = Ua.nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Ua,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();
	Ua[mfi].mult(a, bx, 0, ncomp);
	Ua[mfi].saxpy(b, Ub[mfi], bx, bx, 0, 0, ncomp);
	Ua[mfi].saxpy(c, Uc[mfi], bx, bx, 0, 0, ncomp);
    }
}

void
SMC::reset_density (MultiFab& UU)
{
    int ngrow = UU.nGrow();
    for (MFIter mfi(UU); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.validbox();
	reset_rho_3d(bx.loVect(), bx.hiVect(), &ngrow, UU[mfi].dataPtr());
    }
}

