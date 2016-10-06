
#include <WarpX.H>
#include <PICSAR_f.H>

namespace {
    constexpr Real clight = 299792458.;
    constexpr Real mu0 = 1.2566370614359173e-06;
};

void
WarpX::Evolve ()
{
    BL_PROFILE("WPX::Evolve");

    Real t = 0.0;
    Real dt  = 0.95 * geom_arr[0].CellSize(0) * (1.0/clight);

    for (int istep = 0; istep < max_step; ++istep, t += dt)
    {
	// At the beginning, we have B^{n-1/2} and E^{n}.
	// Particles have p^{n-1/2} and x^{n}.

	EvolveBfield(0.5*dt); // We now B^{n}
    

//	DepositCharge();

//	DepositCurrent(dt);
	
//	EvolveEfield(dt);
    }
}

void
WarpX::DepositCharge ()
{
    BL_PROFILE("WPX::DepositCharge");

    {
	BL_PROFILE_VAR("BoxLib::DepositCharge", bldep);
	int lev = 0;
	int ncomp = 1;
	mypc->NodalDepositionSingleLevel(charge,lev,ncomp);
    }

    {
	BL_PROFILE_VAR("PICSAR::DepositCharge", picsardep);
	int lev = 0;
	mypc->ChargeDeposition(charge, lev, charge_deposition_order);
    }
}

void
WarpX::DepositCurrent (Real dt)
{
    BL_PROFILE("WPX::DepositCurrent");
    int lev = 0;
    mypc->CurrentDeposition(current,lev,dt); 
}

void
WarpX::EvolveBfield (Real dt)
{
    BL_PROFILE("WPX::EvolveBfield");

    Real mu_c2_dt = (mu0*clight*clight) * dt;

    const Real* dx = geom_arr[0].CellSize();

    Real dtsdx[3];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	dtsdx[i] = dt / dx[i];
    }

    long norder = 2;
    long nguard = 1;
    long nstart = 0;
    int l_nodal = 0;

    BL_ASSERT(nguard == Efield[0]->nGrow());
    BL_ASSERT(nguard == Efield[1]->nGrow());
    BL_ASSERT(nguard == Efield[2]->nGrow());
    BL_ASSERT(nguard == Bfield[0]->nGrow());
    BL_ASSERT(nguard == Bfield[1]->nGrow());
    BL_ASSERT(nguard == Bfield[2]->nGrow());

    for ( MFIter mfi(*Bfield[0]); mfi.isValid(); ++mfi )
    {
	const Box& box = BoxLib::enclosedCells(mfi.validbox());
	long nx = box.length(0);
	long ny = box.length(1);
	long nz = box.length(2); 

	pxrpush_em3d_bvec_norder( (*Efield[0])[mfi].dataPtr(),
				  (*Efield[1])[mfi].dataPtr(),
				  (*Efield[2])[mfi].dataPtr(),
				  (*Bfield[0])[mfi].dataPtr(),
				  (*Bfield[1])[mfi].dataPtr(),
				  (*Bfield[2])[mfi].dataPtr(), 
				  dtsdx, dtsdx+1, dtsdx+2,
				  &nx, &ny, &nz,
				  &norder, &norder, &norder,
				  &nguard, &nguard, &nguard,
				  &nstart, &nstart, &nstart,
				  &l_nodal );
    }
}

void
WarpX::EvolveEfield (Real dt)
{
    BL_PROFILE("WPX::EvolveEfield");

    Real mu_c2_dt = (mu0*clight*clight) * dt;

    const Real* dx = geom_arr[0].CellSize();

    Real dtsdx_c2[3];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	dtsdx_c2[i] = (dt*clight*clight) / dx[i];
    }

    long norder = 2;
    long nguard = 1;
    long nstart = 0;
    int l_nodal = 0;

    BL_ASSERT(nguard == Efield[0]->nGrow());
    BL_ASSERT(nguard == Efield[1]->nGrow());
    BL_ASSERT(nguard == Efield[2]->nGrow());
    BL_ASSERT(nguard == Bfield[0]->nGrow());
    BL_ASSERT(nguard == Bfield[1]->nGrow());
    BL_ASSERT(nguard == Bfield[2]->nGrow());
    BL_ASSERT(nguard == current[0]->nGrow());
    BL_ASSERT(nguard == current[1]->nGrow());
    BL_ASSERT(nguard == current[2]->nGrow());

    for ( MFIter mfi(*Efield[0]); mfi.isValid(); ++mfi )
    {
	const Box & bx = BoxLib::enclosedCells(mfi.validbox());
	long nx = bx.length(0);
	long ny = bx.length(1);
	long nz = bx.length(2); 

	pxrpush_em3d_evec_norder( (*Efield[0])[mfi].dataPtr(),
				  (*Efield[1])[mfi].dataPtr(),
				  (*Efield[2])[mfi].dataPtr(),
				  (*Bfield[0])[mfi].dataPtr(),
				  (*Bfield[1])[mfi].dataPtr(),
				  (*Bfield[2])[mfi].dataPtr(), 
				  (*current[0])[mfi].dataPtr(),
				  (*current[1])[mfi].dataPtr(),
				  (*current[2])[mfi].dataPtr(),
				  &mu_c2_dt, dtsdx_c2, dtsdx_c2+1, dtsdx_c2+2,
				  &nx, &ny, &nz,
				  &norder, &norder, &norder,
				  &nguard, &nguard, &nguard,
				  &nstart, &nstart, &nstart,
				  &l_nodal );
    }
}
