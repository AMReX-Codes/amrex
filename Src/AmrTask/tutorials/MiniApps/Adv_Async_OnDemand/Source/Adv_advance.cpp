#include <Adv.H>
#include <Adv_F.H>
#include <WorkerThread.H>
#include <RegionGraph.H>
#include <Perilla.H>
#include <string>
#include <fstream>

using namespace perilla;

    Real
Adv::advance (Real time,
	Real dt,
	int  iteration,
	int  ncycle)
{
    if(isMasterThread())
    {
	for (int k = 0; k < NUM_STATE_TYPE; k++) {
	    state[k].allocOldData();
	    state[k].swapTimeLevels(dt);
	}
    }

    MultiFab& S_new = get_new_data(State_Type);

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time = state[State_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    // Get pointers to Flux registers, or set pointer to zero if not there.
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;

    int finest_level = parent->finestLevel();

    if(isMasterThread())
    {
	if (do_reflux && level < finest_level) {
	    fine = &getFluxReg(level+1);
	    fine->setVal(0.0);
	}
    }

    if (do_reflux && level > 0) {
	current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux)
    {
	for (int j = 0; j < BL_SPACEDIM; j++)
	{
	    BoxArray ba = S_new.boxArray();
	    ba.surroundingNodes(j);
	    fluxes[j].define(ba, dmap, NUM_STATE, 0);
	}
    }

    amrex::Vector<AsyncFillPatchIterator*> upperAFPI;
    if(level < parent->finestLevel())
    {
	Adv& upperLevel = getLevel(level+1);
	upperAFPI = upperLevel.SborderFPI;
    }
    FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

#pragma omp parallel
{
    for (RGIter rgi(SborderFPI, upperAFPI, *(Sborder), NUM_GROW, prev_time, State_Type, 0, NUM_STATE, iteration); rgi.isValid(); ++rgi){
	int f = rgi.currentRegion;
	int fid = S_new.IndexArray()[f];
	int fis = Sborder->IndexArray()[f];
	const FArrayBox& statein = (*(Sborder))[fis];
	FArrayBox& stateout      =   S_new[fid];
	MFIter mfi(S_new, true);


	const Box& bx = rgi.tilebox();
	rgi.sync_workers();
	for (int i = 0; i < BL_SPACEDIM ; i++) {
	    const Box& bxtmp = amrex::surroundingNodes(bx,i);
	    {
		flux[i].resize(bxtmp,NUM_STATE);
		uface[i].resize(amrex::grow(bxtmp,1),1);
	    }
	}

	get_face_velocity(level, ctr_time,
		D_DECL(BL_TO_FORTRAN(uface[0]),
		    BL_TO_FORTRAN(uface[1]),
		    BL_TO_FORTRAN(uface[2])),
		dx, prob_lo);
	advect(time, bx.loVect(), bx.hiVect(),
		BL_TO_FORTRAN_3D(statein),
		BL_TO_FORTRAN_3D(stateout),
		D_DECL(BL_TO_FORTRAN_3D(uface[0]),
		    BL_TO_FORTRAN_3D(uface[1]),
		    BL_TO_FORTRAN_3D(uface[2])),
		D_DECL(BL_TO_FORTRAN_3D(flux[0]),
		    BL_TO_FORTRAN_3D(flux[1]),
		    BL_TO_FORTRAN_3D(flux[2])),
		dx, dt);

	if (do_reflux) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		fluxes[i][f].copy(flux[i],mfi.nodaltilebox(i));       
	}
    }
}
    if (do_reflux) {
	if (fine) {
	    if(isMasterThread())
		for (int i = 0; i < BL_SPACEDIM ; i++)
		    fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
	}
    }
    return dt;
}
