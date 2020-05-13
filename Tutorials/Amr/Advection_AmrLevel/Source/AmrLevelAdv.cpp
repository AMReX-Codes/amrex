
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9;
int      AmrLevelAdv::do_reflux       = 1;

int      AmrLevelAdv::NUM_STATE       = 1;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 3;  // number of ghost cells

#ifdef AMREX_PARTICLES
std::unique_ptr<AmrTracerParticleContainer> AmrLevelAdv::TracerPC =  nullptr;
int AmrLevelAdv::do_tracers                       =  0;
#endif

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
    flux_reg = 0;
}

//
//The basic constructor.
//
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
void
AmrLevelAdv::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers and level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

//
//Write a plotfile to specified directory.
//
void
AmrLevelAdv::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)
{

    AmrLevel::writePlotFile (dir,os,how);

#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      TracerPC->Checkpoint(dir, "Tracer", true);
    }
#endif
}

//
//Define data descriptors.
//
void
AmrLevelAdv::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
			   &cell_cons_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	lo_bc[i] = hi_bc[i] = BCType::int_dir;   // periodic boundaries
    }
    
    BCRec bc(lo_bc, hi_bc);

    desc_lst.setComponent(Phi_Type, 0, "phi", bc, 
			  StateDescriptor::BndryFunc(nullfill));
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelAdv::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

          initdata(&level, &cur_time, AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
		   BL_TO_FORTRAN_3D(S_new[mfi]), AMREX_ZFILL(dx),
		   AMREX_ZFILL(prob_lo));
    }

#ifdef AMREX_PARTICLES
    init_particles();
#endif

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                       << " data " << std::endl;
    }
}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
//
void
AmrLevelAdv::init (AmrLevel &old)
{
    AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{
    MultiFab& S_mm = get_new_data(Phi_Type);
    Real maxval = S_mm.max(0);
    Real minval = S_mm.min(0);

    amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
	fine = &getFluxReg(level+1);
	fine->setVal(0.0);
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

    // State with ghost cells
    MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

    // MF to hold the mac velocity
    MultiFab Umac[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      BoxArray ba = S_new.boxArray();
      ba.surroundingNodes(i);
      Umac[i].define(ba, dmap, 1, iteration);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    const FArrayBox& statein = Sborder[mfi];
	    FArrayBox& stateout      =   S_new[mfi];

	    // Allocate fabs for fluxes and Godunov velocities.
	    for (int i = 0; i < BL_SPACEDIM ; i++) {
		const Box& bxtmp = amrex::surroundingNodes(bx,i);
		flux[i].resize(bxtmp,NUM_STATE);
		uface[i].resize(amrex::grow(bxtmp, iteration), 1);
	    }

	    get_face_velocity(&level, &ctr_time,
			      AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
				     BL_TO_FORTRAN(uface[1]),
				     BL_TO_FORTRAN(uface[2])),
			      dx, prob_lo);

	    for (int i = 0; i < BL_SPACEDIM ; i++) {
                const Box& bxtmp = mfi.grownnodaltilebox(i, iteration);
                Umac[i][mfi].copy<RunOn::Host>(uface[i], bxtmp);
	    }
            advect(&time, bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN_3D(statein), 
		   BL_TO_FORTRAN_3D(stateout),
		   AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
			  BL_TO_FORTRAN_3D(uface[1]),
			  BL_TO_FORTRAN_3D(uface[2])),
		   AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
			  BL_TO_FORTRAN_3D(flux[1]), 
			  BL_TO_FORTRAN_3D(flux[2])), 
		   dx, &dt);

	    if (do_reflux) {
		for (int i = 0; i < BL_SPACEDIM ; i++)
		    fluxes[i][mfi].copy<RunOn::Host>(flux[i],mfi.nodaltilebox(i));
	    }
	}
    }

    if (do_reflux) {
	if (current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
	}
	if (fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
	}
    }

#ifdef AMREX_PARTICLES
    if (TracerPC) {
      TracerPC->AdvectWithUmac(Umac, level, dt);
    }
#endif

    return dt;
}

//
//Estimate time step.
//
Real
AmrLevelAdv::estTimeStep (Real)
{
    // This is just a dummy value to start with 
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[Phi_Type].curTime();
    const MultiFab& S_new = get_new_data(Phi_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	FArrayBox uface[BL_SPACEDIM];

	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    for (int i = 0; i < BL_SPACEDIM ; i++) {
		const Box& bx = mfi.nodaltilebox(i);
		uface[i].resize(bx,1);
	    }

	    get_face_velocity(&level, &cur_time,
			      AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
				     BL_TO_FORTRAN(uface[1]),
				     BL_TO_FORTRAN(uface[2])),
			      dx, prob_lo);

	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		Real umax = uface[i].norm<RunOn::Host>(0);
		if (umax > 1.e-100) {
		    dt_est = std::min(dt_est, dx[i] / umax);
		}
	    }
	}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    if (verbose) {
	amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
                       << ":  dt_est = " << dt_est << std::endl;
    }
    
    return dt_est;
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        AmrLevelAdv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) 
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else 
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Do work after timestep().
//
void
AmrLevelAdv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

#ifdef AMREX_PARTICLES    
    if (TracerPC)
      {
        const int ncycle = parent->nCycle(level);
	
        if (iteration < ncycle || level == 0)
	  {
            int ngrow = (level == 0) ? 0 : iteration;

	    TracerPC->Redistribute(level, TracerPC->finestLevel(), ngrow);
	  }
      }
#endif
}

//
//Do work after regrid().
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest) {
#ifdef AMREX_PARTICLES
  if (TracerPC && level == lbase) {
      TracerPC->Redistribute(lbase);
  }
#endif
}

//
//Do work after a restart().
//
void
AmrLevelAdv::post_restart() 
{
#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      BL_ASSERT(TracerPC == 0);
      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->Restart(parent->theRestartFile(), "Tracer");
    }
#endif
}

//
//Do work after init().
//
void
AmrLevelAdv::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;
	
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box&  tilebx  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    
	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
	    tagfab.get_itags(itags, tilebx);
	    
            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebx.loVect();
	    const int*  thi     = tilebx.hiVect();

	    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
			BL_TO_FORTRAN_3D(S_new[mfi]),
			&tagval, &clearval, 
			AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()), 
			AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);
	    //
	    // Now update the tags in the TagBox.
	    //
	    tagfab.tags_and_untags(itags, tilebx);
	}
    }
}

void
AmrLevelAdv::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");   

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // This tutorial code only supports Cartesian coordinates.
    if (! gg->IsCartesian()) {
	amrex::Abort("Please set geom.coord_sys = 0");
    }

    // This tutorial code only supports periodic boundaries.
    if (! gg->isAllPeriodic()) {
	amrex::Abort("Please set geom.is_periodic = 1 1 1");
    }

#ifdef AMREX_PARTICLES
    pp.query("do_tracers", do_tracers);
#endif 

    //
    // read tagging parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

    // use a fortran routine to
    // read in tagging parameters from probin file
    get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

void
AmrLevelAdv::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
    
    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = amrex::second() - strt;
	
        ParallelDescriptor::ReduceRealMax(end,IOProc);
	
        amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
                       << " : time = " << end << std::endl;
    }
}

void
AmrLevelAdv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelAdv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

#ifdef AMREX_PARTICLES
void
AmrLevelAdv::init_particles ()
{
  if (do_tracers and level == 0)
    {
      BL_ASSERT(TracerPC == nullptr);
      
      TracerPC.reset(new AmrTracerParticleContainer(parent));

      AmrTracerParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
    }
}
#endif
