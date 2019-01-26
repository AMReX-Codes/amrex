
#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef FORCING
#include "Forcing.H"
#endif

using namespace amrex;

using std::string;

Real
Nyx::advance (Real time,
              Real dt,
              int  iteration,
              int  ncycle)

  // Arguments:
  //    time      : the current simulation time
  //    dt        : the timestep to advance (e.g., go from time to
  //                time + dt)
  //    iteration : where we are in the current AMR subcycle.  Each
  //                level will take a number of steps to reach the
  //                final time of the coarser level below it.  This
  //                counter starts at 1
  //    ncycle    : the number of subcycles at this level

{
#ifndef NO_HYDRO
    if (do_hydro)
    {
        if (Nyx::theActiveParticles().size() > 0)
        {
#ifndef AGN
           if (do_dm_particles)
#endif
           {
               return advance_hydro_plus_particles(time, dt, iteration, ncycle);
           } 
        }
        else
        {
           return advance_hydro(time, dt, iteration, ncycle);
        }
    }
#endif

#ifdef GRAVITY
    if ( ! do_hydro)
    {
        return advance_particles_only(time, dt, iteration, ncycle);
    }
#else
    else
    {
        amrex::Abort("Nyx::advance -- do_hydro is false but no gravity -- dont know what to do");
    }
#endif
    return 0;
}

//
// This will advance the Nyx AmrLevel. 
// If no subcycling is used, this will be a full multilevel advance at 
//      at level 0 and a finer levels will be skipped over.
// If subcycling is used, this will check if finer levels are subcycled 
//      relative to this level. All levels that are not subcycled will be 
//      advanced in a multilevel advance with this level and be skipped over
//      subsequently. If the next level is subcycled relative to this one, 
//      then this will only be a single level advance.
//
#ifndef NO_HYDRO
Real
Nyx::advance_hydro_plus_particles (Real time,
                                   Real dt,
                                   int  iteration,
                                   int  ncycle)
{
    // A particle in cell (i) can affect cell values in (i-1) to (i+1)
    int stencil_deposition_width = 1;
 
    // A particle in cell (i) may need information from cell values in (i-1) to (i+1)
    //   to update its position (typically via interpolation of the acceleration from the grid)
    int stencil_interpolation_width = 1;
 
    // A particle that starts in cell (i + ncycle) can reach
    //   cell (i) in ncycle number of steps .. after "iteration" steps
    //   the particle has to be within (i + ncycle+1-iteration) to reach cell (i)
    //   in the remaining (ncycle-iteration) steps
 
    // *** ghost_width ***  is used
    //   *) to set how many cells are used to hold ghost particles i.e copies of particles
    //      that live on (level-1) can affect the grid over all of the ncycle steps.
    //      We define ghost cells at the coarser level to cover all iterations so
    //      we can't reduce this number as iteration increases.
 
    int ghost_width = ncycle + stencil_deposition_width;

    // *** where_width ***  is used
    //   *) to set how many cells the Where call in moveKickDrift tests =
    //      ghost_width + (1-iteration) - 1:
    //      the minus 1 arises because this occurs *after* the move

    int where_width =  ghost_width + (1-iteration)  - 1;
 
    // *** grav_n_grow *** is used
    //   *) to determine how many ghost cells we need to fill in the MultiFab from
    //      which the particle interpolates its acceleration
    //   *) to set how many cells the Where call in moveKickDrift tests = (grav.nGrow()-2).
    //   *) the (1-iteration) arises because the ghost particles are created on the coarser
    //      level which means in iteration 2 the ghost particles may have moved 1 additional cell along
 
    int grav_n_grow = ghost_width + (1-iteration) + (iteration-1) +
                      stencil_interpolation_width ;

    BL_PROFILE_REGION_START("R::Nyx::advance_hydro_plus_particles");
    BL_PROFILE("Nyx::advance_hydro_plus_particles()");
    // Sanity checks
    if (!do_hydro)
        amrex::Abort("In `advance_hydro_plus_particles` but `do_hydro` not true");

   if (Nyx::theActiveParticles().size() <= 0)
        amrex::Abort("In `advance_hydro_plus_particles` but no active particles");

#ifdef GRAVITY
    if (!do_grav)
        amrex::Abort("In `advance_hydro_plus_particles` but `do_grav` not true");
#endif
#ifdef FORCING
    if (do_forcing)
        amrex::Abort("Forcing in `advance_hydro_plus_particles` not admissible");
#endif

    const int finest_level = parent->finestLevel();
    int finest_level_to_advance;
    bool nosub = !parent->subCycle();
    
    if (nosub)
    {
        if (level > 0)
            return dt;
            
        finest_level_to_advance = finest_level;
    }
    else
    {
        if (strict_subcycling)
        {
            finest_level_to_advance = level;
        }
        else
        {
            // This level was advanced by a previous multilevel advance.
            if (level > 0 && ncycle == 1)
                return dt;
            
            // Find the finest level to advance
            int lev = level;
            while(lev < finest_level && parent->nCycle(lev+1) == 1)
                lev++;
            finest_level_to_advance = lev;
        }

        // We must setup virtual and Ghost Particles
        //
        // Setup the virtual particles that represent finer level particles
        // 
        setup_virtual_particles();
        //
        // Setup ghost particles for use in finer levels. Note that Ghost particles
        // that will be used by this level have already been created, the
        // particles being set here are only used by finer levels.
        //
        for(int lev = level; lev <= finest_level_to_advance && lev < finest_level; lev++)
        {
           get_level(lev).setup_ghost_particles(ghost_width);
        }
    }

    Real dt_lev;

    //
    // Move current data to previous, clear current.
    // Don't do this if a coarser level has done this already.
    //
    if (level == 0 || iteration > 1)
    {
        for (int lev = level; lev <= finest_level; lev++)
        {
            dt_lev = parent->dtLevel(lev);
            for (int k = 0; k < NUM_STATE_TYPE; k++)
            {
                get_level(lev).state[k].allocOldData();
                get_level(lev).state[k].swapTimeLevels(dt_lev);
            }
        }
    }

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();
    const Real a_old     = get_comoving_a(prev_time);
    const Real a_new     = get_comoving_a(cur_time);

#ifdef GRAVITY
    //
    // We now do a multilevel solve for old Gravity. This goes to the 
    // finest level regardless of subcycling behavior. Consequentially,
    // If we are subcycling we skip this step on the first iteration of
    // finer levels.
    BL_PROFILE_VAR("solve_for_old_phi", solve_for_old_phi);
    if (level == 0 || iteration > 1)
    {
        // fix fluxes on finer grids
        if (do_reflux)
        {
            for (int lev = level; lev < finest_level; lev++)
            {
                gravity->zero_phi_flux_reg(lev + 1);
            }
        }
        
        // swap grav data
        for (int lev = level; lev <= finest_level; lev++)
            get_level(lev).gravity->swap_time_levels(lev);

        //
        // Solve for phi using the previous phi as a guess.
        //
        int use_previous_phi_as_guess = 1;
        int ngrow_for_solve = iteration + stencil_deposition_width;
        gravity->multilevel_solve_for_old_phi(level, finest_level,
                                              ngrow_for_solve,
                                              use_previous_phi_as_guess);
    }
    BL_PROFILE_VAR_STOP(solve_for_old_phi);
    //
    // Advance Particles
    //
    if (Nyx::theActiveParticles().size() > 0)
    {
        // Advance the particle velocities to the half-time and the positions to the new time
        // We use the cell-centered gravity to correctly interpolate onto particle locations
        if (particle_move_type == "Gravitational")
        {
            const Real a_half = 0.5 * (a_old + a_new);

            if (particle_verbose && ParallelDescriptor::IOProcessor())
                std::cout << "moveKickDrift ... updating particle positions and velocity\n"; 

            for (int lev = level; lev <= finest_level_to_advance; lev++)
            {
                // We need grav_n_grow grow cells to track boundary particles
                const auto& ba = get_level(lev).get_new_data(State_Type).boxArray();
                const auto& dm = get_level(lev).get_new_data(State_Type).DistributionMap();
                MultiFab grav_vec_old(ba, dm, BL_SPACEDIM, grav_n_grow);
                get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
                
                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                    Nyx::theActiveParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, where_width);

                // Only need the coarsest virtual particles here.
                if (lev == level && level < finest_level)
                    for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
                       Nyx::theVirtualParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, where_width);

                // Miiiight need all Ghosts
                for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
                   Nyx::theGhostParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half, where_width);
            }
        }
    }

#endif

    //
    // Call the hydro advance at each level to be advanced
    //
    BL_PROFILE_VAR("just_the_hydro", just_the_hydro);
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
#ifdef SDC
        if (sdc_split > 0) 
        { 
           get_level(lev).sdc_hydro(time, dt, a_old, a_new);
        } else { 
           get_level(lev).strang_hydro(time, dt, a_old, a_new);
        } 
#else
           get_level(lev).strang_hydro(time, dt, a_old, a_new);
#endif
    }
    BL_PROFILE_VAR_STOP(just_the_hydro);

    //
    // We must reflux before doing the next gravity solve
    //
    if (do_reflux)
    {
        for (int lev = level; lev < finest_level_to_advance; lev++)
        {
            get_level(lev).reflux();
        }
    }

    // Always average down the new state from finer to coarser.
    for (int lev = finest_level_to_advance-1; lev >= level; lev--)
    {
        get_level(lev).average_down(  State_Type);
        get_level(lev).average_down(DiagEOS_Type);
    }

#ifdef GRAVITY

    //
    // Here we use the "old" phi from the current time step as a guess for this
    // solve
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        MultiFab::Copy(parent->getLevel(lev).get_new_data(PhiGrav_Type),
                       parent->getLevel(lev).get_old_data(PhiGrav_Type),
                       0, 0, 1, 0);
    }

    // Solve for new Gravity
    BL_PROFILE_VAR("solve_for_new_phi", solve_for_new_phi);
    int use_previous_phi_as_guess = 1;
    if (finest_level_to_advance > level)
    {
        // The particle may be as many as "iteration" ghost cells out
        int ngrow_for_solve = iteration + stencil_deposition_width;
        gravity->multilevel_solve_for_new_phi(level, finest_level_to_advance, 
                                              ngrow_for_solve,
                                              use_previous_phi_as_guess);
    }
    else
    {
        int fill_interior = 0;
        gravity->solve_for_new_phi(level,get_new_data(PhiGrav_Type),
                               gravity->get_grad_phi_curr(level),
                               fill_interior, grav_n_grow);
    }
    BL_PROFILE_VAR_STOP(solve_for_new_phi);

    // Reflux
    if (do_reflux)
    {
        for (int lev = level; lev <= finest_level_to_advance; lev++)
        {
            gravity->add_to_fluxes(lev, iteration, ncycle);
        }
    }

    //
    // Now do corrector part of source term update
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        MultiFab& S_old = get_level(lev).get_old_data(State_Type);
        MultiFab& S_new = get_level(lev).get_new_data(State_Type);
        MultiFab& D_new = get_level(lev).get_new_data(DiagEOS_Type);
        MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
        reset_e_src.setVal(0.0);

        const auto& ba = get_level(lev).get_new_data(State_Type).boxArray();
        const auto& dm = get_level(lev).get_new_data(State_Type).DistributionMap();

        // These vectors are only used for the call to correct_gsrc so they 
        //    don't need any ghost cells
        MultiFab grav_vec_old(ba, dm, BL_SPACEDIM, 0);
        MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, 0);

        get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
        get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);

        const Real* dx = get_level(lev).Geom().CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            Real se  = 0;
            Real ske = 0;

            fort_correct_gsrc
                (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(grav_vec_old[mfi]),
                 BL_TO_FORTRAN(grav_vec_new[mfi]), BL_TO_FORTRAN(S_old[mfi]),
                 BL_TO_FORTRAN(S_new[mfi]), &a_old, &a_new, &dt);
        }

        // First reset internal energy before call to compute_temp
        get_level(lev).reset_internal_energy(S_new,D_new,reset_e_src);
        get_level(lev).compute_new_temp(S_new,D_new);
    }

    // Must average down again after doing the gravity correction;
    //      always average down from finer to coarser.
    // Here we average down both the new state and phi and gravity.
    for (int lev = finest_level_to_advance-1; lev >= level; lev--)
        get_level(lev).average_down();

    if (Nyx::theActiveParticles().size() > 0)
    {
        // Advance the particle velocities by dt/2 to the new time. We use the
        // cell-centered gravity to correctly interpolate onto particle
        // locations.
        if (particle_move_type == "Gravitational")
        {
            const Real a_half = 0.5 * (a_old + a_new);

            if (particle_verbose && ParallelDescriptor::IOProcessor())
                std::cout << "moveKick ... updating velocity only\n";

            for (int lev = level; lev <= finest_level_to_advance; lev++)
            {
                const auto& ba = get_level(lev).get_new_data(State_Type).boxArray();
                const auto& dm = get_level(lev).get_new_data(State_Type).DistributionMap();
                MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, grav_n_grow);
                get_level(lev).gravity->get_new_grav_vector(lev, grav_vec_new, cur_time);

                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                    Nyx::theActiveParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);

                // Virtual particles will be recreated, so we need not kick them.

                // Ghost particles need to be kicked except during the final iteration.
                if (iteration != ncycle)
                    for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
                        Nyx::theGhostParticles()[i]->moveKick(grav_vec_new, lev, dt, a_new, a_half);
            }
        }
    }
#endif

    //
    // Synchronize Energies
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        MultiFab& S_new = get_level(lev).get_new_data(State_Type);
        MultiFab& D_new = get_level(lev).get_new_data(DiagEOS_Type);
        MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
	reset_e_src.setVal(0.0);

	get_level(lev).reset_internal_energy(S_new,D_new,reset_e_src);
	
    }

    BL_PROFILE_REGION_STOP("R::Nyx::advance_hydro_plus_particles");

    // Redistribution happens in post_timestep
    return dt;
}

Real
Nyx::advance_hydro (Real time,
                    Real dt,
                    int  iteration,
                    int  ncycle)
{
    syncAllWorkerThreads();
    BL_PROFILE("Nyx::advance_hydro()");
    // sanity checks
    if (!do_hydro)
        amrex::Abort("In `advance_hydro` but `do_hydro` not true");
#ifdef GRAVITY
    if (!do_grav)
        amrex::Abort("In `advance_hydro` with GRAVITY defined but `do_grav` is false");
#endif
#ifdef FORCING
    if (!do_forcing)
        amrex::Abort("In `advance_hydro` with FORCING defined but `do_forcing` is false");
#endif

    Real prev_time;// = state[State_Type].prevTime();
    Real cur_time;//  = state[State_Type].curTime();
    Real a_old;//     = get_comoving_a(prev_time);
    Real a_new;//     = get_comoving_a(cur_time);

if(wtid()==0){
    for (int k = 0; k < NUM_STATE_TYPE; k++)
    {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }
    prev_time = state[State_Type].prevTime();
    cur_time  = state[State_Type].curTime();
    a_old     = get_comoving_a(prev_time);
    a_new     = get_comoving_a(cur_time);

#ifdef GRAVITY
    const int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
    {
        gravity->zero_phi_flux_reg(level + 1);
    }

    gravity->swap_time_levels(level);

#endif

#ifdef FORCING
    if (do_forcing) 
    {
        forcing->evolve(dt);
    }
#endif
}
syncAllWorkerThreads();
    // Call the hydro advance itself
    BL_PROFILE_VAR("just_the_hydro", just_the_hydro);
#ifdef SDC
    if (sdc_split > 0) 
    { 
       sdc_hydro(time, dt, a_old, a_new);
    } else { 
       strang_hydro(time, dt, a_old, a_new);
    } 
#else
       strang_hydro(time, dt, a_old, a_new);
#endif
    BL_PROFILE_VAR_STOP(just_the_hydro);

syncAllWorkerThreads();
if(wtid()==0){
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
    reset_e_src.setVal(0.0);

#ifdef GRAVITY
    if (verbose && ParallelDescriptor::IOProcessor()) 
        std::cout << "\n... new-time level solve at level " << level << '\n';

    MultiFab& S_old = get_old_data(State_Type);

    //
    // Solve for new phi
    // Here we use the "old" phi from the current time step as a guess for this solve
    //
    MultiFab::Copy(get_new_data(PhiGrav_Type),get_old_data(PhiGrav_Type),0,0,1,0);
    int fill_interior = 0;
    gravity->solve_for_new_phi(level,get_new_data(PhiGrav_Type),
                               gravity->get_grad_phi_curr(level),fill_interior);
    if (do_reflux)  
        gravity->add_to_fluxes(level,iteration,ncycle);

    // These are only used in correct_gsrc so don't need any ghost cells.
    MultiFab grav_vec_old(grids,dmap,BL_SPACEDIM,0);
    gravity->get_old_grav_vector(level,grav_vec_old,time);

    MultiFab grav_vec_new(grids,dmap,BL_SPACEDIM,0);
    gravity->get_new_grav_vector(level,grav_vec_new,cur_time);

    // Now do corrector part of source term update
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        fort_correct_gsrc
            (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(grav_vec_old[mfi]),
             BL_TO_FORTRAN(grav_vec_new[mfi]), BL_TO_FORTRAN(S_old[mfi]),
             BL_TO_FORTRAN(S_new[mfi]), &a_old, &a_new, &dt);
    }

#endif /*GRAVITY*/

    // First reset internal energy before call to compute_temp
    reset_internal_energy(S_new,D_new,reset_e_src);
    compute_new_temp(S_new,D_new);
}//wtid==0
syncAllWorkerThreads();

    return dt;
}
#endif

