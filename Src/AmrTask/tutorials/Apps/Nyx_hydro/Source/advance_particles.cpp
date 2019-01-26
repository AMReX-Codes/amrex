#ifdef GRAVITY
 
#include <AMReX_BLProfiler.H>
#include "Nyx.H"
#include "Nyx_F.H"
#include "Gravity.H"
#include <AMReX_Particles_F.H>
#include <Gravity_F.H>

using namespace amrex;
 
using std::string;

Real
Nyx::advance_particles_only (Real time,
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
     BL_PROFILE("Nyx::advance_particles_only()");

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
 
    // *** grav_n_grow *** is used
    //   *) to determine how many ghost cells we need to fill in the MultiFab from
    //      which the particle interpolates its acceleration
    //   *) to set how many cells the Where call in moveKickDrift tests = (grav.nGrow()-2).
 
    int grav_n_grow = ghost_width + (1-iteration) +
                      stencil_interpolation_width ;

    // Sanity checks
    if (do_hydro)
        amrex::Abort("In `advance_particles_only` but `do_hydro` is true");

    if (!do_grav)
        amrex::Abort("In `advance_particles_only` but `do_grav` not true");
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
#ifndef NO_HYDRO
            MultiFab& S_old = get_level(lev).get_old_data(State_Type);
            MultiFab& S_new = get_level(lev).get_new_data(State_Type);
            MultiFab::Copy(S_new, S_old, 0, 0, S_old.nComp(), 0);
#endif
        }
    }

    const Real prev_time = state[PhiGrav_Type].prevTime();
    const Real cur_time  = state[PhiGrav_Type].curTime();

    const Real a_old     = get_comoving_a(prev_time);
    const Real a_new     = get_comoving_a(cur_time);

    //
    // We now do a multilevel solve for old Gravity. This goes to the 
    // finest level regardless of subcycling behavior. Consequentially,
    // If we are subcycling we skip this step on the first iteration of
    // finer levels.
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
        // Solve for phi
        // If a single-level calculation we can still use the previous phi as a guess.
        // TODO: Check this.
        int use_previous_phi_as_guess = 1;
        gravity->multilevel_solve_for_old_phi(level, finest_level,
                                              use_previous_phi_as_guess);
    }
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
                const auto& ba = get_level(lev).get_new_data(PhiGrav_Type).boxArray();
                const auto& dm = get_level(lev).get_new_data(PhiGrav_Type).DistributionMap();
                MultiFab grav_vec_old(ba, dm, BL_SPACEDIM, grav_n_grow);
                get_level(lev).gravity->get_old_grav_vector(lev, grav_vec_old, time);
                
                for (int i = 0; i < Nyx::theActiveParticles().size(); i++)
                    Nyx::theActiveParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_old, a_half);

                // Only need the coarsest virtual particles here.
                if (lev == level && level < finest_level)
                    for (int i = 0; i < Nyx::theVirtualParticles().size(); i++)
                        Nyx::theVirtualParticles()[i]->moveKickDrift(grav_vec_old, level, dt, a_old, a_half);

                // Miiiight need all Ghosts
                for (int i = 0; i < Nyx::theGhostParticles().size(); i++)
                    Nyx::theGhostParticles()[i]->moveKickDrift(grav_vec_old, lev, dt, a_new, a_half);
            }
        }
    }

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
    int use_previous_phi_as_guess = 1;
    if (finest_level_to_advance > level)
    {
        gravity->multilevel_solve_for_new_phi(level, finest_level_to_advance, 
                                              use_previous_phi_as_guess);
    }
    else
    {
        int fill_interior = 0;
        gravity->solve_for_new_phi(level,get_new_data(PhiGrav_Type),
                               gravity->get_grad_phi_curr(level),
                               fill_interior, grav_n_grow);
    }

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
                const auto& ba = get_level(lev).get_new_data(PhiGrav_Type).boxArray();
                const auto& dm = get_level(lev).get_new_data(PhiGrav_Type).DistributionMap();
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

    // Redistribution happens in post_timestep
    return dt;
}
#endif
