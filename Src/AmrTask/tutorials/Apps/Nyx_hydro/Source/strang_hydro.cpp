
#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>

#ifdef GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;
using namespace perilla;

using std::string;

#ifndef NO_HYDRO
void
Nyx::strang_hydro (Real time,
                   Real dt,
                   Real a_old,
                   Real a_new)
{
    BL_PROFILE("Nyx::strang_hydro()");

    BL_ASSERT(NUM_GROW == 4);

    const Real prev_time    = state[State_Type].prevTime();
    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();

    MultiFab&  S_old        = get_old_data(State_Type);
    MultiFab&  S_new        = get_new_data(State_Type);

    MultiFab&  D_old        = get_old_data(DiagEOS_Type);
    MultiFab&  D_new        = get_new_data(DiagEOS_Type);

#ifndef NDEBUG
    if (std::abs(time-prev_time) > (1.e-10*cur_time) )
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "strang_hydro:  prev_time = " << prev_time << std::endl;
            std::cout << "strang_hydro:       time = " <<      time << std::endl;
        }
        amrex::Abort("time should equal prev_time in strang_hydro!");
    }

    if (S_old.contains_nan(Density, S_old.nComp(), 0))
    {
        for (int i = 0; i < S_old.nComp(); i++)
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "strang_hydro: testing component " << i << " for NaNs" << std::endl;
            if (S_old.contains_nan(Density+i,1,0))
                amrex::Abort("S_old has NaNs in this component");
        }
    }
#endif

    // It's possible for interpolation to create very small negative values for
    // species so we make sure here that all species are non-negative after this
    // point
    enforce_nonnegative_species(S_old);

syncAllWorkerThreads();
if(wtid()==0){
    //MultiFab ext_src_old(grids, dmap, NUM_STATE, NUM_GROW);
    //ext_src_old.setVal(0.);

    //if (add_ext_src)
       //get_old_source(prev_time, dt, ext_src_old);

    // Define the gravity vector 
    //MultiFab grav_vector(grids, dmap, BL_SPACEDIM, NUM_GROW);
   // grav_vector.setVal(0.);

#ifdef GRAVITY
    gravity->get_old_grav_vector(level, grav_vector, time);
    grav_vector.FillBoundary(geom.periodicity());
#endif

    // Create FAB for extended grid values (including boundaries) and fill.
    //MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
    FillPatch(*this, *S_old_tmp, NUM_GROW, time, State_Type, 0, NUM_STATE);

    //MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);
    FillPatch(*this, *D_old_tmp, NUM_GROW, time, DiagEOS_Type, 0, D_old.nComp());

    //MultiFab hydro_src(grids, dmap, NUM_STATE, 0);
    //hydro_src.setVal(0.);

    //MultiFab divu_cc(grids, dmap, 1, 0);
    //divu_cc.setVal(0.);

#ifdef HEATCOOL
    strang_first_step(time,dt,*S_old_tmp,*D_old_tmp);
#endif

    bool   init_flux_register = true;
    bool add_to_flux_register = true;
    compute_hydro_sources(time,dt,a_old,a_new,*S_old_tmp,*D_old_tmp,
                          *ext_src_old,*hydro_src,*grav_vector,*divu_cc,
                          init_flux_register, add_to_flux_register);

    update_state_with_sources(*S_old_tmp,S_new,
                              *ext_src_old,*hydro_src,*grav_vector,*divu_cc,
                              dt,a_old,a_new);

    // We copy old Temp and Ne to new Temp and Ne so that they can be used
    //    as guesses when we next need them.
    MultiFab::Copy(D_new,D_old,0,0,D_old.nComp(),0);

    //grav_vector->clear();

    if (add_ext_src)
    {
        get_old_source(prev_time, dt, *ext_src_old);
        // Must compute new temperature in case it is needed in the source term evaluation
        compute_new_temp(S_new,D_new);

        // Compute source at new time (no ghost cells needed)
        MultiFab ext_src_new(grids, dmap, NUM_STATE, 0);
        ext_src_new.setVal(0);

        get_new_source(prev_time, cur_time, dt, ext_src_new);

        time_center_source_terms(S_new, *ext_src_old, ext_src_new, dt);

        compute_new_temp(S_new,D_new);
    } // end if (add_ext_src)


#ifndef NDEBUG
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
        amrex::Abort("S_new has NaNs before the second strang call");
#endif

#ifdef HEATCOOL
    // This returns updated (rho e), (rho E), and Temperature
    strang_second_step(cur_time,dt,S_new,D_new);
#endif

#ifndef NDEBUG
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
        amrex::Abort("S_new has NaNs after the second strang call");
#endif

}//wtid==0
syncAllWorkerThreads();

}
#endif
