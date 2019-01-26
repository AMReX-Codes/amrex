#include <iomanip>
#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;

Real Nyx::initial_z             = -1.0;
Real Nyx::final_a               = -1.0;
Real Nyx::final_z               = -1.0;
Real Nyx::relative_max_change_a =  0.01;
Real Nyx::absolute_max_change_a = -1.0;
Real Nyx::dt_binpow             = -1.0;

void
Nyx::read_comoving_params ()
{
    ParmParse pp("nyx");

    pp.query("initial_z", initial_z);
    pp.query("final_a",   final_a);
    pp.query("final_z",   final_z);

    if (final_z >= 0)
    {
        if (final_a > 0)
        {
            std::cerr << "ERROR::dont specify both final_a and final_z\n";
            amrex::Error();
        }
        else
        {
            final_a = 1 / (1 + final_z);
        }
    }

    pp.query("relative_max_change_a", relative_max_change_a);
    pp.query("absolute_max_change_a", absolute_max_change_a);
    pp.query("dt_binpow",             dt_binpow);


    // for shrinking box tests, initial_z < 0 is ok
    if (initial_z < 0)
    {
        std::cerr << "ERROR::Need to specify non-negative initial redshift \n";
        amrex::Error();
    }
}

void
Nyx::comoving_est_time_step (Real& cur_time, Real& estdt)
{
    Real change_allowed = relative_max_change_a;
    Real fixed_da = absolute_max_change_a;
    Real dt             = estdt;
    Real new_dummy_a;
    int  dt_modified;

    if ( std::abs(cur_time - new_a_time) <= 1.e-12 * cur_time)
    {

        // Initial guess -- note that we send in "new_a" because we haven't yet swapped
        // "old_a" and "new_a" -- we can't do that until after we compute dt and then
        // integrate a forward.
        fort_estdt_comoving_a
	  (&new_a, &new_dummy_a, &dt, &change_allowed, &fixed_da, &final_a, &dt_modified);

    if(dt_binpow >= 0)
      {
	if(estdt>=dt)
	  estdt=dt;
	else if(estdt>.5*dt)
	  {
	    estdt=.5*dt;
	    //	    std::cout << "Lavel = 1" <<std::endl;
	  }
	else if(estdt>.25*dt)
	  {
	    estdt=.25*dt;
	    //	    std::cout << "Lavel = 2" <<std::endl;
	  }
	else if(estdt>.125*dt)
	  {
	    estdt=.125*dt;
	    //	    std::cout << "Lavel = 3" <<std::endl;
	  }
	else if(estdt>.0625*dt)
	  {
	    estdt=.0625*dt;
	    //	    std::cout << "Lavel = 4" <<std::endl;
	  }
	else
	  {
	    //dta*(2**(-1*np.ceil( np.log2(dta/dth))))
	    estdt = dt*(pow(2,(-std::ceil( std::log2(dt/estdt)))));
	    //	    std::cout << "Lavel > 4" <<std::endl;
	  }
	fort_integrate_comoving_a(&new_a,&new_dummy_a,&estdt);
      }
    else
      {
	estdt=std::min(estdt,dt);
      }
          

        if (verbose && (dt_modified == 1) && ParallelDescriptor::IOProcessor())
        {
            std::cout << "...estdt after call to comoving   : "
                      << dt
                      << "\n...change in a is "
                      << (new_dummy_a - new_a) / new_a * 100.0
                      << " percent\n";
        }
    } 
    // This only happens with the second call at time = 0 where we want to re-advance
    //      from old_a_time time, not start at new_a_time
    else if ( std::abs(cur_time - old_a_time) <= 1.e-12 * cur_time)
    {

        // Initial guess -- note that we send in "new_a" because we haven't yet swapped
        // "old_a" and "new_a" -- we can't do that until after we compute dt and then
        // integrate a forward.
        fort_estdt_comoving_a
	  (&old_a, &new_dummy_a, &dt, &change_allowed, &fixed_da, &final_a, &dt_modified);
    if(dt_binpow >= 0)
      {
	if(estdt>=dt)
	  estdt=dt;
	else if(estdt>.5*dt)
	  {
	    estdt=.5*dt;
	    //	    std::cout << "Lavel = 1" <<std::endl;
	  }
	else if(estdt>.25*dt)
	  {
	    estdt=.25*dt;
	    //	    std::cout << "Lavel = 2" <<std::endl;
	  }
	else if(estdt>.125*dt)
	  {
	    estdt=.125*dt;
	    //	    std::cout << "Lavel = 3" <<std::endl;
	  }
	else if(estdt>.0625*dt)
	  {
	    estdt=.0625*dt;
	    //	    std::cout << "Lavel = 4" <<std::endl;
	  }
	else
	  {
	    //dta*(2**(-1*np.ceil( np.log2(dta/dth))))
	    estdt = dt*(pow(2,(-std::ceil( std::log2(dt/estdt)))));
	    //	    std::cout << "Lavel > 4" <<std::endl;
	  }
	fort_integrate_comoving_a(&old_a,&new_dummy_a,&estdt);
      }
    else
      {
	estdt=std::min(estdt,dt);
      }

        if (verbose && (dt_modified == 1) && ParallelDescriptor::IOProcessor())
        {
            std::cout << "...advancing from old_a_time rather than new_a_time! " << std::endl;
            std::cout << "...estdt after call to comoving   : "
                      << dt
                      << "\n...change in a is "
                      << (new_dummy_a - old_a) / old_a * 100.0
                      << " percent\n";
        }

    } 
    else if ( std::abs(cur_time - old_a_time) <= 1.e-12 * cur_time)
    {
       std::cout << "comoving_est_time_step: DONT KNOW WHAT TIME IT IS " << cur_time << std::endl;
       exit(0);
    } 

    return;
}

Real
Nyx::get_comoving_a (Real time)
{
    const Real eps         = 0.0001 * (new_a_time - old_a_time);

    // Test on whether time == old_a_time == new_a_time -- for example after restart
    //   before a has been integrated for the new time step.
    if ( ( std::abs(time - old_a_time) <= 1.e-12*old_a_time ) &&
         ( std::abs(time - new_a_time) <= 1.e-12*new_a_time ) &&
         (  old_a == new_a ) )
    {
        return old_a;
    }
    else if (time > old_a_time - eps && time < old_a_time + eps)
    {
        return old_a;
    }
    else if (time > new_a_time - eps && time < new_a_time + eps)
    {
        return new_a;
    }
    else if (time > old_a_time && time < new_a_time)
    {
        Real frac = (time - old_a_time) / (new_a_time - old_a_time);
        Real    a = frac*new_a + (1.0-frac)*old_a;
        return a;
    } 
    else
    {
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Invalid:get comoving_a at " << time << std::endl; 
            std::cout << "Old / new a_time  " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new   a     " << old_a      << " " << new_a    << std::endl;
        }
        amrex::Error("get_comoving_a: invalid time");
        return 0;
    }
}

void
Nyx::integrate_comoving_a (Real time,Real dt)
{
    if (level > 0) 
        return;

    bool first;

    if ( std::abs(time-new_a_time) <= (1.e-10 * time) )
    {
       first = true;
    } else {
       first = false;
    } 

    if (first) 
    {

        // Update a
        old_a      = new_a;
        fort_integrate_comoving_a(&old_a, &new_a, &dt);

        // Update the times
        old_a_time = new_a_time;
        new_a_time = old_a_time + dt;
 
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Integrating a from time " << time << " by dt = " << dt << '\n';
            std::cout << "Old / new A time      " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new A           " << old_a      << " " << new_a      << std::endl;
            std::cout << "Old / new z           " << 1./old_a-1.<< " " << 1./new_a-1. << std::endl;
        }
    }
    else if (std::abs(time-old_a_time) <= 1.e-10 * time) 
    {
        // Leave old_a and old_a_time alone -- we have already swapped them
        fort_integrate_comoving_a(&old_a, &new_a, &dt);
            (&old_a, &new_a, &dt);

        // Update the new time only
        new_a_time = old_a_time + dt;

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Re-integrating a from time " << time << " by dt = " << dt << '\n';
            std::cout << "Old / new A time         " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new A              " << old_a      << " " << new_a      << std::endl;
            std::cout << "Old / new z              " << 1./old_a-1.<< " " << 1./new_a-1. << std::endl;
        }
    }
    else 
    {
            std::cout << "Time passed to integrate_comoving_a " << time << std::endl;
            std::cout << "Old / new A time                    " << old_a_time << " " << new_a_time << std::endl;
            std::cout << "Old / new A                         " << old_a      << " " << new_a      << std::endl;
            std::cerr << "ERROR::dont know what to do in integrate_comoving_a" << std::endl;
            amrex::Error();
    }
}

void
Nyx::comoving_a_post_restart (const std::string& restart_file)
{
    if (level > 0)
        return;

    if (ParallelDescriptor::IOProcessor())
    {
        std::string FileName = restart_file + "/comoving_a";
        std::ifstream File;
        File.open(FileName.c_str(),std::ios::in);
        if (!File.good())
            amrex::FileOpenFailed(FileName);
        File >> old_a;
    }
    ParallelDescriptor::Bcast(&old_a, 1, ParallelDescriptor::IOProcessorNumber());

    new_a = old_a;

#ifdef NO_HYDRO
    old_a_time = state[PhiGrav_Type].curTime();
#else
    old_a_time = state[State_Type].curTime();
#endif
    new_a_time = old_a_time;
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "...setting old_a_time to " << old_a_time << std::endl;
    }

#ifdef HEATCOOL
     // Initialize "this_z" in the atomic_rates_module
     if (heat_cool_type == 1 || heat_cool_type == 3 || heat_cool_type == 5 || heat_cool_type == 7) {
         Real old_z = 1.0/old_a - 1.0;
         fort_interp_to_this_z(&old_z);
     }
#endif
}

void
Nyx::plot_z_est_time_step (Real& dt_0, bool& dt_changed)
{
    Real dt = dt_0;
    Real a_old, z_old, a_new, z_new;

    // This is where we are now
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    a_old = get_comoving_a(cur_time);
    z_old = (1. / a_old) - 1.;

    // *****************************************************
    // First test whether we are within dt of a plot_z value
    // *****************************************************

    // This is where we would be if we use the current dt_0
    Real new_time = cur_time + dt_0;
    integrate_comoving_a (cur_time,dt_0);
    a_new = get_comoving_a(new_time);
    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the plot_z_values array
    Real z_value;
    bool found_one = false;
    for (int i = 0; i < plot_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - plot_z_values[i]) * (z_old - plot_z_values[i]) < 0 && !found_one)
        {
            z_value   = plot_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that dt_0 is too big and makes us pass z_value,
    // we must figure out what value of dt < dt_0 makes us exactly reach z_value
    if (found_one)
    {
        fort_integrate_comoving_a_to_z(&old_a, &z_value, &dt);
                

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << dt << std::endl;
            std::cout << " ... in order to write a plotfile at z = " << z_value << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = dt;
        dt_changed = true;
    } 
    else 
    { 

    // *****************************************************
    // If not within one dt, now test whether we are within 2*dt of a plot_z value
    // *****************************************************

    // This is where we would be if we advance by twice the current dt_0
    Real two_dt = 2.0*dt_0;

    fort_integrate_comoving_a(&a_old, &a_new, &two_dt);

    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the plot_z_values array
    Real z_value;
    bool found_one = false;
    for (int i = 0; i < plot_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - plot_z_values[i]) * (z_old - plot_z_values[i]) < 0 && !found_one)
        {
            z_value   = plot_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that 2*dt_0 will make us pass z_value, we set the current dt
    // as half the interval to reach that z_value
    if (found_one)
    {
        Real two_dt = 2.0*dt_0;
        fort_integrate_comoving_a_to_z(&old_a, &z_value, &two_dt);
                

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << 0.5 * two_dt << std::endl;
            std::cout << " ... in order to write a plotfile at z = " << z_value 
                      << " two steps from now " << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = 0.5 * two_dt;
        dt_changed = true;
    }

    }
}

void
Nyx::analysis_z_est_time_step (Real& dt_0, bool& dt_changed)
{
    Real dt = dt_0;
    Real a_old, z_old, a_new, z_new;

    // This is where we are now
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    a_old = get_comoving_a(cur_time);
    z_old = (1. / a_old) - 1.;

    // *****************************************************
    // First test whether we are within dt of a analysis_z value
    // *****************************************************

    // This is where we would be if we use the current dt_0
    Real new_time = cur_time + dt_0;
    integrate_comoving_a (cur_time,dt_0);
    a_new = get_comoving_a(new_time);
    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the analysis_z_values array
    Real z_value;
    bool found_one = false;
    for (int i = 0; i < analysis_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - analysis_z_values[i]) * (z_old - analysis_z_values[i]) < 0 && !found_one)
        {
            z_value   = analysis_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that dt_0 is too big and makes us pass z_value,
    // we must figure out what value of dt < dt_0 makes us exactly reach z_value
    if (found_one)
    {
        fort_integrate_comoving_a_to_z(&old_a, &z_value, &dt);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << dt << std::endl;
            std::cout << " ... in order to do analysis at z = " << z_value << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = dt;
        dt_changed = true;
    } 
    else 
    { 

    // *****************************************************
    // If not within one dt, now test whether we are within 2*dt of a analysis_z value
    // *****************************************************

    // This is where we would be if we advance by twice the current dt_0
    Real two_dt = 2.0*dt_0;

    fort_integrate_comoving_a(&a_old, &a_new, &two_dt); 
    z_new = (1. / a_new) - 1.;

    // Find the relevant entry of the analysis_z_values array
    Real z_value;
    bool found_one = false;
    for (int i = 0; i < analysis_z_values.size(); i++)
    {
        // We have gone from before to after one of the specified values
        if ( (z_new - analysis_z_values[i]) * (z_old - analysis_z_values[i]) < 0 && !found_one)
        {
            z_value   = analysis_z_values[i];
            found_one = true;
        }
    }

    // Now that we know that 2*dt_0 will make us pass z_value, we set the current dt
    // as half the interval to reach that z_value
    if (found_one)
    {
        Real two_dt = 2.0*dt_0;
        fort_integrate_comoving_a_to_z(&old_a, &z_value, &two_dt);

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            std::cout << " " << std::endl;
            std::cout << " ... modifying time step from " << dt_0 << " to " << 0.5 * two_dt << std::endl;
            std::cout << " ... in order to do analysis at z = " << z_value 
                      << " two steps from now " << std::endl;
            std::cout << " " << std::endl;
        }

        // We want to pass this value back out
        dt_0 = 0.5 * two_dt;
        dt_changed = true;
    }

    }
}
