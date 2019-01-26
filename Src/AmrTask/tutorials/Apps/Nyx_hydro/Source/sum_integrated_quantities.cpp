#include <iostream>
#include <iomanip>

#include <Nyx.H>
#include <Nyx_F.H>

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef FORCING
#include "Forcing.H"
#endif


using namespace amrex;

#ifndef NO_HYDRO
void
Nyx::sum_integrated_quantities ()
{
    // This is a static in the Nyx class -- we need to update this every
    // time step
    compute_average_density();

    if ((parent->NumDataLogs() < 2) && (verbose <= 0))
        return;

    if (do_hydro == 0)
        return;

    if (verbose <= 0)
        return;

    int num_global_data_logs = 2;

    const Real box_vol = geom.ProbSize();

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real mass = 0;
    Real xmom = 0;
    Real ymom = 0;
    Real zmom = 0;
    Real rho_E_lev    = 0, rho_E    = 0;
    Real rho_e_lev    = 0, rho_e    = 0;
    Real Temp_lev     = 0, Temp     = 0;
    Real rms_mach_lev = 0, rms_mach = 0;
    Real magvort_lev = 0, magvort = 0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Nyx& nyx_lev = get_level(lev);

        mass += nyx_lev.vol_weight_sum("density", time, true);
        xmom += nyx_lev.vol_weight_sum("xmom", time, true);
        ymom += nyx_lev.vol_weight_sum("ymom", time, true);
        zmom += nyx_lev.vol_weight_sum("zmom", time, true);
        rho_E += nyx_lev.vol_weight_sum("rho_E", time, true);
        rho_e += nyx_lev.vol_weight_sum("rho_e", time, true);
        Temp  += nyx_lev.vol_weight_sum("Temp", time, true);

        magvort      += nyx_lev.vol_weight_squared_sum("magvort", time);
        rms_mach     += nyx_lev.vol_weight_squared_sum("MachNumber", time);

        if (parent->NumDataLogs() >= num_global_data_logs + lev + 1)
        {
  	    if (ParallelDescriptor::IOProcessor())
                std::cout << "Computing statistics for level " << lev << " at time " << time << '\n';

            rho_E_lev = nyx_lev.vol_weight_sum("rho_E", time, false);
            rho_e_lev = nyx_lev.vol_weight_sum("rho_e", time, false);
            Temp_lev  = nyx_lev.vol_weight_sum ("Temp", time, false);

            rms_mach_lev = nyx_lev.vol_weight_squared_sum_level("MachNumber", time);
            magvort_lev  = nyx_lev.vol_weight_squared_sum_level("magvort", time);
        }

        rms_mach_lev = std::sqrt(rms_mach_lev);
         magvort_lev = std::sqrt( magvort_lev);

        if (ParallelDescriptor::IOProcessor() && parent->NumDataLogs() >= num_global_data_logs + lev + 1) 
        {

            std::cout << "Writing level " << lev << " statistics to data log #" << num_global_data_logs + lev << '\n';
            std::ostream& data_log_lev = parent->DataLog(num_global_data_logs + lev);

            if (time == 0)
            {
                data_log_lev << std::setw(14) << "      time    ";
                data_log_lev << std::setw(14) << "        rho_E ";
                data_log_lev << std::setw(14) << "        rho_e ";
                data_log_lev << std::setw(14) << "         Temp ";
                data_log_lev << std::setw(14) << "     rms_mach ";
                data_log_lev << std::setw(14) << "      magvort "
                             << std::endl;
            }

            // Write the quantities at this time      
                                                                                                  
            data_log_lev << std::setw(14) << time;
            data_log_lev << std::setw(14) << std::setprecision(6) << rho_E_lev;
            data_log_lev << std::setw(14) << std::setprecision(6) << rho_e_lev;
            data_log_lev << std::setw(14) << std::setprecision(6) << Temp_lev;
            data_log_lev << std::setw(14) << std::setprecision(6) << rms_mach_lev;
            data_log_lev << std::setw(14) << std::setprecision(6) << magvort_lev
                         << std::endl;
        }
    }

    // normalization by box volume
    mass /= box_vol;
    xmom /= box_vol;
    ymom /= box_vol;
    zmom /= box_vol;
    rho_E    /= box_vol;
    rho_e    /= box_vol;
    Temp     /= box_vol;
    rms_mach /= box_vol;
    rms_mach = std::sqrt(rms_mach);
    magvort /= box_vol;
    magvort = std::sqrt(magvort);

    if (ParallelDescriptor::IOProcessor())
    {
      if (verbose > 0) 
      {
          std::cout << '\n';
          std::cout << "BOX VOLUME= " << box_vol << '\n';
          std::cout << "TIME= " << time << " MASS        = " << mass << '\n';
          std::cout << "TIME= " << time << " XMOM        = " << xmom << '\n';
          std::cout << "TIME= " << time << " YMOM        = " << ymom << '\n';
          std::cout << "TIME= " << time << " ZMOM        = " << zmom << '\n';
          std::cout << "TIME= " << time << " RHO*E       = " << rho_E << '\n';
      }

      if (parent->NumDataLogs() >= 2)
      {
	std::ostream& data_log1 = parent->DataLog(1);

        if (time == 0)
        {
            data_log1 << std::setw(14) << "      time    ";
#ifdef FORCING
            data_log1 << std::setw(14) << "    rms_force ";
#endif
            data_log1 << std::setw(14) << "         xmom ";
            data_log1 << std::setw(14) << "         ymom ";
            data_log1 << std::setw(14) << "         zmom ";
            data_log1 << std::setw(14) << "        rho_E ";
            data_log1 << std::setw(14) << "        rho_e ";
            data_log1 << std::setw(14) << "         Temp ";
            data_log1 << std::setw(14) << "     rms_mach ";
            data_log1 << std::setw(14) << "      magvort " << std::endl;
        }

        // Write the quantities at this time                                                                                            
        data_log1 << std::setw(14) << time;
#ifdef FORCING
        data_log1 << std::setw(14) << std::setprecision(6) << forcing->rms();
#endif
        data_log1 << std::setw(14) << std::setprecision(6) << xmom;
        data_log1 << std::setw(14) << std::setprecision(6) << ymom;
        data_log1 << std::setw(14) << std::setprecision(6) << zmom;
        data_log1 << std::setw(14) << std::setprecision(6) << rho_E;
        data_log1 << std::setw(14) << std::setprecision(6) << rho_e;
        data_log1 << std::setw(14) << std::setprecision(6) << Temp;
        data_log1 << std::setw(14) << std::setprecision(6) << rms_mach;
        data_log1 << std::setw(14) << std::setprecision(6) << magvort
                  << std::endl;
      }
      std::cout << std::endl;
    }
}

void
Nyx::compute_average_density ()
{
    int             finest_level = parent->finestLevel();
    Real            time         = state[State_Type].curTime();
    const Geometry& crse_geom    = parent->Geom(0);

    // This is a static in the Nyx class
    average_gas_density      = 0;
    average_dm_density       = 0;
    average_neutr_density    = 0;
    average_total_density    = 0;

    // Add up the baryonic density
    if (do_hydro == 1)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            Nyx& nyx_lev = get_level(lev);
            average_gas_density += nyx_lev.vol_weight_sum("density", time, true);
        }
    }
 
#ifdef GRAVITY
    // Define the dark matter density on all levels.
    if (Nyx::theDMPC())
    {
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensity(particle_mf);

        // Note that we don't need to call the average_down routine because the 
        //   vol_weight_sum routine zeroes out the data under the finer grids
 
        // Add up the dark matter density; here masked is true because we don't want 
        //     to double count the coarse values under fine grids.
        for (int lev = 0; lev <= finest_level; lev++)
        {
            Nyx& nyx_lev = get_level(lev);
            average_dm_density += nyx_lev.vol_weight_sum(*particle_mf[lev],true);
        }
    }
#ifdef NEUTRINO_PARTICLES
    if (Nyx::theNPC())
    {
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensity(particle_mf);

        // Note that we don't need to call the average_down routine because the 
        //   vol_weight_sum routine zeroes out the data under the finer grids
 
        // Add up the neutrino density
        for (int lev = 0; lev <= finest_level; lev++)
        {
            Nyx& nyx_lev = get_level(lev);
            average_neutr_density += nyx_lev.vol_weight_sum(*particle_mf[lev],true);
        }
    }
#endif
#endif

    // Divide by physical volume of domain.
    if (do_hydro == 1)
        average_gas_density  /= crse_geom.ProbSize();
    average_dm_density       /= crse_geom.ProbSize();
    average_neutr_density    /= crse_geom.ProbSize();

    // Define the total density = gas density + dark matter density
    if (do_hydro == 1)
    {
        average_total_density = average_gas_density + average_dm_density + average_neutr_density;
    }
    else
    {
        average_total_density = average_dm_density + average_neutr_density;
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        if (do_hydro == 1)
        {
            std::cout << "Average     gas density " << average_gas_density << '\n';
        }
        std::cout << "Average       dm density " << average_dm_density << '\n';
#ifdef NEUTRINO_PARTICLES
        std::cout << "Average neutrino density " << average_neutr_density << '\n';
#endif
        std::cout << "Average    total density " << average_total_density << '\n';
    }
}

void
Nyx::compute_average_temperature (Real& average_temperature)
{
    int             finest_level = parent->finestLevel();
    Real            time         = state[State_Type].curTime();
    const Geometry& crse_geom    = parent->Geom(0);

    // Add up the temperature -- this is just volume-weighted, not mass-weighted
    average_temperature = 0;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        Nyx& nyx_lev = get_level(lev);
        MultiFab& S_new = nyx_lev.get_new_data(State_Type);
        MultiFab& D_new = nyx_lev.get_new_data(DiagEOS_Type);

	MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
	reset_e_src.setVal(0.0);

        nyx_lev.reset_internal_energy(S_new,D_new,reset_e_src);
        nyx_lev.compute_new_temp     (S_new,D_new);

        average_temperature += nyx_lev.vol_weight_sum("Temp",time,true);
    }
 
    // Divide by physical volume of domain.
    average_temperature = average_temperature / crse_geom.ProbSize();

    if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
        std::cout << "Average temperature " << average_temperature << '\n';
    }
}

void
Nyx::compute_average_species (int          nspec,
                              int          naux,
                              Vector<Real>& average_species)
{
    if (use_const_species == 1)
    {
       average_species[0] = h_species;
       average_species[1] = he_species;
       if (verbose > 0 && ParallelDescriptor::IOProcessor())
       {
           std::cout << "Average species " << 0 << ": " << average_species[0] << '\n';
           std::cout << "Average species " << 1 << ": " << average_species[1] << '\n';
       }
    }
    else
    {
        Real                 time = state[State_Type].curTime();
        const Geometry& crse_geom = parent->Geom(0);
        //
        // Get the species names from the network model.
        //
        for (int i = 0; i < nspec+naux; i++)
        {
            int len = 20;
            Vector<int> int_names(len);
            //
            // This call return the actual length of each string in "len"
            //
            char* name = new char[len+1];

            if (i < nspec)
            {
               fort_get_spec_names(int_names.dataPtr(), &i, &len);
            }
            else
            {
               fort_get_aux_names(int_names.dataPtr(), &i, &len);
            }
    
            for (int j = 0; j < len; j++)
                name[j] = int_names[j];
            name[len] = '\0';

            std::string spec_string;

            if (i < nspec)
            {
               spec_string = "X(";
               spec_string += name;
               spec_string += ')';
            }
            else
            {
               spec_string = name;
            }

            // Add up the species -- this is just volume-weighted, not mass-weighted
            for (int lev = 0; lev <= parent->finestLevel(); lev++)
            {
                Nyx& nyx_lev = get_level(lev);
                average_species[i] += nyx_lev.vol_weight_sum(spec_string, time, true);
            }
            delete [] name;
 
           // Divide by physical volume of domain.
           average_species[i] = average_species[i] / crse_geom.ProbSize();

           if (verbose > 0 && ParallelDescriptor::IOProcessor())
               std::cout << "Average species " << i << ": " << average_species[i] << '\n';
        }
    } // end if not use_const_species
}
#endif
