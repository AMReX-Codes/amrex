#include <Nyx.H>
#ifdef GRAVITY
#include <Gravity.H>
#endif

using namespace amrex;

std::unique_ptr<MultiFab>
Nyx::particle_derive (const std::string& name, Real time, int ngrow)
{
    if (Nyx::theDMPC() && name == "particle_count")
    {
	std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
        MultiFab temp_dat(grids, dmap, 1, 0);
        temp_dat.setVal(0);
        Nyx::theDMPC()->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
#ifdef AGN
    else if (Nyx::theAPC() && name == "agn_particle_count")
    {
	std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
        MultiFab temp_dat(grids, dmap, 1, 0);
        temp_dat.setVal(0);
        Nyx::theAPC()->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
#endif
#ifdef NEUTRINO_PARTICLES
    else if (Nyx::theNPC() && name == "neutrino_particle_count")
    {
	std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids, dmap, 1, 0));
        MultiFab temp_dat(grids, dmap, 1, 0);
        temp_dat.setVal(0);
        Nyx::theNPC()->Increment(temp_dat, level);
        MultiFab::Copy(*derive_dat, temp_dat, 0, 0, 1, 0);
        return derive_dat;
    }
#endif
    else if (Nyx::theDMPC() && name == "total_particle_count")
    {
        //
        // We want the total particle count at this level or higher.
        //
	std::unique_ptr<MultiFab> derive_dat = particle_derive("particle_count", time, ngrow);
        IntVect trr(D_DECL(1, 1, 1));

        // @todo: level vs. lev
        for (int lev = level + 1; lev <= parent->finestLevel(); lev++)
        {
            auto ba = parent->boxArray(lev);
            const auto& dm = parent->DistributionMap(lev);
            MultiFab temp_dat(ba, dm, 1, 0);

            trr *= parent->refRatio(lev - 1);

            ba.coarsen(trr);
            MultiFab ctemp_dat(ba, dm, 1, 0);

            temp_dat.setVal(0);
            ctemp_dat.setVal(0);

            Nyx::theDMPC()->Increment(temp_dat, lev);

            for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
            {
                const FArrayBox& ffab = temp_dat[mfi];
                FArrayBox& cfab = ctemp_dat[mfi];
                const Box& fbx = ffab.box();

                BL_ASSERT(cfab.box() == amrex::coarsen(fbx, trr));

                for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
                {
                    const Real val = ffab(p);
                    if (val > 0)
                        cfab(amrex::coarsen(p, trr)) += val;
                }
            }

            temp_dat.clear();

            MultiFab dat(grids, dmap, 1, 0);
            dat.setVal(0);
            dat.copy(ctemp_dat);

            MultiFab::Add(*derive_dat, dat, 0, 0, 1, 0);
        }

        return derive_dat;
    }
#ifdef GRAVITY
    else if (Nyx::theDMPC() && name == "particle_mass_density")
    {
	std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        MultiFab::Copy(*derive_dat, *particle_mf[level], 0, 0, 1, 0);

        return derive_dat;
    }
#ifdef AGN
    else if (Nyx::theAPC() && name == "agn_mass_density")
    {
	std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theAPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        MultiFab::Copy(*derive_dat, *particle_mf[level], 0, 0, 1, 0);

        return derive_dat;
    }
#endif
#ifdef NEUTRINO_PARTICLES
    else if (Nyx::theNPC() && name == "neutrino_mass_density")
    {
	std::unique_ptr<MultiFab> derive_dat(new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theNPC()->AssignDensity(particle_mf);

        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        MultiFab::Copy(*derive_dat, *particle_mf[level], 0, 0, 1, 0);

        return derive_dat;
    }
#endif
#endif
    else if (name == "total_density")
    {
#ifdef GRAVITY
      if (Nyx::theDMPC())
      {
        std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids,dmap,1,0));

        // We need to do the multilevel `assign_density` even though we're only
        // asking for one level's worth because otherwise we don't get the
        // coarse-fine distribution of particles correct.
        Vector<std::unique_ptr<MultiFab> > particle_mf;
        Nyx::theDMPC()->AssignDensity(particle_mf);
       
        for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
        {
            amrex::average_down(*particle_mf[lev+1], *particle_mf[lev], 
                                 parent->Geom(lev+1), parent->Geom(lev), 0, 1, 
                                 parent->refRatio(lev));
        }

        MultiFab::Copy(*derive_dat, *particle_mf[level], 0, 0, 1, 0);

#ifndef NO_HYDRO
	std::unique_ptr<MultiFab> gas_density = derive("density",time,0);
        MultiFab::Add(*derive_dat,*gas_density, 0, 0, 1, 0);
#endif
        return derive_dat; 
      }
#ifndef NO_HYDRO
      else 
      {
        return derive("density",time,0);
      }
#endif

#else
        return derive("density",time,0);
#endif
    }
    else
    {
        return AmrLevel::derive(name, time, ngrow);
    }
}
