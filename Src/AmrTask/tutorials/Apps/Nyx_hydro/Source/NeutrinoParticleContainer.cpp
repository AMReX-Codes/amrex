#include "NyxParticleContainer.H"

using namespace amrex;

#ifdef NEUTRINO_PARTICLES
void
NeutrinoParticleContainer::AssignDensity (Vector<std::unique_ptr<MultiFab> >& mf, int lev_min, int ncomp, int finest_level) const 
{
    AssignRelativisticDensity (mf,lev_min,ncomp,finest_level);
}

void
NeutrinoParticleContainer::AssignRelativisticDensity (Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
                                                      int               lev_min,
                                                      int               ncomp,
                                                      int               finest_level) const
{
    BL_PROFILE("NeutrinoParticleContainer::AssignRelativisticDensity()");
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);
    BL_ASSERT(m_csq > 0.0);

    if (finest_level == -1)
    {
        finest_level = m_gdb->finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level))
    {
        finest_level--;
    }
    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level, 
    // lev_index is the corresponding index for mf. 
    //

    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1-lev_min);
    for (int lev = lev_min; lev <= finest_level; lev++)
    { 
        const int lev_index = lev - lev_min;
        mf_to_be_filled[lev].reset(new MultiFab(m_gdb->boxArray(lev), 
						m_gdb->DistributionMap(lev), 
						ncomp, 1));
	mf_to_be_filled[lev]->setVal(0.0);
    }

    // Test whether the grid structure of the boxArray is the same
    //       as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++) {
        if (!OnSameGrids(lev, *mf_to_be_filled[lev-lev_min])) {
	    all_grids_the_same = false;
	    break;
	}
    }

    Vector<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1-lev_min);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_part[lev_index].reset(new MultiFab(m_gdb->ParticleBoxArray(lev),
						  m_gdb->ParticleDistributionMap(lev),
						  ncomp, 1));
	    mf_part[lev_index].setVal(0.0);
        }
    }

    auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;

    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        AssignRelativisticDensitySingleLevel(*mf[0],0,ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0]->copy(*mf[0],0,0,ncomp);
	}
        return;
    }
    
    const bool sub_cycle = m_gdb->subCycle();
    //
    // This'll hold all the info I need for parallel.
    // // What I'll use: m_lev, m_grid, m_cell & m_data[0..ncomp-1].
    //
    // This is the "data" needed by other MPI procs.
    //
    PMap data;

    const Real stime = ParallelDescriptor::second();
    //
    // Minimum M required.
    //
    const int M = D_TERM(2,+2,+4);

    Vector<int>     cgrid(M);
    Vector<int>    cwhich(M),  fwhich(M);
    Vector<Real>    fracs(M),  cfracs(M);
    Vector<IntVect> cells(M),  ccells(M), cfshifts(M);

    ParticleType pb;
    //
    // I'm going to allocate these badboys here & pass'm into routines that use'm.
    // This should greatly cut down on memory allocation/deallocation.
    //
    Vector<IntVect>                    pshifts(27);
    std::vector< std::pair<int,Box> > isects;
    Vector<int>                        fgrid(M);
    Vector<Real>                       ffracs(M);
    Vector<IntVect>                    fcells;
    //
    // "fvalid" contains all the valid region of the MultiFab at this level, together
    // with any ghost cells lying outside the domain, that can be periodically shifted into the
    // valid region.  "compfvalid" is the complement of the "fvalid", while "compfvalid_grown" is 
    // "compfvalid" grown by one.  Using these we can figure out whether or not a cell is in the
    // valid region of our MultiFab as well as whether or not we're at a Fine->Crse boundary.
    //
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const Geometry& gm        = m_gdb->Geom(lev);
        const Geometry& gm_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1) : gm;
        const Geometry& gm_coarse = (lev > 0) ? m_gdb->Geom(lev-1) : gm;
        const Box&      dm        = gm.Domain();
        const Real*     dx        = gm.CellSize();
        const Real*     plo       = gm.ProbLo();
        const Real*     dx_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1).CellSize() : dx;
        const Real*     dx_coarse = (lev > 0) ? m_gdb->Geom(lev-1).CellSize() : dx;
        const int       lev_index = lev - lev_min;
        const BoxArray& grids     = mf[lev_index]->boxArray();
        const int       dgrow     = (lev == 0) ? 1 : m_gdb->MaxRefRatio(lev-1);

        BoxArray compfvalid, compfvalid_grown, fvalid = mf[lev_index]->boxArray();
        //
        // Do we have Fine->Crse overlap on a periodic boundary?
        // We want to add all ghost cells that can be shifted into valid region.
        //
        BoxList valid;

        for (int i = 0; i < grids.size(); i++)
        {
            if (gm.isAnyPeriodic())
            {
                const Box& dest = amrex::grow(grids[i],dgrow);

                if ( ! dm.contains(dest))
                {
                    for (int j = 0; j < grids.size(); j++)
                    {
                        BL_ASSERT(dm.contains(grids[j]));

                        gm.periodicShift(dest, grids[j], pshifts);

                        for (int k = 0; k < pshifts.size(); k++)
                        {
                            const Box& sbx = grids[j] + pshifts[k];
                            const Box& dbx = dest & sbx;

                            BL_ASSERT(dbx.ok());

                            valid.push_back(dbx);
                        }
                    }
                }
            }
        }
        if (valid.isNotEmpty())
        {
            //
            // We've got some Fine->Crse periodic overlap.
            // Don't forget to add the valid boxes too.
            //
            for (int i = 0; i < grids.size(); i++) {
                valid.push_back(grids[i]);
	    }
            fvalid = BoxArray(valid);
            fvalid.removeOverlap();
        }
        //
        // If we're at a lev < finestLevel, this is the coarsened fine BoxArray.
        // We use this for figuring out Crse->Fine issues.
        //
        BoxArray ccba;
        if (lev > 0)
        {
            ccba = m_gdb->boxArray(lev);
            ccba.coarsen(m_gdb->refRatio(lev-1));
        }
        BoxArray cfba;
        if (lev < finest_level)
        {
            cfba = m_gdb->boxArray(lev+1);
            cfba.coarsen(m_gdb->refRatio(lev));

            BL_ASSERT(mf[lev_index]->boxArray().contains(cfba));
        }
        //
        // This is cfba with any shifted ghost cells.
        //
        BoxArray cfvalid = cfba;

        if (lev < finest_level)
        {
            BoxList cvalid;

            const BoxArray& cgrids = mf[lev_index]->boxArray();

            for (int i = 0; i < cfba.size(); i++)
            {
                if (gm.isAnyPeriodic())
                {
                    const Box& dest = amrex::grow(cfba[i],mf[lev_index]->nGrow());

                    if ( ! dm.contains(dest))
                    {
                        for (int j = 0; j < cgrids.size(); j++)
                        {
                            BL_ASSERT(dm.contains(cgrids[j]));

                            gm.periodicShift(dest, cgrids[j], pshifts);

                            for (int k = 0; k < pshifts.size(); k++)
                            {
                                const Box& sbx = cfba[i] - pshifts[k];

                                cvalid.push_back(sbx);
                            }
                        }
                    }
                }
            }
            if (cvalid.isNotEmpty())
            {
                //
                // We've got some Fine->Crse periodic overlap.
                // Don't forget to add the valid boxes too.
                //
                for (int i = 0; i < cfba.size(); i++) {
                    cvalid.push_back(cfba[i]);
		}
                cfvalid = BoxArray(cvalid);
                cfvalid.removeOverlap();
            }
        }
        //
        // The "+1" is so we enclose the valid region together with any
        //  ghost cells that can be periodically shifted into valid.
        //
        BoxList bl_compfvalid(dm.ixType());
        bl_compfvalid.complementIn(amrex::grow(dm,dgrow+1), fvalid);
        compfvalid = BoxArray(std::move(bl_compfvalid));

        compfvalid_grown = compfvalid;
        compfvalid_grown.grow(1);
        compfvalid_grown.removeOverlap();
            
        if (gm.isAnyPeriodic() && ! gm.isAllPeriodic())
        {
            amrex::Error("AssignDensity: problem must be periodic in no or all directions");
        }
        //
        // If we're at a lev > 0, this is the coarsened BoxArray.
        // We use this for figuring out Fine->Crse issues.
        //
        BoxArray cba;
        if (lev > 0)
        {
            cba = m_gdb->boxArray(lev);
            cba.coarsen(m_gdb->refRatio(lev-1));
        }
        //
        // Do the grids at this level cover the full domain? If they do
        // there can be no Fine->Crse interactions at this level.
        //
        const bool GridsCoverDomain = fvalid.contains(m_gdb->Geom(lev).Domain());
        
        for (typename PMap::const_iterator pmap_it = m_particles[lev].begin(),
                 pmapEnd = m_particles[lev].end();
             pmap_it != pmapEnd;
             ++pmap_it)
        {
            const PBox& pbx = pmap_it->second;
            FArrayBox&  fab = (*mf[lev_index])[pmap_it->first];

            for (typename PBox::const_iterator it = pbx.begin(), End = pbx.end();
                 it != End;
                 ++it)
            {
                const ParticleType& p = *it;

                if (p.m_id <= 0) {
		  continue;
		}
                //
                // Get "fracs" and "cells" for the particle "p" at this level.
                //
                const int M = ParticleBase::CIC_Cells_Fracs(p, plo, dx, fracs, cells);
                //
                // If this is not fully periodic then we have to be careful that no
                // particle's support leaves the domain. We test this by checking the low
                // and high corners respectively.
                //
                if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                    if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                        amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		    }
		}
                //
                // This section differs based on whether we subcycle.
                // Without subcycling we use the "stretchy" support for particles.
                // With subcycling a particles support is strictly defined 
                // by its resident level.
                //
                if (sub_cycle)
                {
                    bool isFiner    = false;
                    bool isBoundary = false;
                    //
                    // First sum the mass in the valid region
                    //
                    for (int i = 0; i < M; i++)
                    {
                        if (cfvalid.contains(cells[i]))
                        {
                            //
                            // Some part of the particle's mass lies in a 
                            // finer region; we'll deal with it shortly.
                            //
                            isFiner    = true;
                            isBoundary = true;
                            continue;
                        }
                        if ( ! fvalid.contains(cells[i]))
                        {
                            //
                            // We're out of the valid region.
                            //
                            isBoundary = true;
                            continue;
                        }
                        //
                        // Sum up mass in first component.
                        //
                        if (m_relativistic)
                        {
                            Real vsq = 0.0;
                            for (int n = 1; n < ncomp; n++) {
                               vsq += p.m_data[n] * p.m_data[n];
			    }
                            Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                            fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                        }
                        else 
                        {
                            fab(cells[i],0) += p.m_data[0] * fracs[i];
                        }
                        //
                        // Sum up momenta in next components.
                        //

                        // If the domain is not periodic and we want to let particles
                        //    live near the boundary but "throw away" the contribution that 
                        //    does not fall into the domain ...
                        if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
			     ! gm.Domain().contains(cells[i]))
			{
			  continue;
			}

                        for (int n = 1; n < ncomp; n++) {
                            fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
			}
                    }
                    //
                    // Deal with mass that doesn't belong at this level.
                    // Here we assume proper nesting so that only one special case can
                    // be true for a given particle.
                    //
                    if (isBoundary)
                    {
                        if (isFiner)
                        {
                            BL_ASSERT(lev < finest_level);
                            //
                            // We're at a coarse->fine interface
                            //
                            // get fine cells/fracs
                            //
                            const int MF = ParticleBase::CIC_Cells_Fracs(p, plo, dx_fine ,dx, ffracs, fcells);

                            for (int j = 0; j < MF; j++)
                            {
                                //
                                // Make sure this fine cell is valid. Check for periodicity.
                                //
                                const Box bx(fcells[j],fcells[j]);
                                gm_fine.periodicShift(bx, gm_fine.Domain(), pshifts);
                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(pshifts.size() == 1);
                                    fcells[j] = fcells[j] - pshifts[0];
                                }
                                mf[lev_index + 1]->boxArray().intersections(Box(fcells[j],fcells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first; 
                                const int who  = mf[lev_index+1]->DistributionMap()[grid];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        (*mf[lev_index+1])[grid](fcells[j],0) += p.m_data[0] * ffracs[j] * gamma;
                                    }
                                    else 
                                    {
                                        (*mf[lev_index+1])[grid](fcells[j],0) += p.m_data[0] * ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index+1])[grid](fcells[j],n) += p.m_data[n] * p.m_data[0] * ffracs[j];
				    }
                                }
                                else
                                {
                                    pb.m_lev  = lev+1;
                                    pb.m_grid = grid;
                                    pb.m_cell = fcells[j];
                                    //
                                    // Sum up mass in first component.
                                    //
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        pb.m_data[0] = p.m_data[0] *  ffracs[j] * gamma;
                                    }
                                    else 
                                    {
                                        pb.m_data[0] = p.m_data[0] *  ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        pb.m_data[n] = p.m_data[n] * p.m_data[0] * ffracs[j];
				    }
                                    data[who].push_back(pb);
                                }
                            }
                        }
                        else if (lev_index > 0)
                        {
                            //
                            // We must be at a fine->coarse interface.
                            //
                            const int MC = ParticleBase::CIC_Cells_Fracs(p, plo, dx_coarse, dx, cfracs, ccells);
                            for (int j = 0; j < MC; j++)
                            {
                                //
                                // Make sure this coarse cell isn't in this level's valid region.
                                // This may not matter.
                                //
                                if (cba.contains(ccells[j]))
                                    continue;
                                //
                                // Check for periodicity.
                                //
                                const Box bx(ccells[j],ccells[j]);
                                gm_coarse.periodicShift(bx, gm_coarse.Domain(), pshifts);

                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(pshifts.size() == 1);
                                    ccells[j] = ccells[j] - pshifts[0]; 
                                }
                                //
                                // Find its resident grid.
                                //
                                mf[lev_index - 1]->boxArray().intersections(Box(ccells[j],ccells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first;
                                const int who  = mf[lev_index-1]->DistributionMap()[grid];
                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        (*mf[lev_index-1])[grid](ccells[j],0) += p.m_data[0] * cfracs[j] * gamma;
                                    }
                                    else 
                                    {
                                        (*mf[lev_index-1])[grid](ccells[j],0) += p.m_data[0] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[grid](ccells[j],n) += p.m_data[n] * p.m_data[0] * cfracs[j];
				    }
                                }
                                else
                                {
                                    pb.m_lev  = lev-1;
                                    pb.m_grid = grid;
                                    pb.m_cell = ccells[j];
                                    //
                                    // Sum up mass in first component.
                                    //
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        pb.m_data[0] = p.m_data[0] * cfracs[j] * gamma;
                                    }
                                    else 
                                    {
                                        pb.m_data[0] = p.m_data[0] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        pb.m_data[n] = p.m_data[n] * p.m_data[0] * cfracs[j];
				    }

                                    data[who].push_back(pb);
                                }
                            }
                        }
                        else
                        {
                            // The mass is below levels we care about. Ignore it.
                        }
                    }
                }
                else 
                {
                    bool AnyCrseToFine = false;
                    if (lev < finest_level) {
                        AnyCrseToFine = ParticleBase::CrseToFine(cfba,cells,cfshifts,gm,cwhich,pshifts);
		    }
                    //
                    // lev_index > 0 means that we don't do F->C for lower levels
                    // This may mean that the mass fraction is off.
                    //
                    bool AnyFineToCrse = false;
                    if (lev_index > 0 && !GridsCoverDomain)
                        AnyFineToCrse = ParticleBase::FineToCrse(p,lev,m_gdb,cells,fvalid,compfvalid_grown,ccells,cfracs,fwhich,cgrid,pshifts,isects);

                    BL_ASSERT(!(AnyCrseToFine && AnyFineToCrse));

                    if ( ! AnyCrseToFine && ! AnyFineToCrse)
                    {
                        //
                        // By far the most common case.  Just do it!
                        //
                        for (int i = 0; i < M; i++)
                        {

                            // If the domain is not periodic and we want to let particles
                            //    live near the boundary but "throw away" the contribution that 
                            //    does not fall into the domain ...
                            if (! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i]))
			    {
			      continue;
			    }
                            //
                            // Sum up mass in first component.
                            //
                            if (m_relativistic)
                            {
                                Real vsq = 0.0;
                                for (int n = 1; n < ncomp; n++) {
                                   vsq += p.m_data[n] * p.m_data[n];
				}
                                Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                            }
                            else 
                            {
                                fab(cells[i],0) += p.m_data[0] * fracs[i];
                            }
                            //
                            // Sum up momenta in next components.
                            //
                            for (int n = 1; n < ncomp; n++) {
                                fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
			    }
                        }
                    }
                    else if (AnyFineToCrse)
                    {
                        Real sum_crse = 0, sum_fine = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // We're at a Fine->Crse boundary.
                                //
                                BL_ASSERT(cgrid[i] >= 0);
                                BL_ASSERT(cgrid[i] < mf[lev_index-1]->size());
                                //
                                // Here we need to update the crse region.  The coarse
                                // region is always going to be updated if we have a
                                // particle in a cell bordering a Fine->Crse boundary.
                                //
                                const int who = mf[lev_index-1]->DistributionMap()[cgrid[i]];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    if ( ! (*mf[lev_index-1])[cgrid[i]].box().contains(ccells[i])) {
				      continue;
				    }

                                    // If the domain is not periodic and we want to let particles
                                    //    live near the boundary but "throw away" the contribution that 
                                    //    does not fall into the domain ...
                                    if (! gm_coarse.isAllPeriodic() && allow_particles_near_boundary &&
				        ! gm_coarse.Domain().contains(ccells[i]))
				    {
				      continue;
				    }

                                    //
                                    // Sum up mass in first component.
                                    //
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],0) += p.m_data[0] * cfracs[i] * gamma;
                                    }
                                    else 
                                    {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],0) += p.m_data[0] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],n) += p.m_data[n] * p.m_data[0] * cfracs[i];
				    }
                                }
                                else
                                {
                                    pb.m_lev  = lev-1;
                                    pb.m_grid = cgrid[i];
                                    pb.m_cell = ccells[i];
                                    //
                                    // Sum up mass in first component.
                                    //
                                    if (m_relativistic)
                                    {
                                        Real vsq = 0.0;
                                        for (int n = 1; n < ncomp; n++) {
                                           vsq += p.m_data[n] * p.m_data[n];
					}
                                        Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                        pb.m_data[0] = p.m_data[0] * cfracs[i] * gamma;
                                    }
                                    else 
                                    {
                                        pb.m_data[0] = p.m_data[0] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        pb.m_data[n] = p.m_data[n] * p.m_data[0] * cfracs[i];
				    }
                                    data[who].push_back(pb);
                                }

                                sum_crse += cfracs[i];
                            }
                        }
                        //
                        // We've updated the Crse cells.  Now we have to update the fine
                        // cells in such a way that the total amount of mass we move
                        // around is precisely p.m_data[0]. In other words, the fractions
                        // we use at crse and fine have to sum to zero.  In the fine
                        // case, we have to account for the case where one or more of the
                        // cell indices is not in the valid region of the box containing 
                        // the particle.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {
                            //
                            // Reusing "fwhich" to indicate fine cells that need massaging.
                            //
                            fwhich[i] = true;

                            if ( ! compfvalid_grown.contains(cells[i]))
                            {
                                //
                                // Go ahead and add the full correct amount to these cells.
                                // They can't touch a Fine->Crse boundary.
                                //
                                sum_fine += fracs[i];
                                //
                                // Sum up mass in first component.
                                //
                                if (m_relativistic)
                                {
                                    Real vsq = 0.0;
                                    for (int n = 1; n < ncomp; n++) {
                                       vsq += p.m_data[n] * p.m_data[n];
				    }
                                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                                }
                                else 
                                {
                                    fab(cells[i],0) += p.m_data[0] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
				}
                                fwhich[i] = false;
                            }
                            else if (compfvalid.contains(cells[i]))
                            {
                                fwhich[i] = false;
                            }
                        }

                        const Real sum_so_far = sum_crse + sum_fine; 

                        BL_ASSERT(sum_so_far > 0);
                        BL_ASSERT(sum_so_far < 1);

                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {       
                            if (fwhich[i])
                                //
                                // Got to weight cells in this direction differently.
                                //
                                sum_fine += fracs[i];
                        }

                        const Real mult = (1 - sum_so_far) / sum_fine;
                        //
                        // Now add the weighted amount to the fine cells touching the c-f interface.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // Sum up mass in first component.
                                //
                                if (m_relativistic)
                                {
                                    Real vsq = 0.0;
                                    for (int n = 1; n < ncomp; n++) {
                                       vsq += p.m_data[n] * p.m_data[n];
				    }
                                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * mult * gamma;
                                }
                                else 
                                {
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * mult;
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i] * mult;
				}

                                sum_fine += fracs[i] * mult;
                            }
                        }

                        BL_ASSERT(std::abs(1-(sum_fine+sum_so_far)) < 1.e-9);
                    }
                    else if (AnyCrseToFine)
                    {
                        Real sum = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (!cwhich[i])
                            {
                                // If the domain is not periodic and we want to let particles
                                //    live near the boundary but "throw away" the contribution that 
                                //    does not fall into the domain ...
                                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
				     ! gm.Domain().contains(ccells[i]))
				{
				  continue;
				}
                                //
                                // Sum up mass in first component.
                                //
                                if (m_relativistic)
                                {
                                    Real vsq = 0.0;
                                    for (int n = 1; n < ncomp; n++) {
                                       vsq += p.m_data[n] * p.m_data[n];
				    }
                                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                    fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                                }
                                else 
                                {
                                    fab(cells[i],0) += p.m_data[0] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
				}

                                sum += fracs[i];
                            }
                            else
                            {
                                //
                                // We're at a Crse->Fine boundary.
                                //
                                ParticleBase::FineCellsToUpdateFromCrse(p,lev,m_gdb,cells[i],cfshifts[i],fgrid,ffracs,fcells,isects);

                                for (int j = 0; j < fcells.size(); j++)
                                {
                                    const int who = mf[lev_index+1]->DistributionMap()[fgrid[j]];

                                    if (who == ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
                                        if (m_relativistic)
                                        {
                                            Real vsq = 0.0;
                                            for (int n = 1; n < ncomp; n++) {
                                               vsq += p.m_data[n] * p.m_data[n];
					    }
                                            Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],0) += p.m_data[0] * fracs[i] * ffracs[j] * gamma;
                                        }
                                        else 
                                        {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],0) += p.m_data[0] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],n) += p.m_data[n] * p.m_data[0] * fracs[i] * ffracs[j];
					}
                                    }
                                    else
                                    {
                                        pb.m_lev  = lev+1;
                                        pb.m_grid = fgrid[j];
                                        pb.m_cell = fcells[j];
                                        //
                                        // Sum up mass in first component.
                                        //
                                        if (m_relativistic)
                                        {
                                            Real vsq = 0.0;
                                            for (int n = 1; n < ncomp; n++) {
                                               vsq += p.m_data[n] * p.m_data[n];
					    }
                                            Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
                                            pb.m_data[0] = p.m_data[0] * fracs[i] * ffracs[j] * gamma;
                                        }
                                        else 
                                        {
                                            pb.m_data[0] = p.m_data[0] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            pb.m_data[n] = p.m_data[n] * p.m_data[0] * fracs[i] * ffracs[j];
					}

                                        data[who].push_back(pb);
                                    }

                                    sum += fracs[i] * ffracs[j];
                                }
                            }
                        }

                        BL_ASSERT(std::abs(1-sum) < 1.e-9);
                    }
                }
            }
        }
    }
    //
    // Send any needed data to other MPI processes.
    // This "may" touch ghost cells so we want to do it before
    // the SumBoundary() stuff.
    //
    AssignDensityDoit(0, mf,data,ncomp,lev_min);

    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const int       lev_index = lev - lev_min;
        const Geometry& gm        = m_gdb->Geom(lev);
        const Real*     dx        = gm.CellSize();
        const Real      vol       = D_TERM(dx[0], *dx[1], *dx[2]);

        mf[lev_index]->SumBoundary(gm.periodicity());
        //
        // If ncomp > 1, first divide the momenta (component n) 
        // by the mass (component 0) in order to get velocities.
        // Be careful not to divide by zero.
        //
        for (int n = 1; n < ncomp; n++)
        {
            for (MFIter mfi(*mf[lev_index]); mfi.isValid(); ++mfi)
            {
                (*mf[lev_index])[mfi].protected_divide((*mf[lev_index])[mfi],0,n,1);
            }
        }
        //
        // Only multiply the first component by (1/vol) because this converts mass
        // to density. If there are additional components (like velocity), we don't
        // want to divide those by volume.
        //
        mf[lev_index]->mult(1/vol,0,1);
    }

    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level,  
    // lev_index is the corresponding index for mf. 
    //
    // I believe that we don't need any information in ghost cells so we don't copy those.
    //
    if ( ! all_grids_the_same)
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_to_be_filled[lev_index]->copy(*mf_part[lev_index],0,0,1);
        }
    
    if (m_verbose > 1)
    {
        Real etime = ParallelDescriptor::second() - stime;

        ParallelDescriptor::ReduceRealMax(etime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "NeutrinoParticleContainer::AssignRelativisticDensity(multi-level) time: " << etime << '\n';
        }
    }
}

//
// This is the single-level version for cell-centered density
//
void
NeutrinoParticleContainer::AssignRelativisticDensitySingleLevel (MultiFab& mf_to_be_filled,
                                                         int       lev,
                                                         int       ncomp,
                                                         int       particle_lvl_offset) const
{
    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), 
				  m_gdb->ParticleDistributionMap(lev),
				  ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread its effect to an
    //    adjacent grid by first putting the value into ghost cells of its own grid.  The mf->sumBoundary call then
    //    adds the value from one grid's ghost cell to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //
    // TODO: implement tiling with OpenMP in this grid loop.
    Vector<int>         pgrd(ngrids);
    Vector<const PBox*> pbxs(ngrids);

    int j = 0;
    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it, ++j)
    {
        pgrd[j] =   pmap_it->first;
        pbxs[j] = &(pmap_it->second);
    }

    for (int j = 0; j < ngrids; j++)
    {
        const PBox& pbx = *pbxs[j];
        FArrayBox&  fab = (*mf_pointer)[pgrd[j]];

        Vector<Real>    fracs;
        Vector<IntVect> cells;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(fracs,cells) shared(plo,dx,dx_particle,gm,fab,ncomp,pbx)
#endif
        for (typename PBox::const_iterator it = pbx.begin(); it < pbx.end(); ++it)
        {
            const ParticleType& p = *it;

            if (p.m_id <= 0) {
	      continue;
	    }

            const int M = ParticleBase::CIC_Cells_Fracs(p, plo, dx, dx_particle, fracs, cells);
            //
            // If this is not fully periodic then we have to be careful that the
            // particle's support leaves the domain unless we specifically want to ignore
            // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
            // We test this by checking the low and high corners respectively.
            //
            if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                    amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		}
	    }

            for (int i = 0; i < M; i++)
            {
                if ( ! fab.box().contains(cells[i])) {
		  continue;
		}

                // If the domain is not periodic and we want to let particles
                //    live near the boundary but "throw away" the contribution that 
                //    does not fall into the domain ...
                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i])) {
		  continue;
		}
                //
                // Sum up mass in first component.
                //
                if (m_relativistic)
                {
                    Real vsq = 0.0;
                    for (int n = 1; n < ncomp; n++) {
                       vsq += p.m_data[n] * p.m_data[n];
		    }
                    Real gamma = 1.0 / sqrt(1.0 - vsq / m_csq);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    fab(cells[i],0) += p.m_data[0] * fracs[i] * gamma;
                }
                else 
                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    fab(cells[i],0) += p.m_data[0] * fracs[i];
                }
                // 
                // Sum up momenta in next components.
                //
                for (int n = 1; n < ncomp; n++)
#ifdef _OPENMP
#pragma omp atomic
#endif
                   fab(cells[i],n) += p.m_data[n] * p.m_data[0] * fracs[i];
            }
        }
    }

    mf_pointer->SumBoundary(gm.periodicity());
    //
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    //
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }
    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "NeutrinoParticleContainer::AssignRelativisticDensitySingleLevel time: " << stoptime << '\n';
        }
    }
}
#endif

