/* Copyright 2019 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Maxence Thevenet, Michael Rowan
 * Remi Lehe, Weiqun Zhang, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <WarpX.H>
#include <WarpXAlgorithmSelection.H>
#include <AMReX_BLProfiler.H>

using namespace amrex;

void
WarpX::LoadBalance ()
{
    WARPX_PROFILE_REGION("LoadBalance");
    WARPX_PROFILE("WarpX::LoadBalance()");

    AMREX_ALWAYS_ASSERT(costs[0] != nullptr);

    if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic)
    {
        WarpX::ComputeCostsHeuristic(costs);
    }

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev)
    {
#ifdef AMREX_USE_MPI
        // Parallel reduce the costs_heurisitc
        amrex::Vector<Real>::iterator it = (*costs[lev]).begin();
        amrex::Real* itAddr = &(*it);
        ParallelAllReduce::Sum(itAddr,
                               costs[lev]->size(),
                               ParallelContext::CommunicatorSub());
#endif
        const amrex::Real nboxes = costs[lev]->size();
        const amrex::Real nprocs = ParallelContext::NProcsSub();
        const int nmax = static_cast<int>(std::ceil(nboxes/nprocs*load_balance_knapsack_factor));

        const DistributionMapping newdm = (load_balance_with_sfc)
            ? DistributionMapping::makeSFC(*costs[lev], boxArray(lev), false)
            : DistributionMapping::makeKnapSack(*costs[lev], nmax);

        RemakeLevel(lev, t_new[lev], boxArray(lev), newdm);
    }
    mypc->Redistribute();
}


void
WarpX::RemakeLevel (int lev, Real /*time*/, const BoxArray& ba, const DistributionMapping& dm)
{
    if (ba == boxArray(lev))
    {
        if (ParallelDescriptor::NProcs() == 1) return;

        // Fine patch
        for (int idim=0; idim < 3; ++idim)
        {
            {
                const IntVect& ng = Bfield_fp[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_fp[lev][idim]->boxArray(),
                                                                  dm, Bfield_fp[lev][idim]->nComp(), ng));
                pmf->Redistribute(*Bfield_fp[lev][idim], 0, 0, Bfield_fp[lev][idim]->nComp(), ng);
                Bfield_fp[lev][idim] = std::move(pmf);
            }
            {
                const IntVect& ng = Efield_fp[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_fp[lev][idim]->boxArray(),
                                                                  dm, Efield_fp[lev][idim]->nComp(), ng));
                pmf->Redistribute(*Efield_fp[lev][idim], 0, 0, Efield_fp[lev][idim]->nComp(), ng);
                Efield_fp[lev][idim] = std::move(pmf);
            }
            {
                const IntVect& ng = current_fp[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_fp[lev][idim]->boxArray(),
                                                                  dm, current_fp[lev][idim]->nComp(), ng));
                current_fp[lev][idim] = std::move(pmf);
            }
            if (current_store[lev][idim])
            {
                const IntVect& ng = current_store[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_store[lev][idim]->boxArray(),
                                                                  dm, current_store[lev][idim]->nComp(), ng));
                // no need to redistribute
                current_store[lev][idim] = std::move(pmf);
            }
        }

        if (F_fp[lev] != nullptr) {
            const IntVect& ng = F_fp[lev]->nGrowVect();
            auto pmf = std::unique_ptr<MultiFab>(new MultiFab(F_fp[lev]->boxArray(),
                                                              dm, F_fp[lev]->nComp(), ng));
            pmf->Redistribute(*F_fp[lev], 0, 0, F_fp[lev]->nComp(), ng);
            F_fp[lev] = std::move(pmf);
        }

        if (rho_fp[lev] != nullptr) {
            const int nc = rho_fp[lev]->nComp();
            const IntVect& ng = rho_fp[lev]->nGrowVect();
            auto pmf = std::unique_ptr<MultiFab>(new MultiFab(rho_fp[lev]->boxArray(),
                                                              dm, nc, ng));
            rho_fp[lev] = std::move(pmf);
        }

        // Aux patch
        if (lev == 0 && Bfield_aux[0][0]->ixType() == Bfield_fp[0][0]->ixType())
        {
            for (int idim = 0; idim < 3; ++idim) {
                Bfield_aux[lev][idim].reset(new MultiFab(*Bfield_fp[lev][idim], amrex::make_alias, 0, Bfield_aux[lev][idim]->nComp()));
                Efield_aux[lev][idim].reset(new MultiFab(*Efield_fp[lev][idim], amrex::make_alias, 0, Efield_aux[lev][idim]->nComp()));
            }
        } else {
            for (int idim=0; idim < 3; ++idim)
            {
                {
                    const IntVect& ng = Bfield_aux[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_aux[lev][idim]->boxArray(),
                                                                      dm, Bfield_aux[lev][idim]->nComp(), ng));
                    // pmf->Redistribute(*Bfield_aux[lev][idim], 0, 0, Bfield_aux[lev][idim]->nComp(), ng);
                    Bfield_aux[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = Efield_aux[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_aux[lev][idim]->boxArray(),
                                                                      dm, Efield_aux[lev][idim]->nComp(), ng));
                    // pmf->Redistribute(*Efield_aux[lev][idim], 0, 0, Efield_aux[lev][idim]->nComp(), ng);
                    Efield_aux[lev][idim] = std::move(pmf);
                }
            }
        }

        // Coarse patch
        if (lev > 0) {
            for (int idim=0; idim < 3; ++idim)
            {
                {
                    const IntVect& ng = Bfield_cp[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_cp[lev][idim]->boxArray(),
                                                                      dm, Bfield_cp[lev][idim]->nComp(), ng));
                    pmf->Redistribute(*Bfield_cp[lev][idim], 0, 0, Bfield_cp[lev][idim]->nComp(), ng);
                    Bfield_cp[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = Efield_cp[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_cp[lev][idim]->boxArray(),
                                                                      dm, Efield_cp[lev][idim]->nComp(), ng));
                    pmf->Redistribute(*Efield_cp[lev][idim], 0, 0, Efield_cp[lev][idim]->nComp(), ng);
                    Efield_cp[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = current_cp[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>( new MultiFab(current_cp[lev][idim]->boxArray(),
                                                                       dm, current_cp[lev][idim]->nComp(), ng));
                    current_cp[lev][idim] = std::move(pmf);
                }
            }

            if (F_cp[lev] != nullptr) {
                const IntVect& ng = F_cp[lev]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(F_cp[lev]->boxArray(),
                                                                  dm, F_cp[lev]->nComp(), ng));
                pmf->Redistribute(*F_cp[lev], 0, 0, F_cp[lev]->nComp(), ng);
                F_cp[lev] = std::move(pmf);
            }

            if (rho_cp[lev] != nullptr) {
                const int nc = rho_cp[lev]->nComp();
                const IntVect& ng = rho_cp[lev]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(rho_cp[lev]->boxArray(),
                                                                  dm, nc, ng));
                rho_cp[lev] = std::move(pmf);
            }
        }

        if (lev > 0 && (n_field_gather_buffer > 0 || n_current_deposition_buffer > 0)) {
            for (int idim=0; idim < 3; ++idim)
            {
                if (Bfield_cax[lev][idim])
                {
                    const IntVect& ng = Bfield_cax[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_cax[lev][idim]->boxArray(),
                                                                      dm, Bfield_cax[lev][idim]->nComp(), ng));
                    // pmf->ParallelCopy(*Bfield_cax[lev][idim], 0, 0, Bfield_cax[lev][idim]->nComp(), ng, ng);
                    Bfield_cax[lev][idim] = std::move(pmf);
                }
                if (Efield_cax[lev][idim])
                {
                    const IntVect& ng = Efield_cax[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_cax[lev][idim]->boxArray(),
                                                                      dm, Efield_cax[lev][idim]->nComp(), ng));
                    // pmf->ParallelCopy(*Efield_cax[lev][idim], 0, 0, Efield_cax[lev][idim]->nComp(), ng, ng);
                    Efield_cax[lev][idim] = std::move(pmf);
                }
                if (current_buf[lev][idim])
                {
                    const IntVect& ng = current_buf[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_buf[lev][idim]->boxArray(),
                                                                      dm, current_buf[lev][idim]->nComp(), ng));
                    // pmf->ParallelCopy(*current_buf[lev][idim], 0, 0, current_buf[lev][idim]->nComp(), ng, ng);
                    current_buf[lev][idim] = std::move(pmf);
                }
            }
            if (charge_buf[lev])
            {
                const IntVect& ng = charge_buf[lev]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(charge_buf[lev]->boxArray(),
                                                                  dm, charge_buf[lev]->nComp(), ng));
                // pmf->ParallelCopy(*charge_buf[lev][idim], 0, 0, charge_buf[lev]->nComp(), ng, ng);
                charge_buf[lev] = std::move(pmf);
            }
            if (current_buffer_masks[lev])
            {
                const IntVect& ng = current_buffer_masks[lev]->nGrowVect();
                auto pmf = std::unique_ptr<iMultiFab>(new iMultiFab(current_buffer_masks[lev]->boxArray(),
                                                                    dm, current_buffer_masks[lev]->nComp(), ng));
                // pmf->ParallelCopy(*current_buffer_masks[lev], 0, 0, current_buffer_masks[lev]->nComp(), ng, ng);
                current_buffer_masks[lev] = std::move(pmf);
            }
            if (gather_buffer_masks[lev])
            {
                const IntVect& ng = gather_buffer_masks[lev]->nGrowVect();
                auto pmf = std::unique_ptr<iMultiFab>(new iMultiFab(gather_buffer_masks[lev]->boxArray(),
                                                                    dm, gather_buffer_masks[lev]->nComp(), ng));
                // pmf->ParallelCopy(*gather_buffer_masks[lev], 0, 0, gather_buffer_masks[lev]->nComp(), ng, ng);
                gather_buffer_masks[lev] = std::move(pmf);
            }
        }

        if (costs[lev] != nullptr)
        {
            costs[lev].reset(new amrex::Vector<Real>);
            const int nboxes = Efield_fp[lev][0].get()->size();
            costs[lev]->resize(nboxes, 0.0);
        }

        SetDistributionMap(lev, dm);

    } else
    {
        amrex::Abort("RemakeLevel: to be implemented");
    }
}

void
WarpX::ComputeCostsHeuristic (amrex::Vector<std::unique_ptr<amrex::Vector<amrex::Real> > >& a_costs)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto & mypc_ref = WarpX::GetInstance().GetPartContainer();
        auto nSpecies = mypc_ref.nSpecies();

        // Species loop
        for (int i_s = 0; i_s < nSpecies; ++i_s)
        {
            auto & myspc = mypc_ref.GetParticleContainer(i_s);

            // Particle loop
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
            {
                (*a_costs[lev])[pti.index()] += costs_heuristic_particles_wt*pti.numParticles();
            }
        }

        //Cell loop
        MultiFab* Ex = Efield_fp[lev][0].get();
        for (MFIter mfi(*Ex, false); mfi.isValid(); ++mfi)
        {
            const Box& gbx = mfi.growntilebox();
            (*a_costs[lev])[mfi.index()] += costs_heuristic_cells_wt*gbx.numPts();
        }
    }
}

void
WarpX::ResetCosts ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        costs[lev]->assign((*costs[lev]).size(), 0.0);
    }
}
