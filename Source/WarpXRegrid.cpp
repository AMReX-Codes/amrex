
#include <WarpX.H>
#include <AMReX_BLProfiler.H>

using namespace amrex;

void
WarpX::LoadBalance ()
{
    BL_PROFILE_REGION("LoadBalance");
    BL_PROFILE("WarpX::LoadBalance()");

    AMREX_ALWAYS_ASSERT(costs[0] != nullptr);

    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        const Real nboxes = costs[lev]->size();
        const Real nprocs = ParallelDescriptor::NProcs();
        const int nmax = static_cast<int>(std::ceil(nboxes/nprocs*load_balance_knapsack_factor));
        const DistributionMapping newdm = (load_balance_with_sfc)
	  ? DistributionMapping::makeSFC(*costs[lev], false)
            : DistributionMapping::makeKnapSack(*costs[lev], nmax);
        RemakeLevel(lev, t_new[lev], boxArray(lev), newdm);
    }

    mypc->Redistribute();
}

void
WarpX::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    if (ba == boxArray(lev))
    {
        if (ParallelDescriptor::NProcs() == 1) return;

#ifdef WARPX_DO_ELECTROSTATIC        
        AMREX_ALWAYS_ASSERT(masks[lev] == nullptr);
        AMREX_ALWAYS_ASSERT(gather_masks[lev] == nullptr);
#endif // WARPX_DO_ELECTROSTATIC
        
        // Fine patch

        for (int idim=0; idim < 3; ++idim)
        {
            {
                const IntVect& ng = Bfield_fp[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_fp[lev][idim]->boxArray(),
                                                                  dm, 1, ng));
                pmf->Redistribute(*Bfield_fp[lev][idim], 0, 0, 1, ng);
                Bfield_fp[lev][idim] = std::move(pmf);
            }
            {
                const IntVect& ng = Efield_fp[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_fp[lev][idim]->boxArray(),
                                                                  dm, 1, ng));
                pmf->Redistribute(*Efield_fp[lev][idim], 0, 0, 1, ng);
                Efield_fp[lev][idim] = std::move(pmf);
            }
            {
                const IntVect& ng = current_fp[lev][idim]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_fp[lev][idim]->boxArray(),
                                                                  dm, 1, ng));
                // pmf->Redistribute(*current_fp[lev][idim], 0, 0, 1, ng);
                current_fp[lev][idim] = std::move(pmf);
            }
        }

        if (F_fp[lev] != nullptr) {
            const IntVect& ng = F_fp[lev]->nGrowVect();
            auto pmf = std::unique_ptr<MultiFab>(new MultiFab(F_fp[lev]->boxArray(),
                                                              dm, 1, ng));
            pmf->Redistribute(*F_fp[lev], 0, 0, 1, ng);
            F_fp[lev] = std::move(pmf);
        }

        if (rho_fp[lev] != nullptr) {
            const int nc = rho_fp[lev]->nComp();
            const IntVect& ng = rho_fp[lev]->nGrowVect();
            auto pmf = std::unique_ptr<MultiFab>(new MultiFab(rho_fp[lev]->boxArray(),
                                                              dm, nc, ng));
            // pmf->Redistribute(*rho_fp[lev], 0, 0, nc, ng);
            rho_fp[lev] = std::move(pmf);
        }

        // Aux patch

        if (lev == 0)
        {
            for (int idim = 0; idim < 3; ++idim) {
                Bfield_aux[lev][idim].reset(new MultiFab(*Bfield_fp[lev][idim], amrex::make_alias, 0, 1));
                Efield_aux[lev][idim].reset(new MultiFab(*Efield_fp[lev][idim], amrex::make_alias, 0, 1));
            }
        } else {
            for (int idim=0; idim < 3; ++idim)
            {
                {
                    const IntVect& ng = Bfield_aux[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_aux[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->Redistribute(*Bfield_aux[lev][idim], 0, 0, 1, ng);
                    Bfield_aux[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = Efield_aux[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_aux[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->Redistribute(*Efield_aux[lev][idim], 0, 0, 1, ng);
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
                                                                      dm, 1, ng));
                    pmf->Redistribute(*Bfield_cp[lev][idim], 0, 0, 1, ng);
                    Bfield_cp[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = Efield_cp[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_cp[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    pmf->Redistribute(*Efield_cp[lev][idim], 0, 0, 1, ng);
                    Efield_cp[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = current_cp[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_cp[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->Redistribute(*current_cp[lev][idim], 0, 0, 1, ng);
                    current_cp[lev][idim] = std::move(pmf);
                }
            }

            if (F_cp[lev] != nullptr) {
                const IntVect& ng = F_cp[lev]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(F_cp[lev]->boxArray(),
                                                                  dm, 1, ng));
                pmf->Redistribute(*F_cp[lev], 0, 0, 1, ng);
                F_cp[lev] = std::move(pmf);
            }

            if (rho_cp[lev] != nullptr) {
                const int nc = rho_cp[lev]->nComp();
                const IntVect& ng = rho_cp[lev]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(rho_cp[lev]->boxArray(),
                                                                  dm, nc, ng));
                // pmf->Redistribute(*rho_cp[lev], 0, 0, nc, ng);
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
                                                                      dm, 1, ng));
                    // pmf->ParallelCopy(*Bfield_cax[lev][idim], 0, 0, 1, ng, ng);
                    Bfield_cax[lev][idim] = std::move(pmf);
                }
                if (Efield_cax[lev][idim])
                {
                    const IntVect& ng = Efield_cax[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_cax[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->ParallelCopy(*Efield_cax[lev][idim], 0, 0, 1, ng, ng);
                    Efield_cax[lev][idim] = std::move(pmf);
                }
                if (current_buf[lev][idim])
                {
                    const IntVect& ng = current_buf[lev][idim]->nGrowVect();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_buf[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->ParallelCopy(*current_buf[lev][idim], 0, 0, 1, ng, ng);
                    current_buf[lev][idim] = std::move(pmf);
                }
            }
            if (charge_buf[lev])
            {
                const IntVect& ng = charge_buf[lev]->nGrowVect();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(charge_buf[lev]->boxArray(),
                                                                  dm, 1, ng));
                // pmf->ParallelCopy(*charge_buf[lev][idim], 0, 0, 1, ng, ng);
                charge_buf[lev] = std::move(pmf);
            }
            if (current_buffer_masks[lev])
            {
                const IntVect& ng = current_buffer_masks[lev]->nGrowVect();
                auto pmf = std::unique_ptr<iMultiFab>(new iMultiFab(current_buffer_masks[lev]->boxArray(),
                                                                    dm, 1, ng));
                // pmf->ParallelCopy(*current_buffer_masks[lev], 0, 0, 1, ng, ng);
                current_buffer_masks[lev] = std::move(pmf);                
            }
            if (gather_buffer_masks[lev])
            {
                const IntVect& ng = gather_buffer_masks[lev]->nGrowVect();
                auto pmf = std::unique_ptr<iMultiFab>(new iMultiFab(gather_buffer_masks[lev]->boxArray(),
                                                                    dm, 1, ng));
                // pmf->ParallelCopy(*gather_buffer_masks[lev], 0, 0, 1, ng, ng);
                gather_buffer_masks[lev] = std::move(pmf);                
            }
        }

        if (costs[lev] != nullptr) {
            costs[lev].reset(new MultiFab(costs[lev]->boxArray(), dm, 1, 0));
            costs[lev]->setVal(0.0);
        }

        SetDistributionMap(lev, dm);
    }
    else
    {
        amrex::Abort("RemakeLevel: to be implemented");
    }
}
