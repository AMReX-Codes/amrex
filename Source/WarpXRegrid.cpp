
#include <WarpX.H>
#include <AMReX_BLProfiler.H>

using namespace amrex;

void
WarpX::LoadBalance ()
{
    BL_PROFILE("WarpX::LoadBalance()");

    AMREX_ALWAYS_ASSERT(costs[0] != nullptr);

    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        const DistributionMapping& newdm = DistributionMapping::makeKnapSack(*costs[lev]);
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

        AMREX_ALWAYS_ASSERT(masks[lev] == nullptr);
        AMREX_ALWAYS_ASSERT(gather_masks[lev] == nullptr);

        // Fine patch

        for (int idim=0; idim < 3; ++idim)
        {
            {
                int ng = Bfield_fp[lev][idim]->nGrow();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_fp[lev][idim]->boxArray(),
                                                                  dm, 1, ng));
                pmf->ParallelCopy(*Bfield_fp[lev][idim], 0, 0, 1, ng, ng);
                Bfield_fp[lev][idim] = std::move(pmf);
            }
            {
                int ng = Efield_fp[lev][idim]->nGrow();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_fp[lev][idim]->boxArray(),
                                                                  dm, 1, ng));
                pmf->ParallelCopy(*Efield_fp[lev][idim], 0, 0, 1, ng, ng);
                Efield_fp[lev][idim] = std::move(pmf);
            }
            {
                int ng = current_fp[lev][idim]->nGrow();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_fp[lev][idim]->boxArray(),
                                                                  dm, 1, ng));
                // pmf->ParallelCopy(*current_fp[lev][idim], 0, 0, 1, ng, ng);
                current_fp[lev][idim] = std::move(pmf);
            }
        }

        if (F_fp[lev] != nullptr) {
            int ng = F_fp[lev]->nGrow();
            auto pmf = std::unique_ptr<MultiFab>(new MultiFab(F_fp[lev]->boxArray(),
                                                              dm, 1, ng));
            pmf->ParallelCopy(*F_fp[lev], 0, 0, 1, ng, ng);
            F_fp[lev] = std::move(pmf);            
        }

        if (rho_fp[lev] != nullptr) {
            int ng = rho_fp[lev]->nGrow();
            auto pmf = std::unique_ptr<MultiFab>(new MultiFab(rho_fp[lev]->boxArray(),
                                                              dm, 1, ng));
            // pmf->ParallelCopy(*rho_fp[lev], 0, 0, 1, ng, ng);
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
                    int ng = Bfield_aux[lev][idim]->nGrow();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_aux[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->ParallelCopy(*Bfield_aux[lev][idim], 0, 0, 1, ng, ng);
                    Bfield_aux[lev][idim] = std::move(pmf);
                }
                {
                    int ng = Efield_aux[lev][idim]->nGrow();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_aux[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->ParallelCopy(*Efield_aux[lev][idim], 0, 0, 1, ng, ng);
                    Efield_aux[lev][idim] = std::move(pmf);
                }
            }
        }

        // Coarse patch
        if (lev > 0) {
            for (int idim=0; idim < 3; ++idim)
            {
                {
                    int ng = Bfield_cp[lev][idim]->nGrow();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Bfield_cp[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    pmf->ParallelCopy(*Bfield_cp[lev][idim], 0, 0, 1, ng, ng);
                    Bfield_cp[lev][idim] = std::move(pmf);
                }
                {
                    int ng = Efield_cp[lev][idim]->nGrow();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(Efield_cp[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    pmf->ParallelCopy(*Efield_cp[lev][idim], 0, 0, 1, ng, ng);
                    Efield_cp[lev][idim] = std::move(pmf);
                }
                {
                    int ng = current_cp[lev][idim]->nGrow();
                    auto pmf = std::unique_ptr<MultiFab>(new MultiFab(current_cp[lev][idim]->boxArray(),
                                                                      dm, 1, ng));
                    // pmf->ParallelCopy(*current_cp[lev][idim], 0, 0, 1, ng, ng);
                    current_cp[lev][idim] = std::move(pmf);
                }
            }
            
            if (F_cp[lev] != nullptr) {
                int ng = F_cp[lev]->nGrow();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(F_cp[lev]->boxArray(),
                                                                  dm, 1, ng));
                pmf->ParallelCopy(*F_cp[lev], 0, 0, 1, ng, ng);
                F_cp[lev] = std::move(pmf);            
            }
            
            if (rho_cp[lev] != nullptr) {
                int ng = rho_cp[lev]->nGrow();
                auto pmf = std::unique_ptr<MultiFab>(new MultiFab(rho_cp[lev]->boxArray(),
                                                                  dm, 1, ng));
                // pmf->ParallelCopy(*rho_cp[lev], 0, 0, 1, ng, ng);
                rho_cp[lev] = std::move(pmf);            
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

