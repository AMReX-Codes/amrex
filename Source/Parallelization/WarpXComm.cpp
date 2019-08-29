#include <WarpX.H>
#include <WarpX_f.H>
#include <WarpXSumGuardCells.H>

#include <AMReX_FillPatchUtil_F.H>

#include <algorithm>
#include <cstdlib>

using namespace amrex;

void
WarpX::ExchangeWithPmlB (int lev)
{
    if (do_pml && pml[lev]->ok()) {
        pml[lev]->ExchangeB({ Bfield_fp[lev][0].get(),
                              Bfield_fp[lev][1].get(),
                              Bfield_fp[lev][2].get() },
                            { Bfield_cp[lev][0].get(),
                              Bfield_cp[lev][1].get(),
                              Bfield_cp[lev][2].get() },
                              do_pml_in_domain);
    }
}

void
WarpX::ExchangeWithPmlE (int lev)
{
    if (do_pml && pml[lev]->ok()) {
        pml[lev]->ExchangeE({ Efield_fp[lev][0].get(),
                              Efield_fp[lev][1].get(),
                              Efield_fp[lev][2].get() },
                            { Efield_cp[lev][0].get(),
                              Efield_cp[lev][1].get(),
                              Efield_cp[lev][2].get() },
                              do_pml_in_domain);
    }
}

void
WarpX::ExchangeWithPmlF (int lev)
{
    if (do_pml && pml[lev]->ok()) {
        pml[lev]->ExchangeF(F_fp[lev].get(),
                            F_cp[lev].get(),
                            do_pml_in_domain);
    }
}

void
WarpX::UpdateAuxilaryData ()
{
    BL_PROFILE("UpdateAuxilaryData()");

    const int use_limiter = 0;

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& crse_period = Geom(lev-1).periodicity();
        const IntVect& ng = Bfield_cp[lev][0]->nGrowVect();
        const DistributionMapping& dm = Bfield_cp[lev][0]->DistributionMap();

        // B field
        {
            MultiFab dBx(Bfield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dBy(Bfield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dBz(Bfield_cp[lev][2]->boxArray(), dm, 1, ng);
            dBx.setVal(0.0);
            dBy.setVal(0.0);
            dBz.setVal(0.0);
            dBx.ParallelCopy(*Bfield_aux[lev-1][0], 0, 0, 1, ng, ng, crse_period);
            dBy.ParallelCopy(*Bfield_aux[lev-1][1], 0, 0, 1, ng, ng, crse_period);
            dBz.ParallelCopy(*Bfield_aux[lev-1][2], 0, 0, 1, ng, ng, crse_period);
            if (Bfield_cax[lev][0])
            {
                MultiFab::Copy(*Bfield_cax[lev][0], dBx, 0, 0, 1, ng);
                MultiFab::Copy(*Bfield_cax[lev][1], dBy, 0, 0, 1, ng);
                MultiFab::Copy(*Bfield_cax[lev][2], dBz, 0, 0, 1, ng);
            }
            MultiFab::Subtract(dBx, *Bfield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dBy, *Bfield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dBz, *Bfield_cp[lev][2], 0, 0, 1, ng);

            const Real* dx = Geom(lev-1).CellSize();
            const int refinement_ratio = refRatio(lev-1)[0];
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                std::array<FArrayBox,3> bfab;
                for (MFIter mfi(*Bfield_aux[lev][0]); mfi.isValid(); ++mfi)
                {
                    Box ccbx = mfi.fabbox();
                    ccbx.enclosedCells();
                    ccbx.coarsen(refinement_ratio).refine(refinement_ratio); // so that ccbx is coarsenable

                    const FArrayBox& cxfab = dBx[mfi];
                    const FArrayBox& cyfab = dBy[mfi];
                    const FArrayBox& czfab = dBz[mfi];
                    bfab[0].resize(amrex::convert(ccbx,Bx_nodal_flag));
                    bfab[1].resize(amrex::convert(ccbx,By_nodal_flag));
                    bfab[2].resize(amrex::convert(ccbx,Bz_nodal_flag));

#if (AMREX_SPACEDIM == 3)
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(bfab[0]),
                                                 BL_TO_FORTRAN_ANYD(bfab[1]),
                                                 BL_TO_FORTRAN_ANYD(bfab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(cyfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &refinement_ratio,&use_limiter);
#else
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(bfab[0]),
                                                 BL_TO_FORTRAN_ANYD(bfab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &refinement_ratio,&use_limiter);
                    amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(bfab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &refinement_ratio,&use_limiter);
#endif

                    for (int idim = 0; idim < 3; ++idim)
                    {
                        FArrayBox& aux = (*Bfield_aux[lev][idim])[mfi];
                        FArrayBox& fp  =  (*Bfield_fp[lev][idim])[mfi];
                        const Box& bx = aux.box();
                        aux.copy(fp, bx, 0, bx, 0, 1);
                        aux.plus(bfab[idim], bx, bx, 0, 0, 1);
                    }
                }
            }
        }

        // E field
        {
            MultiFab dEx(Efield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dEy(Efield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dEz(Efield_cp[lev][2]->boxArray(), dm, 1, ng);
            dEx.setVal(0.0);
            dEy.setVal(0.0);
            dEz.setVal(0.0);
            dEx.ParallelCopy(*Efield_aux[lev-1][0], 0, 0, 1, ng, ng, crse_period);
            dEy.ParallelCopy(*Efield_aux[lev-1][1], 0, 0, 1, ng, ng, crse_period);
            dEz.ParallelCopy(*Efield_aux[lev-1][2], 0, 0, 1, ng, ng, crse_period);
            if (Efield_cax[lev][0])
            {
                MultiFab::Copy(*Efield_cax[lev][0], dEx, 0, 0, 1, ng);
                MultiFab::Copy(*Efield_cax[lev][1], dEy, 0, 0, 1, ng);
                MultiFab::Copy(*Efield_cax[lev][2], dEz, 0, 0, 1, ng);
            }
            MultiFab::Subtract(dEx, *Efield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dEy, *Efield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dEz, *Efield_cp[lev][2], 0, 0, 1, ng);

            const int refinement_ratio = refRatio(lev-1)[0];
#ifdef _OPEMP
#pragma omp parallel
#endif
            {
                std::array<FArrayBox,3> efab;
                for (MFIter mfi(*Efield_aux[lev][0]); mfi.isValid(); ++mfi)
                {
                    Box ccbx = mfi.fabbox();
                    ccbx.enclosedCells();
                    ccbx.coarsen(refinement_ratio).refine(refinement_ratio); // so that ccbx is coarsenable

                    const FArrayBox& cxfab = dEx[mfi];
                    const FArrayBox& cyfab = dEy[mfi];
                    const FArrayBox& czfab = dEz[mfi];
                    efab[0].resize(amrex::convert(ccbx,Ex_nodal_flag));
                    efab[1].resize(amrex::convert(ccbx,Ey_nodal_flag));
                    efab[2].resize(amrex::convert(ccbx,Ez_nodal_flag));

#if (AMREX_SPACEDIM == 3)
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(efab[0]),
                                        BL_TO_FORTRAN_ANYD(efab[1]),
                                        BL_TO_FORTRAN_ANYD(efab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(cyfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &refinement_ratio,&use_limiter);
#else
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(efab[0]),
                                        BL_TO_FORTRAN_ANYD(efab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &refinement_ratio,&use_limiter);
                    amrex_interp_nd_efield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(efab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &refinement_ratio);
#endif

                    for (int idim = 0; idim < 3; ++idim)
                    {
                        FArrayBox& aux = (*Efield_aux[lev][idim])[mfi];
                        FArrayBox& fp  =  (*Efield_fp[lev][idim])[mfi];
                        const Box& bx = aux.box();
                        aux.copy(fp, bx, 0, bx, 0, 1);
                        aux.plus(efab[idim], bx, bx, 0, 0, 1);
                    }
                }
            }
        }
    }
}

void
WarpX::FillBoundaryB ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryB(lev);
    }
}

void
WarpX::FillBoundaryE ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryE(lev);
    }
}

void
WarpX::FillBoundaryF ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryF(lev);
    }
}

void
WarpX::FillBoundaryE(int lev)
{
    FillBoundaryE(lev, PatchType::fine);
    if (lev > 0) FillBoundaryE(lev, PatchType::coarse);
}

void
WarpX::FillBoundaryE (int lev, PatchType patch_type)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev]->ok())
        {
    	    pml[lev]->ExchangeE(patch_type,
                              { Efield_fp[lev][0].get(),
                                Efield_fp[lev][1].get(),
                                Efield_fp[lev][2].get() },
                                do_pml_in_domain);
    	    pml[lev]->FillBoundaryE(patch_type);
        }

        const auto& period = Geom(lev).periodicity();
        Vector<MultiFab*> mf{Efield_fp[lev][0].get(),Efield_fp[lev][1].get(),Efield_fp[lev][2].get()};
        amrex::FillBoundary(mf, period);
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
        {
        pml[lev]->ExchangeE(patch_type,
                            { Efield_cp[lev][0].get(),
                              Efield_cp[lev][1].get(),
                              Efield_cp[lev][2].get() },
                              do_pml_in_domain);
        pml[lev]->FillBoundaryE(patch_type);
        }

        const auto& cperiod = Geom(lev-1).periodicity();
        Vector<MultiFab*> mf{Efield_cp[lev][0].get(),Efield_cp[lev][1].get(),Efield_cp[lev][2].get()};
        amrex::FillBoundary(mf, cperiod);
    }
}

void
WarpX::FillBoundaryB (int lev)
{
    FillBoundaryB(lev, PatchType::fine);
    if (lev > 0) FillBoundaryB(lev, PatchType::coarse);
}

void
WarpX::FillBoundaryB (int lev, PatchType patch_type)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev]->ok())
        {
	    pml[lev]->ExchangeB(patch_type,
                            { Bfield_fp[lev][0].get(),
                              Bfield_fp[lev][1].get(),
                              Bfield_fp[lev][2].get() },
                              do_pml_in_domain);
        pml[lev]->FillBoundaryB(patch_type);
        }
        const auto& period = Geom(lev).periodicity();
        Vector<MultiFab*> mf{Bfield_fp[lev][0].get(),Bfield_fp[lev][1].get(),Bfield_fp[lev][2].get()};
        amrex::FillBoundary(mf, period);
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
        {
        pml[lev]->ExchangeB(patch_type,
                      { Bfield_cp[lev][0].get(),
                        Bfield_cp[lev][1].get(),
                        Bfield_cp[lev][2].get() },
                        do_pml_in_domain);
        pml[lev]->FillBoundaryB(patch_type);
        }
        const auto& cperiod = Geom(lev-1).periodicity();
        Vector<MultiFab*> mf{Bfield_cp[lev][0].get(),Bfield_cp[lev][1].get(),Bfield_cp[lev][2].get()};
        amrex::FillBoundary(mf, cperiod);
    }
}

void
WarpX::FillBoundaryF (int lev)
{
  FillBoundaryF(lev, PatchType::fine);
  if (lev > 0) FillBoundaryF(lev, PatchType::coarse);
}

void
WarpX::FillBoundaryF (int lev, PatchType patch_type)
{
    if (patch_type == PatchType::fine && F_fp[lev])
    {
        if (do_pml && pml[lev]->ok())
        {
            pml[lev]->ExchangeF(patch_type, F_fp[lev].get(),
                                do_pml_in_domain);
            pml[lev]->FillBoundaryF(patch_type);
        }

        const auto& period = Geom(lev).periodicity();
        F_fp[lev]->FillBoundary(period);
    }
    else if (patch_type == PatchType::coarse && F_cp[lev])
    {
        if (do_pml && pml[lev]->ok())
        {
        pml[lev]->ExchangeF(patch_type, F_cp[lev].get(),
                            do_pml_in_domain);
        pml[lev]->FillBoundaryF(patch_type);
        }

        const auto& cperiod = Geom(lev-1).periodicity();
        F_cp[lev]->FillBoundary(cperiod);
    }
}

void
WarpX::SyncCurrent ()
{
    BL_PROFILE("SyncCurrent()");

    // Restrict fine patch current onto the coarse patch, before
    // summing the guard cells of the fine patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        current_cp[lev][0]->setVal(0.0);
        current_cp[lev][1]->setVal(0.0);
        current_cp[lev][2]->setVal(0.0);

        const IntVect& refinement_ratio = refRatio(lev-1);

        std::array<const MultiFab*,3> fine { current_fp[lev][0].get(),
                                             current_fp[lev][1].get(),
                                             current_fp[lev][2].get() };
        std::array<      MultiFab*,3> crse { current_cp[lev][0].get(),
                                             current_cp[lev][1].get(),
                                             current_cp[lev][2].get() };
        SyncCurrent(fine, crse, refinement_ratio[0]);
    }

    // For each level
    // - apply filter to the coarse patch/buffer of `lev+1` and fine patch of `lev` (same resolution)
    // - add the coarse patch/buffer of `lev+1` into the fine patch of `lev`
    // - sum guard cells of the coarse patch of `lev+1` and fine patch of `lev`
    for (int lev=0; lev <= finest_level; ++lev) {
        AddCurrentFromFineLevelandSumBoundary(lev);
    }
}

/** \brief Fills the values of the current on the coarse patch by
 *  averaging the values of the current of the fine patch (on the same level).
 */
void
WarpX::SyncCurrent (const std::array<const amrex::MultiFab*,3>& fine,
                    const std::array<      amrex::MultiFab*,3>& crse,
                    int refinement_ratio)
{
    BL_ASSERT(refinement_ratio == 2);
    const IntVect& ng = (fine[0]->nGrowVect() + 1) /refinement_ratio;

#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        FArrayBox ffab;
        for (int idim = 0; idim < 3; ++idim)
        {
            for (MFIter mfi(*crse[idim],true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(ng);
                Box fbx = amrex::grow(amrex::refine(bx,refinement_ratio),1);
                ffab.resize(fbx);
                fbx &= (*fine[idim])[mfi].box();
                ffab.setVal(0.0);
                ffab.copy((*fine[idim])[mfi], fbx, 0, fbx, 0, 1);
                WRPX_SYNC_CURRENT(bx.loVect(), bx.hiVect(),
                                   BL_TO_FORTRAN_ANYD((*crse[idim])[mfi]),
                                   BL_TO_FORTRAN_ANYD(ffab),
                                   &idim);
            }
        }
    }
}

void
WarpX::SyncRho ()
{
    if (!rho_fp[0]) return;
    const int ncomp = rho_fp[0]->nComp();

    // Restrict fine patch onto the coarse patch,
    // before summing the guard cells of the fine patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        rho_cp[lev]->setVal(0.0);
        const IntVect& refinement_ratio = refRatio(lev-1);
        SyncRho(*rho_fp[lev], *rho_cp[lev], refinement_ratio[0]);
    }

    // For each level
    // - apply filter to the coarse patch/buffer of `lev+1` and fine patch of `lev` (same resolution)
    // - add the coarse patch/buffer of `lev+1` into the fine patch of `lev`
    // - sum guard cells of the coarse patch of `lev+1` and fine patch of `lev`
    for (int lev=0; lev <= finest_level; ++lev) {
        AddRhoFromFineLevelandSumBoundary(lev, 0, ncomp);
    }
}

/** \brief Fills the values of the charge density on the coarse patch by
 *  averaging the values of the charge density of the fine patch (on the same level).
 */
void
WarpX::SyncRho (const MultiFab& fine, MultiFab& crse, int refinement_ratio)
{
    BL_ASSERT(refinement_ratio == 2);
    const IntVect& ng = (fine.nGrowVect()+1)/refinement_ratio;
    const int nc = fine.nComp();

#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        FArrayBox ffab;
        for (MFIter mfi(crse,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(ng);
            Box fbx = amrex::grow(amrex::refine(bx,refinement_ratio),1);
            ffab.resize(fbx, nc);
            fbx &= fine[mfi].box();
            ffab.setVal(0.0);
            ffab.copy(fine[mfi], fbx, 0, fbx, 0, nc);
            WRPX_SYNC_RHO(bx.loVect(), bx.hiVect(),
                          BL_TO_FORTRAN_ANYD(crse[mfi]),
                          BL_TO_FORTRAN_ANYD(ffab),
                          &nc);
        }
    }
}

/** \brief Fills the values of the current on the coarse patch by
 *  averaging the values of the current of the fine patch (on the same level).
 */
void
WarpX::RestrictCurrentFromFineToCoarsePatch (int lev)
{
    current_cp[lev][0]->setVal(0.0);
    current_cp[lev][1]->setVal(0.0);
    current_cp[lev][2]->setVal(0.0);

    const IntVect& refinement_ratio = refRatio(lev-1);

    std::array<const MultiFab*,3> fine { current_fp[lev][0].get(),
                                         current_fp[lev][1].get(),
                                         current_fp[lev][2].get() };
    std::array<      MultiFab*,3> crse { current_cp[lev][0].get(),
                                         current_cp[lev][1].get(),
                                         current_cp[lev][2].get() };
    SyncCurrent(fine, crse, refinement_ratio[0]);
}

void
WarpX::ApplyFilterandSumBoundaryJ (int lev, PatchType patch_type)
{
    const int glev = (patch_type == PatchType::fine) ? lev : lev-1;
    const auto& period = Geom(glev).periodicity();
    auto& j = (patch_type == PatchType::fine) ? current_fp[lev] : current_cp[lev];
    for (int idim = 0; idim < 3; ++idim) {
        if (use_filter) {
            IntVect ng = j[idim]->nGrowVect();
            ng += bilinear_filter.stencil_length_each_dir-1;
            MultiFab jf(j[idim]->boxArray(), j[idim]->DistributionMap(), 1, ng);
            bilinear_filter.ApplyStencil(jf, *j[idim]);
            WarpXSumGuardCells(*(j[idim]), jf, period);
        } else {
            WarpXSumGuardCells(*(j[idim]), period);
        }
    }
}

/* /brief Update the currents of `lev` by adding the currents from particles
*         that are in the mesh refinement patches at `lev+1`
*
* More precisely, apply filter and sum boundaries for the current of:
* - the fine patch of `lev`
* - the coarse patch of `lev+1` (same resolution)
* - the buffer regions of the coarse patch of `lev+1` (i.e. for particules
* that are within the mesh refinement patch, but do not deposit on the
* mesh refinement patch because they are too close to the boundary)
*
* Then update the fine patch of `lev` by adding the currents for the coarse
* patch (and buffer region) of `lev+1`
*/
void
WarpX::AddCurrentFromFineLevelandSumBoundary (int lev)
{
    ApplyFilterandSumBoundaryJ(lev, PatchType::fine);

    if (lev < finest_level) {
        // When there are current buffers, unlike coarse patch,
        // we don't care about the final state of them.

        const auto& period = Geom(lev).periodicity();
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab mf(current_fp[lev][idim]->boxArray(),
                        current_fp[lev][idim]->DistributionMap(), 1, 0);
            mf.setVal(0.0);
            if (use_filter && current_buf[lev+1][idim])
            {
                // coarse patch of fine level
                IntVect ng = current_cp[lev+1][idim]->nGrowVect();
                ng += bilinear_filter.stencil_length_each_dir-1;
                MultiFab jfc(current_cp[lev+1][idim]->boxArray(),
                             current_cp[lev+1][idim]->DistributionMap(), 1, ng);
                bilinear_filter.ApplyStencil(jfc, *current_cp[lev+1][idim]);

                // buffer patch of fine level
                MultiFab jfb(current_buf[lev+1][idim]->boxArray(),
                             current_buf[lev+1][idim]->DistributionMap(), 1, ng);
                bilinear_filter.ApplyStencil(jfb, *current_buf[lev+1][idim]);

                MultiFab::Add(jfb, jfc, 0, 0, 1, ng);
                mf.ParallelAdd(jfb, 0, 0, 1, ng, IntVect::TheZeroVector(), period);

                WarpXSumGuardCells(*current_cp[lev+1][idim], jfc, period);
            }
            else if (use_filter) // but no buffer
            {
                // coarse patch of fine level
                IntVect ng = current_cp[lev+1][idim]->nGrowVect();
                ng += bilinear_filter.stencil_length_each_dir-1;
                MultiFab jf(current_cp[lev+1][idim]->boxArray(),
                            current_cp[lev+1][idim]->DistributionMap(), 1, ng);
                bilinear_filter.ApplyStencil(jf, *current_cp[lev+1][idim]);
                mf.ParallelAdd(jf, 0, 0, 1, ng, IntVect::TheZeroVector(), period);
                WarpXSumGuardCells(*current_cp[lev+1][idim], jf, period);
            }
            else if (current_buf[lev+1][idim]) // but no filter
            {
                MultiFab::Copy(*current_buf[lev+1][idim],
                               *current_cp [lev+1][idim], 0, 0, 1,
                               current_cp[lev+1][idim]->nGrow());
                mf.ParallelAdd(*current_buf[lev+1][idim], 0, 0, 1,
                               current_buf[lev+1][idim]->nGrowVect(), IntVect::TheZeroVector(),
                               period);
                WarpXSumGuardCells(*(current_cp[lev+1][idim]), period);
            }
            else // no filter, no buffer
            {
                mf.ParallelAdd(*current_cp[lev+1][idim], 0, 0, 1,
                               current_cp[lev+1][idim]->nGrowVect(), IntVect::TheZeroVector(),
                               period);
                WarpXSumGuardCells(*(current_cp[lev+1][idim]), period);
            }
            MultiFab::Add(*current_fp[lev][idim], mf, 0, 0, 1, 0);
        }
        NodalSyncJ(lev+1, PatchType::coarse);
    }
    NodalSyncJ(lev, PatchType::fine);
}

void
WarpX::RestrictRhoFromFineToCoarsePatch (int lev)
{
    if (rho_fp[lev]) {
        rho_cp[lev]->setVal(0.0);
        const IntVect& refinement_ratio = refRatio(lev-1);
        SyncRho(*rho_fp[lev], *rho_cp[lev], refinement_ratio[0]);
    }
}

void
WarpX::ApplyFilterandSumBoundaryRho (int lev, PatchType patch_type, int icomp, int ncomp)
{
    const int glev = (patch_type == PatchType::fine) ? lev : lev-1;
    const auto& period = Geom(glev).periodicity();
    auto& r = (patch_type == PatchType::fine) ? rho_fp[lev] : rho_cp[lev];
    if (r == nullptr) return;
    if (use_filter) {
        IntVect ng = r->nGrowVect();
        ng += bilinear_filter.stencil_length_each_dir-1;
        MultiFab rf(r->boxArray(), r->DistributionMap(), ncomp, ng);
        bilinear_filter.ApplyStencil(rf, *r, icomp, 0, ncomp);
        WarpXSumGuardCells(*r, rf, period, icomp, ncomp );
    } else {
        WarpXSumGuardCells(*r, period, icomp, ncomp);
    }
}

/* /brief Update the charge density of `lev` by adding the charge density from particles
*         that are in the mesh refinement patches at `lev+1`
*
* More precisely, apply filter and sum boundaries for the charge density of:
* - the fine patch of `lev`
* - the coarse patch of `lev+1` (same resolution)
* - the buffer regions of the coarse patch of `lev+1` (i.e. for particules
* that are within the mesh refinement patch, but do not deposit on the
* mesh refinement patch because they are too close to the boundary)
*
* Then update the fine patch of `lev` by adding the charge density for the coarse
* patch (and buffer region) of `lev+1`
*/
void
WarpX::AddRhoFromFineLevelandSumBoundary(int lev, int icomp, int ncomp)
{
    if (!rho_fp[lev]) return;

    ApplyFilterandSumBoundaryRho(lev, PatchType::fine, icomp, ncomp);

    if (lev < finest_level){

        const auto& period = Geom(lev).periodicity();
        MultiFab mf(rho_fp[lev]->boxArray(),
                    rho_fp[lev]->DistributionMap(),
                    ncomp, 0);
        mf.setVal(0.0);
        if (use_filter && charge_buf[lev+1])
        {
            // coarse patch of fine level
            IntVect ng = rho_cp[lev+1]->nGrowVect();
            ng += bilinear_filter.stencil_length_each_dir-1;
            MultiFab rhofc(rho_cp[lev+1]->boxArray(),
                         rho_cp[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rhofc, *rho_cp[lev+1], icomp, 0, ncomp);

            // buffer patch of fine level
            MultiFab rhofb(charge_buf[lev+1]->boxArray(),
                           charge_buf[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rhofb, *charge_buf[lev+1], icomp, 0, ncomp);

            MultiFab::Add(rhofb, rhofc, 0, 0, ncomp, ng);
            mf.ParallelAdd(rhofb, 0, 0, ncomp, ng, IntVect::TheZeroVector(), period);
            WarpXSumGuardCells( *rho_cp[lev+1], rhofc, period, icomp, ncomp );
        }
        else if (use_filter) // but no buffer
        {
            IntVect ng = rho_cp[lev+1]->nGrowVect();
            ng += bilinear_filter.stencil_length_each_dir-1;
            MultiFab rf(rho_cp[lev+1]->boxArray(), rho_cp[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rf, *rho_cp[lev+1], icomp, 0, ncomp);
            mf.ParallelAdd(rf, 0, 0, ncomp, ng, IntVect::TheZeroVector(), period);
            WarpXSumGuardCells( *rho_cp[lev+1], rf, period, icomp, ncomp );
        }
        else if (charge_buf[lev+1]) // but no filter
        {
            MultiFab::Copy(*charge_buf[lev+1],
                           *rho_cp[lev+1], icomp, icomp, ncomp,
                           rho_cp[lev+1]->nGrow());
            mf.ParallelAdd(*charge_buf[lev+1], icomp, 0,
                           ncomp,
                           charge_buf[lev+1]->nGrowVect(), IntVect::TheZeroVector(),
                           period);
            WarpXSumGuardCells(*(rho_cp[lev+1]), period, icomp, ncomp);
        }
        else // no filter, no buffer
        {
            mf.ParallelAdd(*rho_cp[lev+1], icomp, 0, ncomp,
                           rho_cp[lev+1]->nGrowVect(), IntVect::TheZeroVector(),
                           period);
            WarpXSumGuardCells(*(rho_cp[lev+1]), period, icomp, ncomp);
        }
        MultiFab::Add(*rho_fp[lev], mf, 0, icomp, ncomp, 0);
        NodalSyncRho(lev+1, PatchType::coarse, icomp, ncomp);
    }

    NodalSyncRho(lev, PatchType::fine, icomp, ncomp);
}

void
WarpX::NodalSyncJ (int lev, PatchType patch_type)
{
    if (patch_type == PatchType::fine)
    {
        const auto& period = Geom(lev).periodicity();
        current_fp[lev][0]->OverrideSync(period);
        current_fp[lev][1]->OverrideSync(period);
        current_fp[lev][2]->OverrideSync(period);
    }
    else if (patch_type == PatchType::coarse)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        current_cp[lev][0]->OverrideSync(cperiod);
        current_cp[lev][1]->OverrideSync(cperiod);
        current_cp[lev][2]->OverrideSync(cperiod);
    }
}

void
WarpX::NodalSyncRho (int lev, PatchType patch_type, int icomp, int ncomp)
{
    if (patch_type == PatchType::fine && rho_fp[lev])
    {
        const auto& period = Geom(lev).periodicity();
        MultiFab rhof(*rho_fp[lev], amrex::make_alias, icomp, ncomp);
        rhof.OverrideSync(period);
    }
    else if (patch_type == PatchType::coarse && rho_cp[lev])
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        MultiFab rhoc(*rho_cp[lev], amrex::make_alias, icomp, ncomp);
        rhoc.OverrideSync(cperiod);
    }
}
