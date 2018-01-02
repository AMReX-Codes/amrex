
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBFluxRegister_F.H>
#include <AMReX_EBFArrayBox.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

EBFluxRegister::EBFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                const DistributionMapping& fdm, const DistributionMapping& cdm,
                                const Geometry& fgeom, const Geometry& cgeom,
                                const IntVect& ref_ratio, int fine_lev, int nvar)
    : YAFluxRegister(fba,cba,fdm,cdm,fgeom,cgeom,ref_ratio,fine_lev,nvar)
{
    defineExtra(fba, fdm);
}

void
EBFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
                        const DistributionMapping& fdm, const DistributionMapping& cdm,
                        const Geometry& fgeom, const Geometry& cgeom,
                        const IntVect& ref_ratio, int fine_lev, int nvar)
{
    YAFluxRegister::define(fba,cba,fdm,cdm,fgeom,cgeom,ref_ratio,fine_lev,nvar);
    defineExtra(fba, fdm);
}


void
EBFluxRegister::defineExtra (const BoxArray& fba, const DistributionMapping& fdm)
{
    BoxArray cfba = fba;
    cfba.coarsen(m_ratio);
    m_cfp_inside_mask.define(cfba, fdm, 1, 0, MFInfo(),DefaultFabFactory<IArrayBox>());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(m_cfp_inside_mask); mfi.isValid(); ++mfi)
    {
        IArrayBox& ifab = m_cfp_inside_mask[mfi];
        ifab.setVal(0);

        const int li = mfi.LocalIndex();
        const auto& cfp_fabs = m_cfp_fab[li];
        for (const FArrayBox* cfp : cfp_fabs)
        {
            const Box& bx = amrex::grow(cfp->box(), 1);
            const Box& ibx = bx & ifab.box();
            ifab.setVal(1,ibx);  // cells just inside crse/fine boundary
        }
    }
}

void
EBFluxRegister::CrseAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt,
                         const FArrayBox& volfrac,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac)
{
    BL_ASSERT(m_crse_data.nComp() == flux[0]->nComp());

    if (m_crse_fab_flag[mfi.LocalIndex()] == crse_cell) {
        return;  // this coarse fab is not close to fine fabs.
    }

    FArrayBox& fab = m_crse_data[mfi];
    const Box& bx = mfi.tilebox();
    const int nc = fab.nComp();

    const IArrayBox& amrflag = m_crse_flag[mfi];

    amrex_eb_flux_reg_crseadd_va(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(fab),
                                 BL_TO_FORTRAN_ANYD(amrflag),
                                 AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*flux[0]),
                                              BL_TO_FORTRAN_ANYD(*flux[1]),
                                              BL_TO_FORTRAN_ANYD(*flux[2])),
                                 BL_TO_FORTRAN_ANYD(volfrac),
                                 AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*areafrac[0]),
                                              BL_TO_FORTRAN_ANYD(*areafrac[1]),
                                              BL_TO_FORTRAN_ANYD(*areafrac[2])),
                                 dx, &dt,&nc);
}


void
EBFluxRegister::FineAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& a_flux,
                         const Real* dx, Real dt,
                         const FArrayBox& volfrac,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac,
                         const FArrayBox& dm)
{
    BL_ASSERT(m_cfpatch.nComp() == a_flux[0]->nComp());

    const int li = mfi.LocalIndex();
    Vector<FArrayBox*>& cfp_fabs = m_cfp_fab[li];
    if (cfp_fabs.empty()) return;

    const int nc = m_cfpatch.nComp();

    const Box& tbx = mfi.tilebox();
    BL_ASSERT(tbx.cellCentered());
    const Box& cbx = amrex::coarsen(tbx, m_ratio);
    const Box& fbx = amrex::refine(cbx, m_ratio);
    
    std::array<FArrayBox const*,AMREX_SPACEDIM> flux{AMREX_D_DECL(a_flux[0],a_flux[1],a_flux[2])};
    std::array<FArrayBox,AMREX_SPACEDIM> ftmp;
    if (fbx != tbx) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Box& b = amrex::surroundingNodes(fbx,idim);
            ftmp[idim].resize(b,nc);
            ftmp[idim].setVal(0.0);
            ftmp[idim].copy(*a_flux[idim]);
            flux[idim] = &ftmp[idim];
        }
    }

    FArrayBox cvol;

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
    {
        const Box& lobx = amrex::adjCellLo(cbx, idim);
        const Box& hibx = amrex::adjCellHi(cbx, idim);
        for (FArrayBox* cfp : cfp_fabs)
        {
            {
                const Box& lobx_is = lobx & cfp->box();
                const int side = 0;
                if (lobx_is.ok())
                {
                    cvol.resize(lobx_is);
                    amrex_eb_flux_reg_fineadd_va(BL_TO_FORTRAN_BOX(lobx_is),
                                                 BL_TO_FORTRAN_ANYD(*cfp),
                                                 BL_TO_FORTRAN_ANYD(*flux[idim]),
                                                 BL_TO_FORTRAN_ANYD(cvol),
                                                 BL_TO_FORTRAN_ANYD(volfrac),
                                                 AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*areafrac[0]),
                                                              BL_TO_FORTRAN_ANYD(*areafrac[1]),
                                                              BL_TO_FORTRAN_ANYD(*areafrac[2])),
                                                 dx, &dt, &nc, &idim, &side,
                                                 m_ratio.getVect());
                }
            }
            {
                const Box& hibx_is = hibx & cfp->box();
                const int side = 1;
                if (hibx_is.ok())
                {
                    cvol.resize(hibx_is);
                    amrex_eb_flux_reg_fineadd_va(BL_TO_FORTRAN_BOX(hibx_is),
                                                 BL_TO_FORTRAN_ANYD(*cfp),
                                                 BL_TO_FORTRAN_ANYD(*flux[idim]),
                                                 BL_TO_FORTRAN_ANYD(cvol),
                                                 BL_TO_FORTRAN_ANYD(volfrac),
                                                 AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*areafrac[0]),
                                                              BL_TO_FORTRAN_ANYD(*areafrac[1]),
                                                              BL_TO_FORTRAN_ANYD(*areafrac[2])),
                                                 dx, &dt, &nc, &idim, &side,
                                                 m_ratio.getVect());
                }
            }
        }
    }

    FArrayBox dmgrow(amrex::grow(tbx,m_ratio),nc);
    dmgrow.setVal(0.0);
    const Box& tbxg1 = amrex::grow(tbx,1);
    dmgrow.copy(dm,tbxg1,0,tbxg1,0,nc);

    const Box& cbxg1 = amrex::grow(cbx,1);

    for (FArrayBox* cfp : cfp_fabs)
    {
        const Box& wbx = cbxg1 & cfp->box();
        if (wbx.ok())
        {
            cvol.resize(wbx);
            amrex_eb_flux_reg_fineadd_dm(BL_TO_FORTRAN_BOX(wbx),
                                         BL_TO_FORTRAN_ANYD(*cfp),
                                         BL_TO_FORTRAN_ANYD(dmgrow),
                                         BL_TO_FORTRAN_ANYD(cvol),
                                         BL_TO_FORTRAN_ANYD(volfrac),
                                         dx, &nc, m_ratio.getVect());
        }
    }
}


void
EBFluxRegister::Reflux (MultiFab& crse_state, const amrex::MultiFab& crse_vfrac,
                        MultiFab& fine_state, const amrex::MultiFab& fine_vfrac)
{
    if (!m_cfp_mask.empty())
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(m_cfpatch); mfi.isValid(); ++mfi)
        {
            for (int i = 0; i < m_ncomp; ++i) {
                m_cfpatch[mfi].mult(m_cfp_mask[mfi],0,i);
            }
        }
    }

    m_crse_data.ParallelCopy(m_cfpatch, m_crse_geom.periodicity(), FabArrayBase::ADD);

    {
        MultiFab grown_crse_data(m_crse_data.boxArray(), m_crse_data.DistributionMap(),
                                 m_ncomp, 1, MFInfo(), FArrayBoxFactory());
        MultiFab::Copy(grown_crse_data, m_crse_data, 0, 0, m_ncomp, 0);
        grown_crse_data.FillBoundary(m_crse_geom.periodicity());
        
        m_crse_data.setVal(0.0);
        
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox fab;
            for (MFIter mfi(m_crse_data, MFItInfo().EnableTiling().SetDynamic(true));
                 mfi.isValid(); ++mfi)
            {
                if (m_crse_fab_flag[mfi.LocalIndex()] == fine_cell) // fab has crse/fine cells
                {
                    const Box& bx = mfi.tilebox();
                    
                    auto& sfab = dynamic_cast<EBFArrayBox&>(crse_state[mfi]);
                    const auto& ebflag = sfab.getEBCellFlagFab();
                    
                    if (ebflag.getType(bx) != FabType::covered) {
                        if (ebflag.getType(amrex::grow(bx,1)) == FabType::regular)
                        {
                            // no re-reflux or re-re-redistribution
                            m_crse_data[mfi].plus(grown_crse_data[mfi],bx,0,0,m_ncomp);
                        }
                        else
                        {
                            fab.resize(amrex::grow(bx,2),m_ncomp);
                            fab.setVal(0.0);
                            amrex_eb_rereflux_from_crse(BL_TO_FORTRAN_BOX(bx),
                                                        BL_TO_FORTRAN_ANYD(fab),
                                                        BL_TO_FORTRAN_ANYD(grown_crse_data[mfi]),
                                                        BL_TO_FORTRAN_ANYD(m_crse_flag[mfi]),
                                                        BL_TO_FORTRAN_ANYD(ebflag),
                                                        BL_TO_FORTRAN_ANYD(crse_vfrac[mfi]),
                                                        &m_ncomp);
                            m_crse_data[mfi].plus(fab,bx,0,0,m_ncomp);
                        }
                    }
                }
            }
        }
    }

    MultiFab::Add(crse_state, m_crse_data, 0, 0, m_ncomp, 0);

    // The fine-covered cells of m_crse_data contain the data that should go to the fine level
    BoxArray ba = fine_state.boxArray();
    ba.coarsen(m_ratio);
    MultiFab cf(ba, fine_state.DistributionMap(), m_ncomp, 0, MFInfo(), FArrayBoxFactory());
    cf.ParallelCopy(m_crse_data);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cf,true); mfi.isValid(); ++mfi)
    {
        const Box& cbx = mfi.tilebox();
        const Box& fbx = amrex::refine(cbx, m_ratio);
        
        auto& sfab = dynamic_cast<EBFArrayBox&>(fine_state[mfi]);
        const auto& ebflag = sfab.getEBCellFlagFab();
        
        if (ebflag.getType(fbx) != FabType::covered)
        {
            amrex_eb_rereflux_to_fine(BL_TO_FORTRAN_BOX(cbx),
                                      BL_TO_FORTRAN_ANYD(fine_state[mfi]),
                                      BL_TO_FORTRAN_ANYD(cf[mfi]),
                                      BL_TO_FORTRAN_ANYD(m_cfp_inside_mask[mfi]),
                                      &m_ncomp, m_ratio.getVect());
        }
    }
}


}
