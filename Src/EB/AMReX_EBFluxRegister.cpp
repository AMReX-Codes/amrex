
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBFluxRegister_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

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

    const IArrayBox& flag = m_crse_flag[mfi];

    amrex_eb_flux_reg_crseadd_va(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(fab),
                                 BL_TO_FORTRAN_ANYD(flag),
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
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt,
                         const FArrayBox& volfrac,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac,
                         const FArrayBox& dm)
{
    BL_ASSERT(m_cfpatch.nComp() == flux[0]->nComp());

    const int li = mfi.LocalIndex();
    Array<FArrayBox*>& fabs = m_cfp_fab[li];
    if (fabs.empty()) return;

    const Box& tbx = mfi.tilebox();

    BL_ASSERT(tbx.cellCentered());
    AMREX_ALWAYS_ASSERT(tbx.coarsenable(m_ratio));
    const Box& cbx = amrex::coarsen(tbx, m_ratio);
    const int nc = m_cfpatch.nComp();

    FArrayBox cvol;

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
    {
        const Box& lobx = amrex::adjCellLo(cbx, idim);
        const Box& hibx = amrex::adjCellHi(cbx, idim);
        for (FArrayBox* cfp : m_cfp_fab[li])
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

    for (FArrayBox* cfp : m_cfp_fab[li])
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

}
