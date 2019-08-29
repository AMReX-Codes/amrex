
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBFluxRegister_C.H>
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_cfp_inside_mask); mfi.isValid(); ++mfi)
    {
        const Box& ifabbox = mfi.fabbox();
        auto const& ifab = m_cfp_inside_mask.array(mfi);
        AMREX_HOST_DEVICE_FOR_3D(ifabbox, i, j, k,
        {
            ifab(i,j,k) = 0;
        });

        const int li = mfi.LocalIndex();
        const auto& cfp_fabs = m_cfp_fab[li];
        for (const FArrayBox* cfp : cfp_fabs)
        {
            const Box& bx = amrex::grow(cfp->box(), 1);
            const Box& ibx = bx & ifabbox;
            AMREX_HOST_DEVICE_FOR_3D(ibx, i, j, k,
            {
                ifab(i,j,k) = 1; // cells just inside crse/fine boundary
            });
        }
    }
}

void
EBFluxRegister::CrseAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& flux,
                         const Real* dx, Real dt,
                         const FArrayBox& volfrac,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac,
                         RunOn runon)
{
    BL_ASSERT(m_crse_data.nComp() == flux[0]->nComp());

    if (m_crse_fab_flag[mfi.LocalIndex()] == crse_cell) {
        return;  // this coarse fab is not close to fine fabs.
    }

    Array4<Real> const& fab = m_crse_data.array(mfi);
    const Box& bx = mfi.tilebox();
    const int nc = fab.nComp();

    Array4<int const> const& amrflag = m_crse_flag.array(mfi);

    AMREX_D_TERM(Real dtdx = dt/dx[0];,
                 Real dtdy = dt/dx[1];,
                 Real dtdz = dt/dx[2];);
    AMREX_D_TERM(Array4<Real const> const& fx = flux[0]->const_array();,
                 Array4<Real const> const& fy = flux[1]->const_array();,
                 Array4<Real const> const& fz = flux[2]->const_array(););
    AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array();,
                 Array4<Real const> const& apy = areafrac[1]->const_array();,
                 Array4<Real const> const& apz = areafrac[2]->const_array(););
    Array4<Real const> const& vfrac = volfrac.const_array();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu, bx, nc, i, j, k, n,
    {
        eb_flux_reg_crseadd_va(i,j,k,n,fab,amrflag,AMREX_D_DECL(fx,fy,fz),
                               vfrac,AMREX_D_DECL(apx,apy,apz),
                               AMREX_D_DECL(dtdx,dtdy,dtdz));
    });
}


void
EBFluxRegister::FineAdd (const MFIter& mfi,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& a_flux,
                         const Real* dx, Real dt,
                         const FArrayBox& volfrac,
                         const std::array<FArrayBox const*, AMREX_SPACEDIM>& areafrac,
                         const FArrayBox& dm,
                         RunOn runon)
{
    BL_ASSERT(m_cfpatch.nComp() == a_flux[0]->nComp());

    const int li = mfi.LocalIndex();
    Vector<FArrayBox*>& cfp_fabs = m_cfp_fab[li];
    if (cfp_fabs.empty()) return;

    bool run_on_gpu = (runon == RunOn::Gpu and Gpu::inLaunchRegion());

    const int nc = m_cfpatch.nComp();

    const Box& tbx = mfi.tilebox();
    BL_ASSERT(tbx.cellCentered());
    const Box& cbx = amrex::coarsen(tbx, m_ratio);
    const Box& fbx = amrex::refine(cbx, m_ratio);

    AMREX_D_TERM(Array4<Real const> const& fx = a_flux[0]->const_array();,
                 Array4<Real const> const& fy = a_flux[1]->const_array();,
                 Array4<Real const> const& fz = a_flux[2]->const_array(););

    Array4<Real const> const& vfrac = volfrac.const_array();
    AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array();,
                 Array4<Real const> const& apy = areafrac[1]->const_array();,
                 Array4<Real const> const& apz = areafrac[2]->const_array(););

    Dim3 ratio = m_ratio.dim3();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        Real fac = dt/dx[idim];
        const Box& lobx = amrex::adjCellLo(cbx, idim);
        const Box& hibx = amrex::adjCellHi(cbx, idim);
        for (FArrayBox* cfp : cfp_fabs)
        {
            Array4<Real> const& cfa = cfp->array();
            const Box& lobx_is = lobx & cfp->box();
            if (lobx_is.ok()) {
                if (idim == 0)
                {
                    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu,lobx_is,nc,i,j,k,n,
                    {
                        eb_flux_reg_fineadd_va_xlo(i,j,k,n, cfa, fx, vfrac, apx, fac, ratio);
                    });
                }
#if (AMREX_SPACEDIM >= 2)
                else if (idim == 1)
                {
                    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu,lobx_is,nc,i,j,k,n,
                    {
                        eb_flux_reg_fineadd_va_ylo(i,j,k,n, cfa, fy, vfrac, apy, fac, ratio);
                    });
                }
#if (AMREX_SPACEDIM == 3)
                else
                {
                    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu,lobx_is,nc,i,j,k,n,
                    {
                        eb_flux_reg_fineadd_va_zlo(i,j,k,n, cfa, fz, vfrac, apz, fac, ratio);
                    });
                }
#endif
#endif
            }
            const Box& hibx_is = hibx & cfp->box();
            if (hibx_is.ok()) {
                if (idim == 0)
                {
                    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu,hibx_is,nc,i,j,k,n,
                    {
                        eb_flux_reg_fineadd_va_xhi(i,j,k,n, cfa, fx, vfrac, apx, fac, ratio);
                    });
                }
#if (AMREX_SPACEDIM >= 2)
                else if (idim == 1)
                {
                    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu,hibx_is,nc,i,j,k,n,
                    {
                        eb_flux_reg_fineadd_va_yhi(i,j,k,n, cfa, fy, vfrac, apy, fac, ratio);
                    });
                }
#if (AMREX_SPACEDIM == 3)
                else
                {
                    AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu,hibx_is,nc,i,j,k,n,
                    {
                        eb_flux_reg_fineadd_va_zhi(i,j,k,n, cfa, fz, vfrac, apz, fac, ratio);
                    });
                }
#endif
#endif
            }
        }
    }

    Real threshold = amrex_eb_get_reredistribution_threshold()*(AMREX_D_TERM(ratio.x,*ratio.y,*ratio.z));
    const Box& tbxg1 = amrex::grow(tbx,1);
    const Box& cbxg1 = amrex::grow(cbx,1);
    Array4<Real const> const& dma = dm.const_array();
    for (FArrayBox* cfp : cfp_fabs)
    {
        const Box& wbx = cbxg1 & cfp->box();
        if (wbx.ok())
        {
            Array4<Real> const& cfa = cfp->array();
            AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu, wbx, nc, i, j, k, n,
            {
                eb_flux_reg_fineadd_dm(i,j,k,n,tbxg1, cfa, dma, vfrac, ratio, threshold);
            });
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
        
        auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(crse_state.Factory());
        auto const& flags = factory.getMultiEBCellFlagFab();

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
                    
                    const auto& ebflag = flags[mfi];
                    
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

    auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(fine_state.Factory());
    auto const& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cf,true); mfi.isValid(); ++mfi)
    {
        const Box& cbx = mfi.tilebox();
        const Box& fbx = amrex::refine(cbx, m_ratio);
        
        const auto& ebflag = flags[mfi];
        
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
