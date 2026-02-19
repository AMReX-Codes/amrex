#include <AMReX_EdgeFluxRegister.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

#if (AMREX_SPACEDIM == 3)
static constexpr std::array<IntVect,3> E_ixtype{IntVect(0,1,1),IntVect(1,0,1),IntVect(1,1,0)};
#endif

EdgeFluxRegister::EdgeFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                    const DistributionMapping& fdm, const DistributionMapping& cdm,
                                    const Geometry& fgeom, const Geometry& cgeom,
                                    int nvar)
{
    define(fba, cba, fdm, cdm, fgeom, cgeom, nvar);
}

void EdgeFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
                               const DistributionMapping& fdm, const DistributionMapping& cdm,
                               const Geometry& fgeom, const Geometry& cgeom,
                               int nvar)
{
    m_fine_geom = fgeom;
    m_crse_geom = cgeom;
    m_ratio = fgeom.Domain().size() / cgeom.Domain().size();
    AMREX_ALWAYS_ASSERT(fgeom.Domain() == amrex::refine(cgeom.Domain(), m_ratio));
    m_ncomp = nvar;

    m_has_cf.define(cba, cdm);
    for (MFIter mfi(m_has_cf, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
        m_has_cf[mfi] = 0;
    }

#if (AMREX_SPACEDIM == 3)

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_E_crse[idim].define(amrex::convert(cba,E_ixtype[idim]), cdm, nvar, 0);
    }
    for (OrientationIter oit; oit.isValid(); ++oit) {
        auto face = oit();
        const int direction = face.coordDir();
        int count = 0;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim != direction) {
                BoxArray bba(amrex::convert(fba,IntVect(0)),
                             BATransformer(face, IndexType(E_ixtype[idim]), 0, 1, 0));
                bba.coarsen(m_ratio);
                m_E_fine[face][count++].define(bba, fdm, nvar, 0);
            }
        }
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Array<MultiFab const*,4> fmf;
        int count = 0;
        for (OrientationIter oit; oit.isValid(); ++oit) {
            auto face = oit();
            const int direction = face.coordDir();
            if (direction != idim) {
                // For x-direction, we store Ey and then Ez in m_E_fine.
                // For y-direction, we store Ex and then Ez in m_E_fine.
                // For z-direction, we store Ey and then Ez in m_E_fine.
                const int m = (idim < direction) ? idim : idim-1;
                fmf[count++] = & m_E_fine[face][m];
            }
        }

        for (int m = 0; m < 4; ++m) {
            LayoutData<int> tmp_has_cf;
            // We use IntVect(1) as ref ratio because fmf has already be coarsened
            auto tmp_mask = makeFineMask(m_E_crse[idim], *fmf[m], IntVect(0), IntVect(1),
                                         m_crse_geom.periodicity(), 0, 1, tmp_has_cf);
            if (m == 0) {
                m_fine_mask[idim] = std::move(tmp_mask);
            } else {
                amrex::Add(m_fine_mask[idim], tmp_mask, 0, 0, 1, 0);
            }
            for (MFIter mfi(m_has_cf, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
                m_has_cf[mfi] += tmp_has_cf[mfi];
            }
        }
    }

#else /* 2D */

    m_E_crse.define(amrex::convert(cba,IntVect(1)), cdm, nvar, 0);

    for (OrientationIter oit; oit.isValid(); ++oit) {
        auto face = oit();
        BoxArray bba(amrex::convert(fba,IntVect(0)),
                     BATransformer(face, IndexType::TheNodeType(), 0, 1, 0));
        bba.coarsen(m_ratio);
        m_E_fine[face].define(bba, fdm, nvar, 0);
    }

    for (OrientationIter oit; oit.isValid(); ++oit) {
        auto face = oit();
        LayoutData<int> tmp_has_cf;
        // We use IntVect(1) as ref ratio because fmf has already be coarsened
        auto tmp_mask = makeFineMask(m_E_crse, m_E_fine[face], IntVect(0), IntVect(1),
                                     m_crse_geom.periodicity(), 0, 1, tmp_has_cf);
        if (int(face) == 0) {
            m_fine_mask = std::move(tmp_mask);
        } else {
            amrex::Add(m_fine_mask, tmp_mask, 0, 0, 1, 0);
        }
        for (MFIter mfi(m_has_cf, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
            m_has_cf[mfi] += tmp_has_cf[mfi];
        }
    }

#endif
}

void EdgeFluxRegister::reset ()
{
#if (AMREX_SPACEDIM == 3)

    for (auto& mf : m_E_crse) {
        auto const& ma = mf.arrays();
        ParallelFor(mf, IntVect(0), mf.nComp(),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
        {
            ma[bno](i,j,k,n) = Real(0.0);
        });
    }
    for (auto& a : m_E_fine) {
        for (auto& mf : a) {
#ifdef AMREX_USE_GPU
            auto const& ma = mf.arrays();
            ParallelFor(mf, IntVect(0), mf.nComp(),
            [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
            {
                ma[bno](i,j,k,n) = Real(0.0);
            });
#else
            // Due to its special BoxArray, it's not safe do tiling
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                mf[mfi].template setVal<RunOn::Host>(Real(0.0));
            }
#endif
        }
    }

#else /* 2D */

    {
        auto const& ma = m_E_crse.arrays();
        ParallelFor(m_E_crse, IntVect(0), m_E_crse.nComp(),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
        {
            ma[bno](i,j,k,n) = Real(0.0);
        });
    }
    for (auto& mf : m_E_fine) {
#ifdef AMREX_USE_GPU
        auto const& ma = mf.arrays();
        ParallelFor(mf, IntVect(0), mf.nComp(),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
        {
            ma[bno](i,j,k,n) = Real(0.0);
        });
#else
        // Due to its special BoxArray, it's not safe do tiling
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            mf[mfi].template setVal<RunOn::Host>(Real(0.0));
        }
#endif
    }

#endif

    Gpu::synchronize();
}

#if (AMREX_SPACEDIM == 3)

void EdgeFluxRegister::CrseAdd (MFIter const& mfi, const Array<FArrayBox const*,3>& E_crse,
                                Real dt_crse)
{
    AMREX_ASSERT(mfi.validbox() == mfi.tilebox());
    if (m_has_cf[mfi]) {
        auto const& dst0 = m_E_crse[0].array(mfi);
        auto const& dst1 = m_E_crse[1].array(mfi);
        auto const& dst2 = m_E_crse[2].array(mfi);
        auto const& src0 = E_crse[0]->const_array();
        auto const& src1 = E_crse[1]->const_array();
        auto const& src2 = E_crse[2]->const_array();
        amrex::ParallelFor
            (amrex::convert(mfi.validbox(),E_ixtype[0]), m_ncomp,
             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
             {
                 dst0(i,j,k,n) += src0(i,j,k,n) * dt_crse;
             },
             amrex::convert(mfi.validbox(),E_ixtype[1]), m_ncomp,
             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
             {
                 dst1(i,j,k,n) += src1(i,j,k,n) * dt_crse;
             },
             amrex::convert(mfi.validbox(),E_ixtype[2]), m_ncomp,
             [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
             {
                 dst2(i,j,k,n) += src2(i,j,k,n) * dt_crse;
             });
    }
}

void EdgeFluxRegister::FineAdd (MFIter const& mfi, const Array<FArrayBox const*,3>& E_fine,
                                Real dt_fine)
{
    AMREX_ASSERT(mfi.validbox() == mfi.tilebox());
    auto const ratio = m_ratio;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        auto const& src = E_fine[idim]->const_array();
        for (OrientationIter oit; oit.isValid(); ++oit) {
            auto face = oit();
            const int direction = face.coordDir();
            if (direction != idim) {
                // For x-direction, we store Ey and then Ez in m_E_fine.
                // For y-direction, we store Ex and then Ez in m_E_fine.
                // For z-direction, we store Ey and then Ez in m_E_fine.
                const int m = (idim < direction) ? idim : idim-1;
                auto const& dst = m_E_fine[face][m].array(mfi);
                AMREX_ASSERT(E_fine[idim]->box().ixType() == m_E_fine[face][m].ixType());
                auto offset = IntVect::TheDimensionVector(idim).dim3();
                auto dt2 = dt_fine / Real(ratio[idim]);
                ParallelFor(Box(dst), m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    int ii = ratio[0]*i;
                    int jj = ratio[1]*j;
                    int kk = ratio[2]*k;
                    for (int rr = 0; rr < ratio[idim]; ++rr) {
                        dst(i,j,k,n) += src(ii+offset.x*rr,jj+offset.y*rr,kk+offset.z*rr,n)*dt2;
                    }
                });
            }
        }
    }
}

#else /* 2D */

void EdgeFluxRegister::CrseAdd (MFIter const& mfi, FArrayBox const& E_crse, Real dt_crse)
{
    AMREX_ASSERT(mfi.validbox() == mfi.tilebox());
    if (m_has_cf[mfi]) {
        auto const& dst = m_E_crse.array(mfi);
        auto const& src = E_crse.const_array();
        amrex::ParallelFor(amrex::convert(mfi.validbox(),IntVect(1)), m_ncomp,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                           {
                               dst(i,j,k,n) += src(i,j,k,n) * dt_crse;
                           });
    }
}

void EdgeFluxRegister::FineAdd (MFIter const& mfi, FArrayBox const& E_fine, Real dt_fine)
{
    AMREX_ASSERT(mfi.validbox() == mfi.tilebox());
    auto const ratio = m_ratio;
    auto const& src = E_fine.const_array();
    for (OrientationIter oit; oit.isValid(); ++oit) {
        auto face = oit();
        auto const& dst = m_E_fine[face].array(mfi);
        ParallelFor(Box(dst), m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int, int n)
        {
            int ii = ratio[0]*i;
            int jj = ratio[1]*j;
            dst(i,j,0,n) += src(ii,jj,0,n) * dt_fine;
        });
    }

}

#endif

void EdgeFluxRegister::Reflux (Array<MultiFab*,AMREX_SPACEDIM> const& B_crse) const
{
#if (AMREX_SPACEDIM == 3)

    Array<MultiFab,AMREX_SPACEDIM> E_cfine;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        E_cfine[idim].define(m_E_crse[idim].boxArray(), m_E_crse[idim].DistributionMap(), m_ncomp, 0);
        for (OrientationIter oit; oit.isValid(); ++oit) {
            auto face = oit();
            const int direction = face.coordDir();
            if (direction != idim) {
                // For x-direction, we store Ey and then Ez in m_E_fine.
                // For y-direction, we store Ex and then Ez in m_E_fine.
                // For z-direction, we store Ey and then Ez in m_E_fine.
                const int m = (idim < direction) ? idim : idim-1;
                E_cfine[idim].ParallelCopy(m_E_fine[face][m], m_crse_geom.periodicity());
            }
        }
    }

    Real dxi = m_crse_geom.InvCellSize(0);
    Real dyi = m_crse_geom.InvCellSize(1);
    Real dzi = m_crse_geom.InvCellSize(2);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*B_crse[0]); mfi.isValid(); ++mfi) {
        if (m_has_cf[mfi]) {
            Box xbx = amrex::convert(mfi.validbox(), B_crse[0]->ixType());
            Box ybx = amrex::convert(mfi.validbox(), B_crse[1]->ixType());
            Box zbx = amrex::convert(mfi.validbox(), B_crse[2]->ixType());
            auto const& xmsk = m_fine_mask[0].const_array(mfi);
            auto const& ymsk = m_fine_mask[1].const_array(mfi);
            auto const& zmsk = m_fine_mask[2].const_array(mfi);
            auto const& Bx = B_crse[0]->array(mfi);
            auto const& By = B_crse[1]->array(mfi);
            auto const& Bz = B_crse[2]->array(mfi);
            auto const& cEx = m_E_crse[0].const_array(mfi);
            auto const& cEy = m_E_crse[1].const_array(mfi);
            auto const& cEz = m_E_crse[2].const_array(mfi);
            auto const& fEx = E_cfine[0].const_array(mfi);
            auto const& fEy = E_cfine[1].const_array(mfi);
            auto const& fEz = E_cfine[2].const_array(mfi);
            ParallelFor
                (xbx, m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                 {
                     Real dEym = ymsk(i,j  ,k  ) ? (fEy(i,j  ,k  ,n)-cEy(i,j  ,k  ,n)) : Real(0.0);
                     Real dEyp = ymsk(i,j  ,k+1) ? (fEy(i,j  ,k+1,n)-cEy(i,j  ,k+1,n)) : Real(0.0);
                     Real dEzm = zmsk(i,j  ,k  ) ? (fEz(i,j  ,k  ,n)-cEz(i,j  ,k  ,n)) : Real(0.0);
                     Real dEzp = zmsk(i,j+1,k  ) ? (fEz(i,j+1,k  ,n)-cEz(i,j+1,k  ,n)) : Real(0.0);
                     Bx(i,j,k,n) -= (dEzp-dEzm)*dyi - (dEyp-dEym)*dzi;
                 },
                 ybx, m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                 {
                     Real dExm = xmsk(i  ,j,k  ) ? (fEx(i  ,j,k  ,n)-cEx(i  ,j,k  ,n)) : Real(0.0);
                     Real dExp = xmsk(i  ,j,k+1) ? (fEx(i  ,j,k+1,n)-cEx(i  ,j,k+1,n)) : Real(0.0);
                     Real dEzm = zmsk(i  ,j,k  ) ? (fEz(i  ,j,k  ,n)-cEz(i  ,j,k  ,n)) : Real(0.0);
                     Real dEzp = zmsk(i+1,j,k  ) ? (fEz(i+1,j,k  ,n)-cEz(i+1,j,k  ,n)) : Real(0.0);
                     By(i,j,k,n) -= (dExp-dExm)*dzi - (dEzp-dEzm)*dxi;
                 },
                 zbx, m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                 {
                     Real dExm = xmsk(i  ,j  ,k) ? (fEx(i  ,j  ,k,n)-cEx(i  ,j  ,k,n)) : Real(0.0);
                     Real dExp = xmsk(i  ,j+1,k) ? (fEx(i  ,j+1,k,n)-cEx(i  ,j+1,k,n)) : Real(0.0);
                     Real dEym = ymsk(i  ,j  ,k) ? (fEy(i  ,j  ,k,n)-cEy(i  ,j  ,k,n)) : Real(0.0);
                     Real dEyp = ymsk(i+1,j  ,k) ? (fEy(i+1,j  ,k,n)-cEy(i+1,j  ,k,n)) : Real(0.0);
                     Bz(i,j,k,n) -= (dEyp-dEym)*dxi - (dExp-dExm)*dyi;
                 });
        }
    }

#else /* 2D */

    MultiFab E_cfine(m_E_crse.boxArray(), m_E_crse.DistributionMap(), m_ncomp, 0);

    for (OrientationIter oit; oit.isValid(); ++oit) {
        auto face = oit();
        E_cfine.ParallelCopy(m_E_fine[face], m_crse_geom.periodicity());
    }

    Real dxi = m_crse_geom.InvCellSize(0);
    Real dyi = m_crse_geom.InvCellSize(1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*B_crse[0]); mfi.isValid(); ++mfi) {
        if (m_has_cf[mfi]) {
            Box xbx = amrex::convert(mfi.validbox(), B_crse[0]->ixType());
            Box ybx = amrex::convert(mfi.validbox(), B_crse[1]->ixType());
            auto const& zmsk = m_fine_mask.const_array(mfi);
            auto const& Bx = B_crse[0]->array(mfi);
            auto const& By = B_crse[1]->array(mfi);
            auto const& cEz = m_E_crse.const_array(mfi);
            auto const& fEz = E_cfine.const_array(mfi);
            ParallelFor
                (xbx, m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                 {
                     Real dEzm = zmsk(i,j  ,k  ) ? (fEz(i,j  ,k  ,n)-cEz(i,j  ,k  ,n)) : Real(0.0);
                     Real dEzp = zmsk(i,j+1,k  ) ? (fEz(i,j+1,k  ,n)-cEz(i,j+1,k  ,n)) : Real(0.0);
                     Bx(i,j,k,n) -= (dEzp-dEzm)*dyi;
                 },
                 ybx, m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                 {
                     Real dEzm = zmsk(i  ,j,k  ) ? (fEz(i  ,j,k  ,n)-cEz(i  ,j,k  ,n)) : Real(0.0);
                     Real dEzp = zmsk(i+1,j,k  ) ? (fEz(i+1,j,k  ,n)-cEz(i+1,j,k  ,n)) : Real(0.0);
                     By(i,j,k,n) -= - (dEzp-dEzm)*dxi;
                 });
        }
    }

#endif
}

}
