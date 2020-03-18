
#include <AMReX_MultiMask.H>
#include <AMReX_BndryData.H>

namespace amrex {

MultiMask::MultiMask (const BoxArray& ba, const DistributionMapping& dm, int ncomp)
    : m_fa(ba, dm, ncomp, 0, MFInfo(), DefaultFabFactory<Mask>())
{ }

MultiMask::MultiMask (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
		      Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval)
{
    define(regba, dm, geom, face, in_rad, out_rad, extent_rad, ncomp, initval);
}

void
MultiMask::define (const BoxArray& ba, const DistributionMapping& dm, int ncomp)
{
    BL_ASSERT(m_fa.size() == 0);
    m_fa.define(ba,dm,ncomp,0,MFInfo(),DefaultFabFactory<Mask>());
}

void
MultiMask::define (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
		   Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval)
{
    BL_ASSERT(m_fa.size() == 0);

    BoxArray mskba(regba, BATransformer(face,IndexType::TheCellType(),in_rad,out_rad,extent_rad));
    m_fa.define(mskba, dm, ncomp, 0, MFInfo(), DefaultFabFactory<Mask>());

    if (initval)
    {
        const int bndrydata_outside_domain = BndryData::outside_domain;
        const int bndrydata_not_covered = BndryData::not_covered;
        const int bndrydata_covered = BndryData::covered;

        int ngrow = std::max(out_rad, extent_rad);
        Box domain = geom.Domain();
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (geom.isPeriodic(i)) {
                domain.grow(i, ngrow);
            }
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_fa); mfi.isValid(); ++mfi)
        {
            auto const& fab = m_fa.array(mfi);
            Box const fbx{fab};
            AMREX_HOST_DEVICE_FOR_3D(fbx, i, j, k,
            {
                if (domain.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    fab(i,j,k) = bndrydata_not_covered;
                } else {
                    fab(i,j,k) = bndrydata_outside_domain;
                }
            });
        }

        FabArray<Mask> regmf(regba, dm, 1, 0, MFInfo().SetAlloc(false));
        const FabArrayBase::CPC& cpc = m_fa.getCPC(IntVect::TheZeroVector(),
                                                   regmf,
                                                   IntVect::TheZeroVector(),
                                                   geom.periodicity());
        m_fa.setVal(bndrydata_covered, cpc, 0, 1);
    }
}

void 
MultiMask::Copy (MultiMask& dst, const MultiMask& src)
{
    BL_ASSERT(dst.nComp() == src.nComp());
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.DistributionMap() == src.DistributionMap());
    const int ncomp = dst.nComp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst.m_fa); mfi.isValid(); ++mfi) {
        auto const srcfab = src.m_fa.array(mfi);
        auto       dstfab = dst.m_fa.array(mfi);
        const Box& bx = dst.m_fa[mfi].box();
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dstfab(i,j,k,n) = srcfab(i,j,k,n);
        });
    }
}

}
