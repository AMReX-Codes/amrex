#include <AMReX_MLEBTensorOp.H>
#include <AMReX_MLTensor_K.H>
#include <AMReX_MLEBTensor_K.H>

namespace amrex {

void
MLEBTensorOp::applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                             BCMode bc_mode, StateMode, const MLMGBndry* bndry) const
{
    const int inhomog = bc_mode == BCMode::Inhomogeneous;
    const int imaxorder = maxorder;
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];
    const auto& maskvals = m_maskvals[amrlev][mglev];

    FArrayBox foofab(Box::TheUnitBox(),3);
    const auto& foo = foofab.array();

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vel, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();

        auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;
        if (fabtyp != FabType::covered)
        {
            const auto& velfab = vel.array(mfi);

            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            GpuArray<BoundCond,2*AMREX_SPACEDIM*AMREX_SPACEDIM> bct;
            GpuArray<Real,2*AMREX_SPACEDIM*AMREX_SPACEDIM> bcl;
            for (OrientationIter face; face; ++face) {
                Orientation ori = face();
                const int iface = ori;
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    bct[iface*AMREX_SPACEDIM+icomp] = bdcv[icomp][ori];
                    bcl[iface*AMREX_SPACEDIM+icomp] = bdlv[icomp][ori];
                }
            }

            const auto& mxlo = maskvals[Orientation(0,Orientation::low )].array(mfi);
            const auto& mylo = maskvals[Orientation(1,Orientation::low )].array(mfi);
            const auto& mxhi = maskvals[Orientation(0,Orientation::high)].array(mfi);
            const auto& myhi = maskvals[Orientation(1,Orientation::high)].array(mfi);

            const auto& bvxlo = (bndry != nullptr) ?
                (*bndry)[Orientation(0,Orientation::low )].array(mfi) : foo;
            const auto& bvylo = (bndry != nullptr) ?
                (*bndry)[Orientation(1,Orientation::low )].array(mfi) : foo;
            const auto& bvxhi = (bndry != nullptr) ?
                (*bndry)[Orientation(0,Orientation::high)].array(mfi) : foo;
            const auto& bvyhi = (bndry != nullptr) ?
                (*bndry)[Orientation(1,Orientation::high)].array(mfi) : foo;

#if (AMREX_SPACEDIM == 2)

            AMREX_HOST_DEVICE_FOR_1D ( 4, icorner,
            {
                mltensor_fill_corners(icorner, vbx, velfab,
                                      mxlo, mylo, mxhi, myhi,
                                      bvxlo, bvylo, bvxhi, bvyhi,
                                      bct, bcl, inhomog, imaxorder,
                                      dxinv, domain);
            });
#else
            const auto& mzlo = maskvals[Orientation(2,Orientation::low )].array(mfi);
            const auto& mzhi = maskvals[Orientation(2,Orientation::high)].array(mfi);

            const auto& bvzlo = (bndry != nullptr) ?
                (*bndry)[Orientation(2,Orientation::low )].array(mfi) : foo;
            const auto& bvzhi = (bndry != nullptr) ?
                (*bndry)[Orientation(2,Orientation::high)].array(mfi) : foo;

            AMREX_HOST_DEVICE_FOR_1D ( 12, iedge,
            {
                mltensor_fill_edges(iedge, vbx, velfab,
                                    mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                    bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                    bct, bcl, inhomog, imaxorder,
                                    dxinv, domain);
            });

            AMREX_HOST_DEVICE_FOR_1D ( 8, icorner,
            {
                mltensor_fill_corners(icorner, vbx, velfab,
                                      mxlo, mylo, mzlo, mxhi, myhi, mzhi,
                                      bvxlo, bvylo, bvzlo, bvxhi, bvyhi, bvzhi,
                                      bct, bcl, inhomog, imaxorder,
                                      dxinv, domain);
            });
#endif
        }
    }

    // Notet that it is incorrect to call EnforcePeriodicity on vel.
}

}
