#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLTensor_K.H>

namespace amrex {

MLTensorOp::MLTensorOp (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLTensorOp::~MLTensorOp ()
{}

void
MLTensorOp::define (const Vector<Geometry>& a_geom,
                    const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap,
                    const LPInfo& a_info,
                    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLTensorOp::define()");

    MLABecLaplacian::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    setScalars(1.0,1.0);

    m_gradeta.clear();
    m_gradeta.resize(NAMRLevels());
    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
        m_gradeta[amrlev].define(m_grids[amrlev][0], m_dmap[amrlev][0], AMREX_SPACEDIM, 0,
                                    MFInfo(), *m_factory[amrlev][0]);
    }
}

void
MLTensorOp::prepareForSolve ()
{
    MLABecLaplacian::prepareForSolve();
    for (int amrlev = 0; amrlev < NAMRLevels(); ++amrlev) {
        amrex::computeGradient(m_gradeta[amrlev], getBCoeffs(amrlev,0), m_geom[amrlev][0]);
    }
}

void
MLTensorOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                           const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLTensorOp::solutionResidual()");

    MLABecLaplacian::solutionResidual(amrlev, resid, x, b, crse_bcdata);

    const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resid, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const rfab = resid.array(mfi);
        Array4<Real const> const vfab = x.array(mfi);
        Array4<Real const> const gradetafab = m_gradeta[amrlev].array(mfi);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mltensor_resid_add_extra(bx, rfab, vfab, gradetafab, dxinv);
        });
    }
}

}
