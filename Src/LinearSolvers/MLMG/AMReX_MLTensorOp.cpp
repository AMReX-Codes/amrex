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

    MLABecLaplacian::setScalars(1.0,1.0);

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
MLTensorOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                   StateMode s_mode, const MLMGBndry* bndry) const
{
    BL_PROFILE("MLTensorOp::apply()");

    MLABecLaplacian::apply(amrlev, mglev, out, in, bc_mode, s_mode, bndry);

    // Note that we cannot have inhomog bc_mode at mglev>0.
    if (mglev == 0 && bc_mode == BCMode::Inhomogeneous && s_mode == StateMode::Solution)
    {
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const axfab = out.array(mfi);
            Array4<Real const> const xfab = in.array(mfi);
            Array4<Real const> const gradetafab = m_gradeta[amrlev].array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mltensor_adotx_add_extra(bx, axfab, xfab, gradetafab, dxinv);
            });
        }
    }
}

}
