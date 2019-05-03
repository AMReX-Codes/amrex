#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFabUtil.H>

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
}

void
MLTensorOp::prepareForSolve ()
{
    MLABecLaplacian::prepareForSolve();
}

}
