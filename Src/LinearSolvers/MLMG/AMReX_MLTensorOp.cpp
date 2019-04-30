#include <AMReX_MLTensorOp.H>

namespace amrex {

MLTensorOp::MLTensorOp (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FArrayBox> const*>& a_factory)
{

}

MLTensorOp::~MLTensorOp ()
{
}

void
MLTensorOp::define (const Vector<Geometry>& a_geom,
                    const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap,
                    const LPInfo& a_info,
                    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
}

void
MLTensorOp::setACoeffs (int amrlev, const MultiFab& alpha)
{
}

void
MLTensorOp::setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
}

void
MLTensorOp::prepareForSolve ()
{
}

void
MLTensorOp::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
}

void
MLTensorOp::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
}

void
MLTensorOp::FFlux (int amrlev, const MFIter& mfi,
                   const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                   const FArrayBox& sol, Location /* loc */,
                   const int face_only) const
{

}

void
MLTensorOp::normalize (int amrlev, int mglev, MultiFab& mf) const
{
}

void
MLTensorOp::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                           Vector<Array<MultiFab,AMREX_SPACEDIM> >& b)
{
}

void
MLTensorOp::averageDownCoeffs ()
{
}

void
MLTensorOp::averageDownCoeffsToCoarseAmrLevel (int flev)
{
}

}
