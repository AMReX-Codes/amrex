
#include <AMReX_MLNodeLinOp.H>

namespace amrex {

MLNodeLinOp::MLNodeLinOp ()
{
    m_ixtype = IntVect::TheNodeVector();
}

MLNodeLinOp::~MLNodeLinOp () {}

void
MLNodeLinOp::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info)
{
    MLLinOp::define(a_geom, a_grids, a_dmap, a_info);
}

void
MLNodeLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                               const MultiFab* crse_bcdata)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(crse_bcdata == nullptr, "solutionResidual not fully implemented");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrlev == 0, "solutionResidual not fully implemented");

    const int mglev = 0;
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous);
    MultiFab::Xpay(resid, -1.0, b, 0, 0, 1, 0);
}

void
MLNodeLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                    const MLMGBndry*) const
{
    in.FillBoundary(m_geom[amrlev][mglev].periodicity());
    Fapply(amrlev, mglev, out, in);
}

}

