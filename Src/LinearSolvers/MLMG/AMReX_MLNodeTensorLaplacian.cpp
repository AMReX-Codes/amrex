#include <AMReX_MLNodeTensorLaplacian.H>

namespace amrex {

MLNodeTensorLaplacian::MLNodeTensorLaplacian (const Vector<Geometry>& a_geom,
                                              const Vector<BoxArray>& a_grids,
                                              const Vector<DistributionMapping>& a_dmap,
                                              const LPInfo& a_info)
{}

MLNodeTensorLaplacian::~MLNodeTensorLaplacian ()
{}

void
MLNodeTensorLaplacian::define (const Vector<Geometry>& a_geom,
                               const Vector<BoxArray>& a_grids,
                               const Vector<DistributionMapping>& a_dmap,
                               const LPInfo& a_info)
{}

void
MLNodeTensorLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
}

void
MLNodeTensorLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine,
                                      const MultiFab& crse) const
{

}

void
MLNodeTensorLaplacian::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                               const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    amrex::Abort("MLNodeTensorLaplacian::averageDownSolutionRHS: TODO");
}

void
MLNodeTensorLaplacian::reflux (int crse_amrlev,
                               MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                               MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
    amrex::Abort("MLNodeTensorLaplacian::reflux: TODO");
}

void
MLNodeTensorLaplacian::prepareForSolve ()
{}

void
MLNodeTensorLaplacian::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode s_mode,
                                bool skip_fillboundary) const
{
}

void
MLNodeTensorLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
}

void
MLNodeTensorLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
}

void
MLNodeTensorLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
}

void
MLNodeTensorLaplacian::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    amrex::Abort("MLNodeTensorLaplacian::fixUpResidualMask: TODO");
}

}
