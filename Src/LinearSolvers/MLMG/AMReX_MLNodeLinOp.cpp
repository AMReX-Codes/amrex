
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLap_F.H>
#include <AMReX_MLNodeLap_K.H>

#ifdef _OPENMP
#include <omp.h>
#endif

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
                     const LPInfo& a_info,
                     const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    m_owner_mask.resize(m_num_amr_levels);
    m_dirichlet_mask.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_owner_mask[amrlev].resize(m_num_mg_levels[amrlev]);
        m_dirichlet_mask[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_owner_mask[amrlev][mglev] = makeOwnerMask(m_grids[amrlev][mglev],
                                                        m_dmap[amrlev][mglev],
                                                        m_geom[amrlev][mglev]);
            m_dirichlet_mask[amrlev][mglev].reset
                (new iMultiFab(amrex::convert(m_grids[amrlev][mglev],IntVect::TheNodeVector()),
                               m_dmap[amrlev][mglev], 1, 0));
        }
    }


    m_cc_fine_mask.resize(m_num_amr_levels);
    m_nd_fine_mask.resize(m_num_amr_levels);
    m_has_fine_bndry.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (amrlev < m_num_amr_levels-1)
        {
            m_cc_fine_mask[amrlev].reset(new iMultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1));
            m_nd_fine_mask[amrlev].reset(new iMultiFab(amrex::convert(m_grids[amrlev][0],IntVect::TheNodeVector()),
                                                       m_dmap[amrlev][0], 1, 0));
        }
        m_has_fine_bndry[amrlev].reset(new LayoutData<int>(m_grids[amrlev][0], m_dmap[amrlev][0]));
    }
}

std::unique_ptr<iMultiFab>
MLNodeLinOp::makeOwnerMask (const BoxArray& a_ba, const DistributionMapping& dm,
                            const Geometry& geom)
{
    const BoxArray& ba = amrex::convert(a_ba, IntVect::TheNodeVector());
    MultiFab foo(ba,dm,1,0, MFInfo().SetAlloc(false));
    return foo.OwnerMask(geom.periodicity());
}

void
MLNodeLinOp::nodalSync (int amrlev, int mglev, MultiFab& mf) const
{
    mf.OverrideSync(*m_owner_mask[amrlev][mglev], m_geom[amrlev][mglev].periodicity());
}

void
MLNodeLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                               const MultiFab* crse_bcdata)
{
    const int mglev = 0;
    const int ncomp = b.nComp();
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][0];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resid, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& res = resid.array(mfi);
        Array4<Real const> const& bb = b.const_array(mfi);
        Array4<int const> const& dd = dmsk.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            if (dd(i,j,k)) {
                res(i,j,k,n) = 0.0;
            } else {
                res(i,j,k,n) = bb(i,j,k,n) - res(i,j,k,n);
            }
        });
    }
}

void
MLNodeLinOp::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                 BCMode bc_mode, const MultiFab* crse_bcdata)
{
    apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
    int ncomp = b.nComp();
    MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 0);
}

void
MLNodeLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                    StateMode s_mode, const MLMGBndry*) const
{
    applyBC(amrlev, mglev, in, bc_mode, s_mode);
    Fapply(amrlev, mglev, out, in);
}

void
MLNodeLinOp::smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                     bool skip_fillboundary) const
{
    if (!skip_fillboundary) {
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Solution);
    }
    Fsmooth(amrlev, mglev, sol, rhs);
}

Real
MLNodeLinOp::xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const
{
    AMREX_ASSERT(amrlev==0);
    AMREX_ASSERT(mglev+1==m_num_mg_levels[0] || mglev==0);
    const auto& mask = (mglev+1 == m_num_mg_levels[0]) ? m_bottom_dot_mask : m_coarse_dot_mask;
    const int ncomp = y.nComp();
    const int nghost = 0;
    MultiFab tmp(x.boxArray(), x.DistributionMap(), ncomp, 0);
    MultiFab::Copy(tmp, x, 0, 0, ncomp, nghost);
    for (int i = 0; i < ncomp; i++) {
        MultiFab::Multiply(tmp, mask, 0, i, 1, nghost);
    }
    Real result = MultiFab::Dot(tmp,0,y,0,ncomp,nghost,true);
    if (!local) {
        ParallelAllReduce::Sum(result, Communicator(amrlev, mglev));
    }
    return result;
}

void
MLNodeLinOp::applyInhomogNeumannTerm (int amrlev, MultiFab& rhs) const
{
    int ncomp = rhs.nComp();
    for (int n = 0; n < ncomp; ++n)
    {
        auto itlo = std::find(m_lo_inhomog_neumann[n].begin(),
                              m_lo_inhomog_neumann[n].end(),   1);
        auto ithi = std::find(m_hi_inhomog_neumann[n].begin(),
                              m_hi_inhomog_neumann[n].end(),   1);
        if (itlo != m_lo_inhomog_neumann[n].end() or
            ithi != m_hi_inhomog_neumann[n].end())
        {
            amrex::Abort("Inhomogeneous Neumann not supported for nodal solver");
        }
    }
}

}

