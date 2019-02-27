
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLap_F.H>

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
    const int owner = 1;
    const int nonowner = 0;

    const BoxArray& ba = amrex::convert(a_ba, IntVect::TheNodeVector());
    std::unique_ptr<iMultiFab> p{new iMultiFab(ba,dm,1,0, MFInfo(),
                                               DefaultFabFactory<IArrayBox>())};
    p->setVal(owner);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();
        
        for (MFIter mfi(*p); mfi.isValid(); ++mfi)
        {
            IArrayBox& fab = (*p)[mfi];
            const Box& bx = fab.box();
            const int i = mfi.index();
            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);                    
                for (const auto& is : isects)
                {
                    const int idx_other = is.first;
                    if (idx_other == i && iv == IntVect::TheZeroVector()) {
                        continue;
                    }
                    const Box& bx_other = ba[idx_other];

                    const Box& sect_other = is.second;
                    const Box& sect_this  = sect_other - iv;

                    const IntVect& sect_corner_other = sect_other.smallEnd();
                    const IntVect& sect_corner_this  = sect_this.smallEnd();

                    const long dist_other = bx_other.index(sect_corner_other);
                    const long dist_this  = bx.index(sect_corner_this);

                    bool yield = dist_this < dist_other;
                    if (dist_this == dist_other) {  // gauranteed to be two different boxes
                        yield = idx_other < i;
                    }

                    if (yield) {
                        fab.setVal(nonowner, sect_this, 0, 1);
                    }
                }
            }
        }
    }

    return p;
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
    // todo: gpu
    const int mglev = 0;
    const int ncomp = b.nComp();
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][0];
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(resid, true); mfi.isValid(); ++mfi)
    {
	    const Box& bx = mfi.tilebox();
	    amrex_mlndlap_solution_residual(BL_TO_FORTRAN_BOX(bx),
					    BL_TO_FORTRAN_ANYD(resid[mfi]),
					    BL_TO_FORTRAN_ANYD(b[mfi]),
					    BL_TO_FORTRAN_ANYD(dmsk[mfi]),
					    &ncomp);
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

}

