
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLLinOp_K.H>
#include <AMReX_MLMG_K.H>
#include <AMReX_MultiFabUtil.H>

#ifndef BL_NO_FORT
#include <AMReX_MLLinOp_F.H>
#endif

namespace amrex {

#ifdef AMREX_SOFT_PERF_COUNTERS
// perf_counters
MLCellLinOp::Counters MLCellLinOp::perf_counters;
#endif

namespace {
    // Have to put it here due to CUDA extended lambda limitation
    struct ABCTag {
        Array4<Real> fab;
        Array4<Real const> bcval_lo;
        Array4<Real const> bcval_hi;
        Array4<int const> mask_lo;
        Array4<int const> mask_hi;
        Real bcloc_lo;
        Real bcloc_hi;
        Box bx;
        BoundCond bctype_lo;
        BoundCond bctype_hi;
        int blen;
        int comp;
        int dir;

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box const& box() const noexcept { return bx; }
    };

    struct PSTag {
        Array4<Real> flo;
        Array4<Real> fhi;
        Array4<int const> mlo;
        Array4<int const> mhi;
        Real bcllo;
        Real bclhi;
        Box bx;
        BoundCond bctlo;
        BoundCond bcthi;
        int blen;
        int comp;
        int dir;

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box const& box() const noexcept { return bx; }
    };

#ifdef AMREX_USE_EB
    struct PSEBTag {
        Array4<Real> flo;
        Array4<Real> fhi;
        Array4<Real const> ap;
        Array4<int const> mlo;
        Array4<int const> mhi;
        Real bcllo;
        Real bclhi;
        Box bx;
        BoundCond bctlo;
        BoundCond bcthi;
        int blen;
        int comp;
        int dir;

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box const& box() const noexcept { return bx; }
    };
#endif
}

MLCellLinOp::MLCellLinOp ()
{
    m_ixtype = IntVect::TheCellVector();
}

MLCellLinOp::~MLCellLinOp () {}

void
MLCellLinOp::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
    defineAuxData();
    defineBC();
}

void
MLCellLinOp::defineAuxData ()
{
    BL_PROFILE("MLCellLinOp::defineAuxData()");

    m_undrrelxr.resize(m_num_amr_levels);
    m_maskvals.resize(m_num_amr_levels);
    m_fluxreg.resize(m_num_amr_levels-1);
    m_norm_fine_mask.resize(m_num_amr_levels-1);

    const int ncomp = getNComp();

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_undrrelxr[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_undrrelxr[amrlev][mglev].define(m_grids[amrlev][mglev],
                                              m_dmap[amrlev][mglev],
                                              1, 0, 0, ncomp);
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_maskvals[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            for (OrientationIter oitr; oitr; ++oitr)
            {
                const Orientation face = oitr();
                const int ngrow = 1;
                const int extent = isCrossStencil() ? 0 : 1; // extend to corners
                m_maskvals[amrlev][mglev][face].define(m_grids[amrlev][mglev],
                                                       m_dmap[amrlev][mglev],
                                                       m_geom[amrlev][mglev],
                                                       face, 0, ngrow, extent, 1, true);
            }
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
    {
        const IntVect ratio{m_amr_ref_ratio[amrlev]};
        m_fluxreg[amrlev].define(m_grids[amrlev+1][0], m_grids[amrlev][0],
                                 m_dmap[amrlev+1][0], m_dmap[amrlev][0],
                                 m_geom[amrlev+1][0], m_geom[amrlev][0],
                                 ratio, amrlev+1, ncomp);
        m_norm_fine_mask[amrlev] = std::make_unique<iMultiFab>
            (makeFineMask(m_grids[amrlev][0], m_dmap[amrlev][0], m_grids[amrlev+1][0],
                          ratio, 1, 0));
    }

#if (AMREX_SPACEDIM != 3)
    m_has_metric_term = !m_geom[0][0].IsCartesian() && info.has_metric_term;
#endif
}

void
MLCellLinOp::defineBC ()
{
    BL_PROFILE("MLCellLinOp::defineBC()");

    const int ncomp = getNComp();

    m_bndry_sol.resize(m_num_amr_levels);
    m_crse_sol_br.resize(m_num_amr_levels);

    m_bndry_cor.resize(m_num_amr_levels);
    m_crse_cor_br.resize(m_num_amr_levels);

    m_robin_bcval.resize(m_num_amr_levels);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bndry_sol[amrlev] = std::make_unique<MLMGBndry>(m_grids[amrlev][0], m_dmap[amrlev][0],
                                                          ncomp, m_geom[amrlev][0]);
    }

    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = m_amr_ref_ratio[amrlev-1];
        BoxArray cba = m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        m_crse_sol_br[amrlev] = std::make_unique<BndryRegister>
            (cba, m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
    }

    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = m_amr_ref_ratio[amrlev-1];
        BoxArray cba = m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        m_crse_cor_br[amrlev] = std::make_unique<BndryRegister>
            (cba, m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
        m_crse_cor_br[amrlev]->setVal(0.0);
    }

    // This has be to done after m_crse_cor_br is defined.
    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bndry_cor[amrlev] = std::make_unique<MLMGBndry>(m_grids[amrlev][0], m_dmap[amrlev][0],
                                                          ncomp, m_geom[amrlev][0]);
        MultiFab bc_data(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, 1);
        bc_data.setVal(0.0);

        m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, ncomp,
                                            m_amr_ref_ratio[amrlev-1], BCRec());

        Vector<Array<LinOpBCType,AMREX_SPACEDIM> > bclohi
            (ncomp,Array<LinOpBCType,AMREX_SPACEDIM>{{AMREX_D_DECL(BCType::Dirichlet,
                                                                   BCType::Dirichlet,
                                                                   BCType::Dirichlet)}});
        m_bndry_cor[amrlev]->setLOBndryConds(bclohi, bclohi, m_amr_ref_ratio[amrlev-1], RealVect{});
    }

    m_bcondloc.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bcondloc[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_bcondloc[amrlev][mglev] = std::make_unique<BndryCondLoc>(m_grids[amrlev][mglev],
                                                                       m_dmap[amrlev][mglev],
                                                                       ncomp);
        }
    }
}

void
MLCellLinOp::setLevelBC (int amrlev, const MultiFab* a_levelbcdata, const MultiFab* robinbc_a,
                         const MultiFab* robinbc_b, const MultiFab* robinbc_f)
{
    BL_PROFILE("MLCellLinOp::setLevelBC()");

    AMREX_ALWAYS_ASSERT(amrlev >= 0 && amrlev < m_num_amr_levels);

    const int ncomp = getNComp();

    MultiFab zero;
    IntVect ng(1);
    if (hasHiddenDimension()) ng[hiddenDirection()] = 0;
    if (a_levelbcdata == nullptr) {
        zero.define(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp, ng);
        zero.setVal(0.0);
    } else {
        AMREX_ALWAYS_ASSERT(a_levelbcdata->nGrowVect().allGE(ng));
    }
    const MultiFab& bcdata = (a_levelbcdata == nullptr) ? zero : *a_levelbcdata;

    int br_ref_ratio = -1;

    if (amrlev == 0)
    {
        if (needsCoarseDataForBC())
        {
            AMREX_ALWAYS_ASSERT(!hasHiddenDimension());
            br_ref_ratio = m_coarse_data_crse_ratio > 0 ? m_coarse_data_crse_ratio : 2;
            if (m_crse_sol_br[amrlev] == nullptr && br_ref_ratio > 0)
            {
                const int in_rad = 0;
                const int out_rad = 1;
                const int extent_rad = 2;
                const int crse_ratio = br_ref_ratio;
                BoxArray cba = m_grids[amrlev][0];
                cba.coarsen(crse_ratio);
                m_crse_sol_br[amrlev] = std::make_unique<BndryRegister>
                    (cba, m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
            }
            if (m_coarse_data_for_bc != nullptr) {
                AMREX_ALWAYS_ASSERT(m_coarse_data_crse_ratio > 0);
                const Box& cbx = amrex::coarsen(m_geom[0][0].Domain(), m_coarse_data_crse_ratio);
                m_crse_sol_br[amrlev]->copyFrom(*m_coarse_data_for_bc, 0, 0, 0, ncomp,
                                                m_geom[0][0].periodicity(cbx));
            } else {
                m_crse_sol_br[amrlev]->setVal(0.0);
            }
            m_bndry_sol[amrlev]->setBndryValues(*m_crse_sol_br[amrlev], 0,
                                                bcdata, 0, 0, ncomp,
                                                br_ref_ratio, BCRec());
            br_ref_ratio = m_coarse_data_crse_ratio;
        }
        else
        {
            m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,ncomp,BCRec());
            br_ref_ratio = 1;
        }
    }
    else
    {
        m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,ncomp, m_amr_ref_ratio[amrlev-1], BCRec());
        br_ref_ratio = m_amr_ref_ratio[amrlev-1];
    }

    m_bndry_sol[amrlev]->setLOBndryConds(m_lobc, m_hibc, br_ref_ratio, m_coarse_bc_loc);

    const Real* dx = m_geom[amrlev][0].CellSize();
    for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        m_bcondloc[amrlev][mglev]->setLOBndryConds(m_geom[amrlev][mglev], dx,
                                                   m_lobc, m_hibc,
                                                   br_ref_ratio, m_coarse_bc_loc,
                                                   m_domain_bloc_lo, m_domain_bloc_hi);
    }

    if (hasRobinBC()) {
        AMREX_ASSERT(robinbc_a != nullptr && robinbc_b != nullptr && robinbc_f != nullptr);
        m_robin_bcval[amrlev] = std::make_unique<MultiFab>(m_grids[amrlev][0], m_dmap[amrlev][0],
                                                           ncomp*3, 1);
        const Box& domain = m_geom[amrlev][0].Domain();
        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_robin_bcval[amrlev], mfi_info); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            Array4<Real const> const& ra = robinbc_a->const_array(mfi);
            Array4<Real const> const& rb = robinbc_b->const_array(mfi);
            Array4<Real const> const& rf = robinbc_f->const_array(mfi);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& blo = amrex::adjCellLo(vbx, idim);
                const Box& bhi = amrex::adjCellHi(vbx, idim);
                bool outside_domain_lo = !(domain.contains(blo));
                bool outside_domain_hi = !(domain.contains(bhi));
                if ((!outside_domain_lo) && (!outside_domain_hi)) continue;
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    Array4<Real> const& rbc = (*m_robin_bcval[amrlev])[mfi].array(icomp*3);
                    if (m_lobc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_lo)
                    {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(blo, i, j, k,
                        {
                            rbc(i,j,k,0) = ra(i,j,k,icomp);
                            rbc(i,j,k,1) = rb(i,j,k,icomp);
                            rbc(i,j,k,2) = rf(i,j,k,icomp);
                        });
                    }
                    if (m_hibc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_hi)
                    {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bhi, i, j, k,
                        {
                            rbc(i,j,k,0) = ra(i,j,k,icomp);
                            rbc(i,j,k,1) = rb(i,j,k,icomp);
                            rbc(i,j,k,2) = rf(i,j,k,icomp);
                        });
                    }
                }
            }
        }
    }
}

BoxArray
MLCellLinOp::makeNGrids (int grid_size) const
{
    const Box& dombx = m_geom[0].back().Domain();

    const BoxArray& old_ba = m_grids[0].back();
    const int N = old_ba.size();
    Vector<Box> bv;
    bv.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        Box b = old_ba[i];
        b.coarsen(grid_size);
        b.refine(grid_size);
        IntVect sz = b.size();
        const IntVect nblks {AMREX_D_DECL(sz[0]/grid_size, sz[1]/grid_size, sz[2]/grid_size)};

        IntVect big = b.smallEnd() + grid_size - 1;
        b.setBig(big);

#if (AMREX_SPACEDIM == 3)
        for (int kk = 0; kk < nblks[2]; ++kk) {
#endif
#if (AMREX_SPACEDIM >= 2)
            for (int jj = 0; jj < nblks[1]; ++jj) {
#endif
                for (int ii = 0; ii < nblks[0]; ++ii)
                {
                    IntVect shft{AMREX_D_DECL(ii*grid_size,jj*grid_size,kk*grid_size)};
                    Box bb = amrex::shift(b,shft);
                    bb &= dombx;
                    bv.push_back(bb);
                }
#if (AMREX_SPACEDIM >= 2)
            }
#endif
#if (AMREX_SPACEDIM == 3)
        }
#endif
    }

    std::sort(bv.begin(), bv.end());
    bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

    BoxList bl(std::move(bv));

    return BoxArray{std::move(bl)};
}

void
MLCellLinOp::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    const int ncomp = getNComp();
#ifdef AMREX_SOFT_PERF_COUNTERS
    perf_counters.restrict(crse);
#endif
    IntVect ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[cmglev-1];
    amrex::average_down(fine, crse, 0, ncomp, ratio);
}

void
MLCellLinOp::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
#ifdef AMREX_SOFT_PERF_COUNTERS
    perf_counters.interpolate(fine);
#endif

    const int ncomp = getNComp();

    Dim3 ratio3 = {2,2,2};
    IntVect ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[fmglev];
    AMREX_D_TERM(ratio3.x = ratio[0];,
                 ratio3.y = ratio[1];,
                 ratio3.z = ratio[2];);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && fine.isFusingCandidate()) {
        auto const& finema = fine.arrays();
        auto const& crsema = crse.const_arrays();
        ParallelFor(fine, IntVect(0), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            int ic = amrex::coarsen(i,ratio3.x);
            int jc = amrex::coarsen(j,ratio3.y);
            int kc = amrex::coarsen(k,ratio3.z);
            finema[box_no](i,j,k,n) += crsema[box_no](ic,jc,kc,n);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx    = mfi.tilebox();
            Array4<Real const> const& cfab = crse.const_array(mfi);
            Array4<Real> const& ffab = fine.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                int ic = amrex::coarsen(i,ratio3.x);
                int jc = amrex::coarsen(j,ratio3.y);
                int kc = amrex::coarsen(k,ratio3.z);
                ffab(i,j,k,n) += cfab(ic,jc,kc,n);
            });
        }
    }
}

void
MLCellLinOp::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                     const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto amrrr = AMRRefRatio(camrlev);
    const int ncomp = getNComp();
    amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
    amrex::average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
}

void
MLCellLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                    StateMode s_mode, const MLMGBndry* bndry) const
{
    BL_PROFILE("MLCellLinOp::apply()");
    applyBC(amrlev, mglev, in, bc_mode, s_mode, bndry);
#ifdef AMREX_SOFT_PERF_COUNTERS
    perf_counters.apply(out);
#endif
    Fapply(amrlev, mglev, out, in);
}

void
MLCellLinOp::smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                     bool skip_fillboundary) const
{
    BL_PROFILE("MLCellLinOp::smooth()");
    for (int redblack = 0; redblack < 2; ++redblack)
    {
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Solution,
                nullptr, skip_fillboundary);
#ifdef AMREX_SOFT_PERF_COUNTERS
        perf_counters.smooth(sol);
#endif
        Fsmooth(amrlev, mglev, sol, rhs, redblack);
        skip_fillboundary = false;
    }
}

void
MLCellLinOp::updateSolBC (int amrlev, const MultiFab& crse_bcdata) const
{
    BL_PROFILE("MLCellLinOp::updateSolBC()");

    AMREX_ALWAYS_ASSERT(amrlev > 0);
    const int ncomp = getNComp();
    m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_geom[amrlev-1][0].periodicity());
    m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
}

void
MLCellLinOp::updateCorBC (int amrlev, const MultiFab& crse_bcdata) const
{
    BL_PROFILE("MLCellLinOp::updateCorBC()");
    AMREX_ALWAYS_ASSERT(amrlev > 0);
    const int ncomp = getNComp();
    m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp, m_geom[amrlev-1][0].periodicity());
    m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, ncomp, m_amr_ref_ratio[amrlev-1]);
}

void
MLCellLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                               const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::solutionResidual()");
    const int ncomp = getNComp();
    if (crse_bcdata != nullptr) {
        updateSolBC(amrlev, *crse_bcdata);
    }
    const int mglev = 0;
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution,
          m_bndry_sol[amrlev].get());

    AMREX_ALWAYS_ASSERT(resid.nComp() == b.nComp());
    MultiFab::Xpay(resid, Real(-1.0), b, 0, 0, ncomp, 0);
}

void
MLCellLinOp::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                 BCMode bc_mode, const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::correctionResidual()");
    const int ncomp = getNComp();
    if (bc_mode == BCMode::Inhomogeneous)
    {
        if (crse_bcdata)
        {
            AMREX_ALWAYS_ASSERT(mglev == 0);
            AMREX_ALWAYS_ASSERT(amrlev > 0);
            updateCorBC(amrlev, *crse_bcdata);
        }
        apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Correction,
              m_bndry_cor[amrlev].get());
    }
    else
    {
        AMREX_ALWAYS_ASSERT(crse_bcdata == nullptr);
        apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction, nullptr);
    }

    MultiFab::Xpay(resid, Real(-1.0), b, 0, 0, ncomp, 0);
}

void
MLCellLinOp::applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode,
                      const MLMGBndry* bndry, bool skip_fillboundary) const
{
    BL_PROFILE("MLCellLinOp::applyBC()");
    // No coarsened boundary values, cannot apply inhomog at mglev>0.
    BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);
    BL_ASSERT(bndry != nullptr || bc_mode == BCMode::Homogeneous);

    const int ncomp = getNComp();
    const int cross = isCrossStencil();
    const int tensorop = isTensorOp();
    if (!skip_fillboundary) {
        in.FillBoundary(0, ncomp, m_geom[amrlev][mglev].periodicity(),cross);
    }

    int flagbc = bc_mode == BCMode::Inhomogeneous;
    const int imaxorder = maxorder;

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
    const Real dxi = dxinv[0];
    const Real dyi = (AMREX_SPACEDIM >= 2) ? dxinv[1] : Real(1.0);
    const Real dzi = (AMREX_SPACEDIM == 3) ? dxinv[2] : Real(1.0);

    const auto& maskvals = m_maskvals[amrlev][mglev];
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];

    FArrayBox foofab(Box::TheUnitBox(),ncomp);
    const auto& foo = foofab.const_array();

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cross || tensorop || Gpu::notInLaunchRegion(),
                                     "non-cross stencil not support for gpu");

    const int hidden_direction = hiddenDirection();

#ifdef AMREX_USE_GPU
    if ((cross || tensorop) && Gpu::inLaunchRegion())
    {
        Vector<ABCTag> tags;
        tags.reserve(in.local_size()*AMREX_SPACEDIM*ncomp);

        for (MFIter mfi(in); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            const auto& iofab = in.array(mfi);
            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (idim != hidden_direction) {
                    const Orientation olo(idim,Orientation::low);
                    const Orientation ohi(idim,Orientation::high);
                    const auto& bvlo = (bndry != nullptr) ?
                        bndry->bndryValues(olo).const_array(mfi) : foo;
                    const auto& bvhi = (bndry != nullptr) ?
                        bndry->bndryValues(ohi).const_array(mfi) : foo;
                    for (int icomp = 0; icomp < ncomp; ++icomp) {
                        tags.emplace_back(ABCTag{iofab, bvlo, bvhi,
                                                 maskvals[olo].const_array(mfi),
                                                 maskvals[ohi].const_array(mfi),
                                                 bdlv[icomp][olo], bdlv[icomp][ohi],
                                                 amrex::adjCell(vbx,olo),
                                                 bdcv[icomp][olo], bdcv[icomp][ohi],
                                                 vbx.length(idim), icomp, idim});
                    }
                }
            }
        }

        ParallelFor(tags,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, ABCTag const& tag) noexcept
        {
            if (tag.dir == 0)
            {
                mllinop_apply_bc_x(0, i, j, k, tag.blen, tag.fab,
                                   tag.mask_lo, tag.bctype_lo, tag.bcloc_lo, tag.bcval_lo,
                                   imaxorder, dxi, flagbc, tag.comp);
                mllinop_apply_bc_x(1, i+tag.blen+1, j, k, tag.blen, tag.fab,
                                   tag.mask_hi, tag.bctype_hi, tag.bcloc_hi, tag.bcval_hi,
                                   imaxorder, dxi, flagbc, tag.comp);
            }
#if (AMREX_SPACEDIM > 1)
            else
#if (AMREX_SPACEDIM > 2)
            if (tag.dir == 1)
#endif
            {
                mllinop_apply_bc_y(0, i, j, k, tag.blen, tag.fab,
                                   tag.mask_lo, tag.bctype_lo, tag.bcloc_lo, tag.bcval_lo,
                                   imaxorder, dyi, flagbc, tag.comp);
                mllinop_apply_bc_y(1, i, j+tag.blen+1, k, tag.blen, tag.fab,
                                   tag.mask_hi, tag.bctype_hi, tag.bcloc_hi, tag.bcval_hi,
                                   imaxorder, dyi, flagbc, tag.comp);
            }
#if (AMREX_SPACEDIM > 2)
            else {
                mllinop_apply_bc_z(0, i, j, k, tag.blen, tag.fab,
                                   tag.mask_lo, tag.bctype_lo, tag.bcloc_lo, tag.bcval_lo,
                                   imaxorder, dzi, flagbc, tag.comp);
                mllinop_apply_bc_z(1, i, j, k+tag.blen+1, tag.blen, tag.fab,
                                   tag.mask_hi, tag.bctype_hi, tag.bcloc_hi, tag.bcval_hi,
                                   imaxorder, dzi, flagbc, tag.comp);
            }
#endif
#endif
        });
    } else
#endif
    if (cross || tensorop)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(in, mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& vbx   = mfi.validbox();
            const auto& iofab = in.array(mfi);

            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (hidden_direction == idim) continue;
                const Orientation olo(idim,Orientation::low);
                const Orientation ohi(idim,Orientation::high);
                const Box blo = amrex::adjCellLo(vbx, idim);
                const Box bhi = amrex::adjCellHi(vbx, idim);
                const int blen = vbx.length(idim);
                const auto& mlo = maskvals[olo].array(mfi);
                const auto& mhi = maskvals[ohi].array(mfi);
                const auto& bvlo = (bndry != nullptr) ? bndry->bndryValues(olo).const_array(mfi) : foo;
                const auto& bvhi = (bndry != nullptr) ? bndry->bndryValues(ohi).const_array(mfi) : foo;
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    const BoundCond bctlo = bdcv[icomp][olo];
                    const BoundCond bcthi = bdcv[icomp][ohi];
                    const Real bcllo = bdlv[icomp][olo];
                    const Real bclhi = bdlv[icomp][ohi];
                    if (idim == 0) {
                        mllinop_apply_bc_x(0, blo, blen, iofab, mlo,
                                           bctlo, bcllo, bvlo,
                                           imaxorder, dxi, flagbc, icomp);
                        mllinop_apply_bc_x(1, bhi, blen, iofab, mhi,
                                           bcthi, bclhi, bvhi,
                                           imaxorder, dxi, flagbc, icomp);
                    } else if (idim == 1) {
                        mllinop_apply_bc_y(0, blo, blen, iofab, mlo,
                                           bctlo, bcllo, bvlo,
                                           imaxorder, dyi, flagbc, icomp);
                        mllinop_apply_bc_y(1, bhi, blen, iofab, mhi,
                                           bcthi, bclhi, bvhi,
                                           imaxorder, dyi, flagbc, icomp);
                    } else {
                        mllinop_apply_bc_z(0, blo, blen, iofab, mlo,
                                           bctlo, bcllo, bvlo,
                                           imaxorder, dzi, flagbc, icomp);
                        mllinop_apply_bc_z(1, bhi, blen, iofab, mhi,
                                           bcthi, bclhi, bvhi,
                                           imaxorder, dzi, flagbc, icomp);
                    }
                }
            }
        }
    }
    else
    {
#ifdef BL_NO_FORT
        amrex::Abort("amrex_mllinop_apply_bc not available when BL_NO_FORT=TRUE");
#else
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(in, mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& vbx   = mfi.validbox();

            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            const RealTuple & bdl = bdlv[0];
            const BCTuple   & bdc = bdcv[0];

            for (OrientationIter oitr; oitr; ++oitr)
            {
                const Orientation ori = oitr();

                int  cdr = ori;
                Real bcl = bdl[ori];
                int  bct = bdc[ori];

                const FArrayBox& fsfab = (bndry != nullptr) ? bndry->bndryValues(ori)[mfi] : foofab;

                const Mask& m = maskvals[ori][mfi];

                amrex_mllinop_apply_bc(BL_TO_FORTRAN_BOX(vbx),
                                       BL_TO_FORTRAN_ANYD(in[mfi]),
                                       BL_TO_FORTRAN_ANYD(m),
                                       cdr, bct, bcl,
                                       BL_TO_FORTRAN_ANYD(fsfab),
                                       maxorder, dxinv, flagbc, ncomp, cross);
            }
        }
#endif
    }
}

void
MLCellLinOp::reflux (int crse_amrlev,
                     MultiFab& res, const MultiFab& crse_sol, const MultiFab&,
                     MultiFab&, MultiFab& fine_sol, const MultiFab&) const
{
    BL_PROFILE("MLCellLinOp::reflux()");
    YAFluxRegister& fluxreg = m_fluxreg[crse_amrlev];
    fluxreg.reset();

    const int ncomp = getNComp();

    const int fine_amrlev = crse_amrlev+1;

    Real dt = Real(1.0);
    const Real* crse_dx = m_geom[crse_amrlev][0].CellSize();
    const Real* fine_dx = m_geom[fine_amrlev][0].CellSize();

    const int mglev = 0;
    applyBC(fine_amrlev, mglev, fine_sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[fine_amrlev].get());

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FArrayBox,AMREX_SPACEDIM> flux;
        Array<FArrayBox*,AMREX_SPACEDIM> pflux {{ AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) }};
        Array<FArrayBox const*,AMREX_SPACEDIM> cpflux {{ AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) }};

        for (MFIter mfi(crse_sol, mfi_info);  mfi.isValid(); ++mfi)
        {
            if (fluxreg.CrseHasWork(mfi))
            {
                const Box& tbx = mfi.tilebox();
                AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
                             flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
                             flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
                AMREX_D_TERM(Elixir elifx = flux[0].elixir();,
                             Elixir elify = flux[1].elixir();,
                             Elixir elifz = flux[2].elixir(););
                FFlux(crse_amrlev, mfi, pflux, crse_sol[mfi], Location::FaceCentroid);
                fluxreg.CrseAdd(mfi, cpflux, crse_dx, dt, RunOn::Gpu);
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp barrier
#endif

        for (MFIter mfi(fine_sol, mfi_info);  mfi.isValid(); ++mfi)
        {
            if (fluxreg.FineHasWork(mfi))
            {
                const Box& tbx = mfi.tilebox();
                const int face_only = true;
                AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
                             flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
                             flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
                AMREX_D_TERM(Elixir elifx = flux[0].elixir();,
                             Elixir elify = flux[1].elixir();,
                             Elixir elifz = flux[2].elixir(););
                FFlux(fine_amrlev, mfi, pflux, fine_sol[mfi], Location::FaceCentroid, face_only);
                fluxreg.FineAdd(mfi, cpflux, fine_dx, dt, RunOn::Gpu);
            }
        }
    }

    fluxreg.Reflux(res);
}

void
MLCellLinOp::compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                       MultiFab& sol, Location loc) const
{
    BL_PROFILE("MLCellLinOp::compFlux()");

    const int mglev = 0;
    const int ncomp = getNComp();
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[amrlev].get());

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FArrayBox,AMREX_SPACEDIM> flux;
        Array<FArrayBox*,AMREX_SPACEDIM> pflux {{ AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) }};
        for (MFIter mfi(sol, mfi_info);  mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
                         flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
                         flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
            AMREX_D_TERM(Elixir elifx = flux[0].elixir();,
                         Elixir elify = flux[1].elixir();,
                         Elixir elifz = flux[2].elixir(););
            FFlux(amrlev, mfi, pflux, sol[mfi], loc);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& nbx = mfi.nodaltilebox(idim);
                Array4<Real      > dst = fluxes[idim]->array(mfi);
                Array4<Real const> src =  pflux[idim]->array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, ncomp, i, j, k, n,
                {
                    dst(i,j,k,n) = src(i,j,k,n);
                });
            }
        }
    }
}

void
MLCellLinOp::compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                       MultiFab& sol, Location /*loc*/) const
{
    BL_PROFILE("MLCellLinOp::compGrad()");

    if (sol.nComp() > 1)
      amrex::Abort("MLCellLinOp::compGrad called, but only works for single-component solves");

    const int mglev = 0;
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[amrlev].get());

    const int ncomp = getNComp();

    AMREX_D_TERM(const Real dxi = m_geom[amrlev][mglev].InvCellSize(0);,
                 const Real dyi = m_geom[amrlev][mglev].InvCellSize(1);,
                 const Real dzi = m_geom[amrlev][mglev].InvCellSize(2););
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, TilingIfNotGPU());  mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbx = mfi.nodaltilebox(0);,
                     const Box& ybx = mfi.nodaltilebox(1);,
                     const Box& zbx = mfi.nodaltilebox(2););
        const auto& s = sol.array(mfi);
        AMREX_D_TERM(const auto& gx = grad[0]->array(mfi);,
                     const auto& gy = grad[1]->array(mfi);,
                     const auto& gz = grad[2]->array(mfi););

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( xbx, ncomp, i, j, k, n,
        {
            gx(i,j,k,n) = dxi*(s(i,j,k,n) - s(i-1,j,k,n));
        });
#if (AMREX_SPACEDIM >= 2)
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( ybx, ncomp, i, j, k, n,
        {
            gy(i,j,k,n) = dyi*(s(i,j,k,n) - s(i,j-1,k,n));
        });
#endif
#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( zbx, ncomp, i, j, k, n,
        {
            gz(i,j,k,n) = dzi*(s(i,j,k,n) - s(i,j,k-1,n));
        });
#endif
    }
}

void
MLCellLinOp::prepareForSolve ()
{
    BL_PROFILE("MLCellLinOp::prepareForSolve()");

    const int imaxorder = maxorder;
    const int ncomp = getNComp();
    const int hidden_direction = hiddenDirection();
    for (int amrlev = 0;  amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            const auto& bcondloc = *m_bcondloc[amrlev][mglev];
            const auto& maskvals = m_maskvals[amrlev][mglev];

            const Real dxi = m_geom[amrlev][mglev].InvCellSize(0);
            const Real dyi = (AMREX_SPACEDIM >= 2) ? m_geom[amrlev][mglev].InvCellSize(1) : Real(1.0);
            const Real dzi = (AMREX_SPACEDIM == 3) ? m_geom[amrlev][mglev].InvCellSize(2) : Real(1.0);

            BndryRegister& undrrelxr = m_undrrelxr[amrlev][mglev];
            MultiFab foo(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], ncomp, 0, MFInfo().SetAlloc(false));

#ifdef AMREX_USE_EB
            auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
            const FabArray<EBCellFlagFab>* flags =
                (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
            auto area = (factory) ? factory->getAreaFrac()
                : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
#endif

#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
#ifdef AMREX_USE_EB
                if (factory) {
                    Vector<PSEBTag> tags;
                    tags.reserve(foo.local_size()*AMREX_SPACEDIM*ncomp);

                    for (MFIter mfi(foo); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();

                        const auto & bdlv = bcondloc.bndryLocs(mfi);
                        const auto & bdcv = bcondloc.bndryConds(mfi);

                        auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;

                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                        {
                            if (idim != hidden_direction && fabtyp != FabType::covered) {
                                const Orientation olo(idim,Orientation::low);
                                const Orientation ohi(idim,Orientation::high);
                                auto const& ap = (fabtyp == FabType::singlevalued)
                                    ? area[idim]->const_array(mfi) : Array4<Real const>{};
                                for (int icomp = 0; icomp < ncomp; ++icomp) {
                                    tags.emplace_back(PSEBTag{undrrelxr[olo].array(mfi),
                                                              undrrelxr[ohi].array(mfi),
                                                              ap,
                                                              maskvals[olo].const_array(mfi),
                                                              maskvals[ohi].const_array(mfi),
                                                              bdlv[icomp][olo], bdlv[icomp][ohi],
                                                              amrex::adjCell(vbx,olo),
                                                              bdcv[icomp][olo], bdcv[icomp][ohi],
                                                              vbx.length(idim), icomp, idim});
                                }
                            }
                        }
                    }

                    ParallelFor(tags,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, PSEBTag const& tag) noexcept
                    {
                        if (tag.ap) {
                            if (tag.dir == 0)
                            {
                                mllinop_comp_interp_coef0_x_eb
                                    (0, i           , j, k, tag.blen, tag.flo, tag.mlo, tag.ap,
                                     tag.bctlo, tag.bcllo, imaxorder, dxi, tag.comp);
                                mllinop_comp_interp_coef0_x_eb
                                    (1, i+tag.blen+1, j, k, tag.blen, tag.fhi, tag.mhi, tag.ap,
                                     tag.bcthi, tag.bclhi, imaxorder, dxi, tag.comp);
                            }
#if (AMREX_SPACEDIM > 1)
                            else
#if (AMREX_SPACEDIM > 2)
                            if (tag.dir == 1)
#endif
                            {
                                mllinop_comp_interp_coef0_y_eb
                                    (0, i, j           , k, tag.blen, tag.flo, tag.mlo, tag.ap,
                                     tag.bctlo, tag.bcllo, imaxorder, dyi, tag.comp);
                                mllinop_comp_interp_coef0_y_eb
                                    (1, i, j+tag.blen+1, k, tag.blen, tag.fhi, tag.mhi, tag.ap,
                                     tag.bcthi, tag.bclhi, imaxorder, dyi, tag.comp);
                            }
#if (AMREX_SPACEDIM > 2)
                            else {
                                mllinop_comp_interp_coef0_z_eb
                                    (0, i, j, k           , tag.blen, tag.flo, tag.mlo, tag.ap,
                                     tag.bctlo, tag.bcllo, imaxorder, dzi, tag.comp);
                                mllinop_comp_interp_coef0_z_eb
                                    (1, i, j, k+tag.blen+1, tag.blen, tag.fhi, tag.mhi, tag.ap,
                                     tag.bcthi, tag.bclhi, imaxorder, dzi, tag.comp);
                            }
#endif
#endif
                        } else {
                            if (tag.dir == 0)
                            {
                                mllinop_comp_interp_coef0_x
                                    (0, i           , j, k, tag.blen, tag.flo, tag.mlo,
                                     tag.bctlo, tag.bcllo, imaxorder, dxi, tag.comp);
                                mllinop_comp_interp_coef0_x
                                    (1, i+tag.blen+1, j, k, tag.blen, tag.fhi, tag.mhi,
                                     tag.bcthi, tag.bclhi, imaxorder, dxi, tag.comp);
                            }
#if (AMREX_SPACEDIM > 1)
                            else
#if (AMREX_SPACEDIM > 2)
                            if (tag.dir == 1)
#endif
                            {
                                mllinop_comp_interp_coef0_y
                                    (0, i, j           , k, tag.blen, tag.flo, tag.mlo,
                                     tag.bctlo, tag.bcllo, imaxorder, dyi, tag.comp);
                                mllinop_comp_interp_coef0_y
                                    (1, i, j+tag.blen+1, k, tag.blen, tag.fhi, tag.mhi,
                                     tag.bcthi, tag.bclhi, imaxorder, dyi, tag.comp);
                            }
#if (AMREX_SPACEDIM > 2)
                            else {
                                mllinop_comp_interp_coef0_z
                                    (0, i, j, k           , tag.blen, tag.flo, tag.mlo,
                                     tag.bctlo, tag.bcllo, imaxorder, dzi, tag.comp);
                                mllinop_comp_interp_coef0_z
                                    (1, i, j, k+tag.blen+1, tag.blen, tag.fhi, tag.mhi,
                                     tag.bcthi, tag.bclhi, imaxorder, dzi, tag.comp);
                            }
#endif
#endif
                        }
                    });
                } else
#endif
                {
                    Vector<PSTag> tags;
                    tags.reserve(foo.local_size()*AMREX_SPACEDIM*ncomp);

                    for (MFIter mfi(foo); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();

                        const auto & bdlv = bcondloc.bndryLocs(mfi);
                        const auto & bdcv = bcondloc.bndryConds(mfi);

                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                        {
                            if (idim != hidden_direction) {
                                const Orientation olo(idim,Orientation::low);
                                const Orientation ohi(idim,Orientation::high);
                                for (int icomp = 0; icomp < ncomp; ++icomp) {
                                    tags.emplace_back(PSTag{undrrelxr[olo].array(mfi),
                                                            undrrelxr[ohi].array(mfi),
                                                            maskvals[olo].const_array(mfi),
                                                            maskvals[ohi].const_array(mfi),
                                                            bdlv[icomp][olo], bdlv[icomp][ohi],
                                                            amrex::adjCell(vbx,olo),
                                                            bdcv[icomp][olo], bdcv[icomp][ohi],
                                                            vbx.length(idim), icomp, idim});
                                }
                            }
                        }
                    }

                    ParallelFor(tags,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, PSTag const& tag) noexcept
                    {
                        if (tag.dir == 0)
                        {
                            mllinop_comp_interp_coef0_x
                                (0, i           , j, k, tag.blen, tag.flo, tag.mlo,
                                 tag.bctlo, tag.bcllo, imaxorder, dxi, tag.comp);
                            mllinop_comp_interp_coef0_x
                                (1, i+tag.blen+1, j, k, tag.blen, tag.fhi, tag.mhi,
                                 tag.bcthi, tag.bclhi, imaxorder, dxi, tag.comp);
                        }
#if (AMREX_SPACEDIM > 1)
                        else
#if (AMREX_SPACEDIM > 2)
                        if (tag.dir == 1)
#endif
                        {
                            mllinop_comp_interp_coef0_y
                                (0, i, j           , k, tag.blen, tag.flo, tag.mlo,
                                 tag.bctlo, tag.bcllo, imaxorder, dyi, tag.comp);
                            mllinop_comp_interp_coef0_y
                                (1, i, j+tag.blen+1, k, tag.blen, tag.fhi, tag.mhi,
                                 tag.bcthi, tag.bclhi, imaxorder, dyi, tag.comp);
                        }
#if (AMREX_SPACEDIM > 2)
                        else {
                            mllinop_comp_interp_coef0_z
                                (0, i, j, k           , tag.blen, tag.flo, tag.mlo,
                                 tag.bctlo, tag.bcllo, imaxorder, dzi, tag.comp);
                            mllinop_comp_interp_coef0_z
                                (1, i, j, k+tag.blen+1, tag.blen, tag.fhi, tag.mhi,
                                 tag.bcthi, tag.bclhi, imaxorder, dzi, tag.comp);
                        }
#endif
#endif
                    });
                }
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(foo, MFItInfo{}.SetDynamic(true)); mfi.isValid(); ++mfi)
                {
                    const Box& vbx = mfi.validbox();

                    const auto & bdlv = bcondloc.bndryLocs(mfi);
                    const auto & bdcv = bcondloc.bndryConds(mfi);

#ifdef AMREX_USE_EB
                    auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;
#endif
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        if (idim == hidden_direction) continue;
                        const Orientation olo(idim,Orientation::low);
                        const Orientation ohi(idim,Orientation::high);
                        const Box blo = amrex::adjCellLo(vbx, idim);
                        const Box bhi = amrex::adjCellHi(vbx, idim);
                        const int blen = vbx.length(idim);
                        const auto& mlo = maskvals[olo].array(mfi);
                        const auto& mhi = maskvals[ohi].array(mfi);
                        const auto& flo = undrrelxr[olo].array(mfi);
                        const auto& fhi = undrrelxr[ohi].array(mfi);
                        for (int icomp = 0; icomp < ncomp; ++icomp) {
                            const BoundCond bctlo = bdcv[icomp][olo];
                            const BoundCond bcthi = bdcv[icomp][ohi];
                            const Real bcllo = bdlv[icomp][olo];
                            const Real bclhi = bdlv[icomp][ohi];
#ifdef AMREX_USE_EB
                            if (fabtyp == FabType::singlevalued) {
                                Array4<Real const> const& ap = area[idim]->const_array(mfi);
                                if (idim == 0) {
                                    mllinop_comp_interp_coef0_x_eb
                                        (0, blo, blen, flo, mlo, ap, bctlo, bcllo,
                                         imaxorder, dxi, icomp);
                                    mllinop_comp_interp_coef0_x_eb
                                        (1, bhi, blen, fhi, mhi, ap, bcthi, bclhi,
                                         imaxorder, dxi, icomp);
                                } else if (idim == 1) {
                                    mllinop_comp_interp_coef0_y_eb
                                        (0, blo, blen, flo, mlo, ap, bctlo, bcllo,
                                         imaxorder, dyi, icomp);
                                    mllinop_comp_interp_coef0_y_eb
                                        (1, bhi, blen, fhi, mhi, ap, bcthi, bclhi,
                                         imaxorder, dyi, icomp);
                                } else {
                                    mllinop_comp_interp_coef0_z_eb
                                        (0, blo, blen, flo, mlo, ap, bctlo, bcllo,
                                         imaxorder, dzi, icomp);
                                    mllinop_comp_interp_coef0_z_eb
                                        (1, bhi, blen, fhi, mhi, ap, bcthi, bclhi,
                                         imaxorder, dzi, icomp);
                                }
                            } else if (fabtyp == FabType::regular)
#endif
                            {
                                if (idim == 0) {
                                    mllinop_comp_interp_coef0_x
                                        (0, blo, blen, flo, mlo, bctlo, bcllo,
                                         imaxorder, dxi, icomp);
                                    mllinop_comp_interp_coef0_x
                                        (1, bhi, blen, fhi, mhi, bcthi, bclhi,
                                         imaxorder, dxi, icomp);
                                } else if (idim == 1) {
                                    mllinop_comp_interp_coef0_y
                                        (0, blo, blen, flo, mlo, bctlo, bcllo,
                                         imaxorder, dyi, icomp);
                                    mllinop_comp_interp_coef0_y
                                        (1, bhi, blen, fhi, mhi, bcthi, bclhi,
                                         imaxorder, dyi, icomp);
                                } else {
                                    mllinop_comp_interp_coef0_z
                                        (0, blo, blen, flo, mlo, bctlo, bcllo,
                                         imaxorder, dzi, icomp);
                                    mllinop_comp_interp_coef0_z
                                        (1, bhi, blen, fhi, mhi, bcthi, bclhi,
                                         imaxorder, dzi, icomp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

Real
MLCellLinOp::xdoty (int /*amrlev*/, int /*mglev*/, const MultiFab& x, const MultiFab& y, bool local) const
{
    const int ncomp = getNComp();
    const int nghost = 0;
    Real result = MultiFab::Dot(x,0,y,0,ncomp,nghost,true);
    if (!local) {
        ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    }
    return result;
}

MLCellLinOp::BndryCondLoc::BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm, int ncomp)
    : bcond(ba, dm),
      bcloc(ba, dm),
      bctl(ba, dm),
      bctl_dv(bctl.local_size()*ncomp),
      m_ncomp(ncomp)
{
    auto dp = bctl_dv.data();
    for (MFIter mfi(bcloc); mfi.isValid(); ++mfi) {
        bcond[mfi].resize(ncomp);
        bcloc[mfi].resize(ncomp);
        bctl[mfi] = dp;
        dp += ncomp;
    }
}

void
MLCellLinOp::BndryCondLoc::setLOBndryConds (const Geometry& geom, const Real* dx,
                                            const Vector<Array<BCType,AMREX_SPACEDIM> >& lobc,
                                            const Vector<Array<BCType,AMREX_SPACEDIM> >& hibc,
                                            int ratio, const RealVect& interior_bloc,
                                            const Array<Real,AMREX_SPACEDIM>& domain_bloc_lo,
                                            const Array<Real,AMREX_SPACEDIM>& domain_bloc_hi)
{
    const Box& domain = geom.Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(bcloc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            RealTuple & bloc  = bcloc[mfi][icomp];
            BCTuple   & bctag = bcond[mfi][icomp];
            MLMGBndry::setBoxBC(bloc, bctag, bx, domain, lobc[icomp], hibc[icomp],
                                dx, ratio, interior_bloc, domain_bloc_lo, domain_bloc_hi,
                                geom.isPeriodicArray());
        }
    }

    Gpu::PinnedVector<GpuArray<BCTL,2*AMREX_SPACEDIM> > hv;
    hv.reserve(bctl_dv.size());
    for (MFIter mfi(bctl); mfi.isValid(); ++mfi)
    {
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            GpuArray<BCTL,2*AMREX_SPACEDIM> tmp;
            for (int m = 0; m < 2*AMREX_SPACEDIM; ++m) {
                tmp[m].type = bcond[mfi][icomp][m];
                tmp[m].location = bcloc[mfi][icomp][m];
            }
            hv.push_back(std::move(tmp));
        }
    }
    Gpu::copyAsync(Gpu::hostToDevice, hv.begin(), hv.end(), bctl_dv.begin());
    Gpu::streamSynchronize();
}

void
MLCellLinOp::applyMetricTerm (int amrlev, int mglev, Any& rhs) const
{
    amrex::ignore_unused(amrlev,mglev,rhs);
#if (AMREX_SPACEDIM != 3)

    if (!m_has_metric_term) return;

    AMREX_ASSERT(rhs.is<MultiFab>());
    applyMetricTermToMF(amrlev, mglev, rhs.get<MultiFab>());
#endif
}

void
MLCellLinOp::applyMetricTermToMF (int amrlev, int mglev, MultiFab& rhs) const
{
    amrex::ignore_unused(amrlev,mglev,rhs);
#if (AMREX_SPACEDIM != 3)

    if (!m_has_metric_term) return;

    const int ncomp = rhs.nComp();

    bool cc = rhs.ixType().cellCentered(0);

    const Geometry& geom = m_geom[amrlev][mglev];
    const Real dx = geom.CellSize(0);
    const Real probxlo = geom.ProbLo(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        Array4<Real> const& rhsarr = rhs.array(mfi);
#if (AMREX_SPACEDIM == 1)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real rc = probxlo + (i+Real(0.5))*dx;
                rhsarr(i,j,k,n) *= rc*rc;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real re = probxlo + i*dx;
                rhsarr(i,j,k,n) *= re*re;
            });
        }
#elif (AMREX_SPACEDIM == 2)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real rc = probxlo + (i+Real(0.5))*dx;
                rhsarr(i,j,k,n) *= rc;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real re = probxlo + i*dx;
                rhsarr(i,j,k,n) *= re;
            });
        }
#endif
    }
#endif
}

void
MLCellLinOp::unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const
{
    amrex::ignore_unused(amrlev,mglev,rhs);
#if (AMREX_SPACEDIM != 3)

    if (!m_has_metric_term) return;

    const int ncomp = rhs.nComp();

    bool cc = rhs.ixType().cellCentered(0);

    const Geometry& geom = m_geom[amrlev][mglev];
    const Real dx = geom.CellSize(0);
    const Real probxlo = geom.ProbLo(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        Array4<Real> const& rhsarr = rhs.array(mfi);
#if (AMREX_SPACEDIM == 1)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real rcinv = Real(1.0)/(probxlo + (i+Real(0.5))*dx);
                rhsarr(i,j,k,n) *= rcinv*rcinv;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real re = probxlo + i*dx;
                Real reinv = (re==Real(0.0)) ? Real(0.0) : Real(1.)/re;
                rhsarr(i,j,k,n) *= reinv*reinv;
            });
        }
#elif (AMREX_SPACEDIM == 2)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real rcinv = Real(1.0)/(probxlo + (i+Real(0.5))*dx);
                rhsarr(i,j,k,n) *= rcinv;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                Real re = probxlo + i*dx;
                Real reinv = (re==Real(0.0)) ? Real(0.0) : Real(1.)/re;
                rhsarr(i,j,k,n) *= reinv;
            });
        }
#endif
    }
#endif
}

void
MLCellLinOp::update ()
{
    if (MLLinOp::needsUpdate()) MLLinOp::update();
}

void
MLCellLinOp::computeVolInv () const
{
    if (!m_volinv.empty()) return;

    m_volinv.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_volinv[amrlev].resize(NMGLevels(amrlev));
    }

    // We don't need to compute for every level

    auto f = [&] (int amrlev, int mglev) {
#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(Factory(amrlev,mglev));
        if (factory)
        {
            const MultiFab& vfrac = factory->getVolFrac();
            m_volinv[amrlev][mglev] = vfrac.sum(0,true);
        }
        else
#endif
        {
            m_volinv[amrlev][mglev]
                = Real(1.0 / compactify(Geom(amrlev,mglev).Domain()).d_numPts());
        }
    };

    // amrlev = 0, mglev = 0
    f(0,0);

    int mgbottom = NMGLevels(0)-1;
    f(0,mgbottom);

#ifdef AMREX_USE_EB
    Real temp1, temp2;
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(Factory(0,0));
    if (factory)
    {
        ParallelAllReduce::Sum<Real>({m_volinv[0][0], m_volinv[0][mgbottom]},
                                     ParallelContext::CommunicatorSub());
        temp1 = Real(1.0)/m_volinv[0][0];
        temp2 = Real(1.0)/m_volinv[0][mgbottom];
    }
    else
    {
        temp1 = m_volinv[0][0];
        temp2 = m_volinv[0][mgbottom];
    }
    m_volinv[0][0] = temp1;
    m_volinv[0][mgbottom] = temp2;
#endif
}

Vector<Real>
MLCellLinOp::getSolvabilityOffset (int amrlev, int mglev, Any const& a_rhs) const
{
    AMREX_ASSERT(a_rhs.is<MultiFab>());
    auto const& rhs = a_rhs.get<MultiFab>();

    computeVolInv();

    const int ncomp = getNComp();
    Vector<Real> offset(ncomp);

#ifdef AMREX_USE_EB
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(Factory(amrlev,mglev));
    if (factory)
    {
        const MultiFab& vfrac = factory->getVolFrac();
        for (int c = 0; c < ncomp; ++c) {
            offset[c] = MultiFab::Dot(rhs, c, vfrac, 0, 1, 0, true) * m_volinv[amrlev][mglev];
        }
    }
    else
#endif
    {
        for (int c = 0; c < ncomp; ++c) {
            offset[c] = rhs.sum(c,true) * m_volinv[amrlev][mglev];
        }
    }

    ParallelAllReduce::Sum(offset.data(), ncomp, ParallelContext::CommunicatorSub());

    return offset;
}

Real
MLCellLinOp::AnyNormInfMask (int amrlev, Any const& a, bool local) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    auto& mf = a.get<MultiFab>();

    const int finest_level = NAMRLevels() - 1;
    Real norm = 0._rt;
#ifdef AMREX_USE_EB
    const int ncomp = getNComp();
    if (! mf.isAllRegular()) {
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(Factory(amrlev));
        const MultiFab& vfrac = factory->getVolFrac();
        if (amrlev == finest_level) {
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
                auto const& ma = mf.const_arrays();
                auto const& vfrac_ma = vfrac.const_arrays();
                norm = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                                 mf, IntVect(0), ncomp,
                                 [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                                     -> GpuTuple<Real>
                                 {
                                     return amrex::Math::abs(ma[box_no](i,j,k,n)
                                                             *vfrac_ma[box_no](i,j,k));
                                 });
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:norm)
#endif
                for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                    Box const& bx = mfi.tilebox();
                    auto const& fab = mf.const_array(mfi);
                    auto const& v = vfrac.const_array(mfi);
                    AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                    {
                        norm = std::max(norm, amrex::Math::abs(fab(i,j,k,n)*v(i,j,k)));
                    });
                }
            }
        } else {
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
                auto const& ma = mf.const_arrays();
                auto const& mask_ma = m_norm_fine_mask[amrlev]->const_arrays();
                auto const& vfrac_ma = vfrac.const_arrays();
                norm = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                                 mf, IntVect(0), ncomp,
                                 [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                                     -> GpuTuple<Real>
                                 {
                                     if (mask_ma[box_no](i,j,k)) {
                                         return amrex::Math::abs(ma[box_no](i,j,k,n)
                                                                 *vfrac_ma[box_no](i,j,k));
                                     } else {
                                         return Real(0.0);
                                     }
                                 });
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:norm)
#endif
                for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                    Box const& bx = mfi.tilebox();
                    auto const& fab = mf.const_array(mfi);
                    auto const& mask = m_norm_fine_mask[amrlev]->const_array(mfi);
                    auto const& v = vfrac.const_array(mfi);
                    AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                    {
                        if (mask(i,j,k)) {
                            norm = std::max(norm, amrex::Math::abs(fab(i,j,k,n)*v(i,j,k)));
                        }
                    });
                }
            }
        }
    } else
#endif
    {
        iMultiFab const* fine_mask = (amrlev == finest_level)
            ? nullptr : m_norm_fine_mask[amrlev].get();
        norm = MFNormInf(mf, fine_mask, true);
    }

    if (!local) ParallelAllReduce::Max(norm, ParallelContext::CommunicatorSub());
    return norm;
}

void
MLCellLinOp::AnyAvgDownResAmr (int clev, Any& cres, Any const& fres) const
{
    AMREX_ASSERT(cres.is<MultiFab>() && fres.is<MultiFab>());
#ifdef AMREX_USE_EB
    amrex::EB_average_down
#else
    amrex::average_down
#endif
        (fres.get<MultiFab>(), cres.get<MultiFab>(), 0, getNComp(), AMRRefRatio(clev));
}

void
MLCellLinOp::AnyInterpolationAmr (int famrlev, Any& a_fine, const Any& a_crse,
                                  IntVect const& /*nghost*/) const
{
    AMREX_ASSERT(a_fine.is<MultiFab>());
    MultiFab& fine = a_fine.get<MultiFab>();
    MultiFab const& crse = a_crse.get<MultiFab>();

    const int ncomp = getNComp();
    const int refratio = AMRRefRatio(famrlev-1);

#ifdef AMREX_USE_EB
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(Factory(famrlev));
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& ff = fine.array(mfi);
        Array4<Real const> const& cc = crse.const_array(mfi);
#ifdef AMREX_USE_EB
        bool call_lincc;
        if (factory)
        {
            const auto& flag = (*flags)[mfi];
            if (flag.getType(amrex::grow(bx,1)) == FabType::regular) {
                call_lincc = true;
            } else {
                Array4<EBCellFlag const> const& flg = flag.const_array();
                switch(refratio) {
                case 2:
                {
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                    {
                        mlmg_eb_cc_interp_r<2>(tbx, ff, cc, flg, ncomp);
                    });
                    break;
                }
                case 4:
                {
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                    {
                        mlmg_eb_cc_interp_r<4>(tbx, ff, cc, flg, ncomp);
                    });
                    break;
                }
                default:
                    amrex::Abort("mlmg_eb_cc_interp: only refratio 2 and 4 are supported");
                }

                call_lincc = false;
            }
        }
        else
        {
            call_lincc = true;
        }
#else
        const bool call_lincc = true;
#endif
        if (call_lincc)
        {
            switch(refratio) {
            case 2:
            {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    mlmg_lin_cc_interp_r2(tbx, ff, cc, ncomp);
                });
                break;
            }
            case 4:
            {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    mlmg_lin_cc_interp_r4(tbx, ff, cc, ncomp);
                });
                break;
            }
            default:
                amrex::Abort("mlmg_lin_cc_interp: only refratio 2 and 4 are supported");
            }
        }
    }
}

void
MLCellLinOp::interpAssign (int amrlev, int fmglev, MultiFab& fine, MultiFab& crse) const
{
    const int ncomp = getNComp();

    const Geometry& crse_geom = Geom(amrlev,fmglev+1);
    const IntVect refratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[fmglev];
    const IntVect ng = crse.nGrowVect();

    MultiFab cfine;
    const MultiFab* cmf;

    if (amrex::isMFIterSafe(crse, fine))
    {
        crse.FillBoundary(crse_geom.periodicity());
        cmf = &crse;
    }
    else
    {
        BoxArray cba = fine.boxArray();
        cba.coarsen(refratio);
        cfine.define(cba, fine.DistributionMap(), ncomp, ng);
        cfine.setVal(0.0);
        cfine.ParallelCopy(crse, 0, 0, ncomp, IntVect(0), ng, crse_geom.periodicity());
        cmf = & cfine;
    }

    bool isEB = fine.hasEBFabFactory();
    ignore_unused(isEB);

#ifdef AMREX_USE_EB
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(fine.Factory()));
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& ff = fine.array(mfi);
        const auto& cc = cmf->array(mfi);
#ifdef AMREX_USE_EB
        bool call_lincc;
        if (isEB)
        {
            const auto& flag = (*flags)[mfi];
            if (flag.getType(amrex::grow(bx,1)) == FabType::regular) {
                call_lincc = true;
            } else {
                Array4<EBCellFlag const> const& flg = flag.const_array();
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    mlmg_eb_cc_interp_r<2>(tbx, ff, cc, flg, ncomp);
                });

                call_lincc = false;
            }
        }
        else
        {
            call_lincc = true;
        }
#else
        const bool call_lincc = true;
#endif
        if (call_lincc)
        {
#if (AMREX_SPACEDIM == 3)
            if (hasHiddenDimension()) {
                Box const& bx_2d = compactify(bx);
                auto const& ff_2d = compactify(ff);
                auto const& cc_2d = compactify(cc);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx_2d, tbx,
                {
                    TwoD::mlmg_lin_cc_interp_r2(tbx, ff_2d, cc_2d, ncomp);
                });
            } else
#endif
            {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
                {
                    mlmg_lin_cc_interp_r2(tbx, ff, cc, ncomp);
                });
            }
        }
    }
}

void
MLCellLinOp::AnyAverageDownAndSync (Vector<Any>& sol) const
{
    AMREX_ASSERT(sol[0].is<MultiFab>());

    int ncomp = getNComp();
    for (int falev = NAMRLevels()-1; falev > 0; --falev)
    {
#ifdef AMREX_USE_EB
        amrex::EB_average_down(sol[falev  ].get<MultiFab>(),
                               sol[falev-1].get<MultiFab>(), 0, ncomp, AMRRefRatio(falev-1));
#else
        amrex::average_down(sol[falev  ].get<MultiFab>(),
                            sol[falev-1].get<MultiFab>(), 0, ncomp, AMRRefRatio(falev-1));
#endif
    }
}

void
MLCellLinOp::fixSolvabilityByOffset (int amrlev, int mglev, Any& a_rhs,
                                     Vector<Real> const& offset) const
{
    amrex::ignore_unused(amrlev, mglev);
    AMREX_ASSERT(a_rhs.is<MultiFab>());
    auto& rhs = a_rhs.get<MultiFab>();

    const int ncomp = getNComp();
    for (int c = 0; c < ncomp; ++c) {
        rhs.plus(-offset[c], c, 1);
    }
#ifdef AMREX_USE_EB
    if (rhs.hasEBFabFactory()) {
        Vector<Real> val(ncomp, 0.0_rt);
        amrex::EB_set_covered(rhs, 0, ncomp, val);
    }
#endif
}

}
