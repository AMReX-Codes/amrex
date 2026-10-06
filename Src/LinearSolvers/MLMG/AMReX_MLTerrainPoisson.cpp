#include <AMReX_MLTerrainPoisson.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <map>
#include <set>

namespace amrex {

MLTerrainPoisson::MLTerrainPoisson (const Vector<Geometry>& a_geom,
                                    const Vector<BoxArray>& a_grids,
                                    const Vector<DistributionMapping>& a_dmap,
                                    const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

void
MLTerrainPoisson::define (const Vector<Geometry>& a_geom,
                          const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap,
                          const LPInfo& a_info)
{
    BL_PROFILE("MLTerrainPoisson::define()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_geom.size() == 1,
                                     "MLTerrainPoisson: only one AMR level is supported");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!a_info.hasHiddenDimension() || a_info.hidden_direction < 2,
                                     "MLTerrainPoisson: z cannot be the hidden direction");

    // The line smoother solves z exactly, so coarsen horizontally only, down
    // to single columns, and keep agglomerated boxes whole in z.
    LPInfo lpinfo = a_info;
    lpinfo.setSemicoarsening(true).setSemicoarseningDirection(2)
        .setMaxSemicoarseningLevel(lpinfo.max_coarsening_level);
    mg_box_min_width = 1;
    mg_domain_min_width = 1;
    mg_agg_no_split_direction = 2;
    mg_odd_coarsening = true;
    mg_independent_coarsening = true;

    MLCellLinOpT<MultiFab>::define(a_geom, a_grids, a_dmap, lpinfo);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_domain_covered[0],
                                     "MLTerrainPoisson: the grids must cover the domain");

    m_is_singular.assign(1, 0);

    const int nmglevs = m_num_mg_levels[0];
    m_zphys.define(amrex::convert(m_grids[0][0], IntVect(1)), m_dmap[0][0], 1, 1);
    m_area.clear();
    m_area.resize(nmglevs);
    m_rx.clear();
    m_ry.clear();
    m_zf.clear();
    m_rx.resize(nmglevs);
    m_ry.resize(nmglevs);
    m_zf.resize(nmglevs);
    m_bc_tags.clear();
    m_bc_tags.resize(nmglevs);
    m_smooth_res.clear();
    m_smooth_res.resize(nmglevs);
    for (int mglev = 0; mglev < nmglevs; ++mglev) {
        auto const& ba = m_grids[0][mglev];
        auto const& dm = m_dmap[0][mglev];
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_area[mglev][idim].define(amrex::convert(ba, IntVect::TheDimensionVector(idim)),
                                       dm, 1, 0);
        }
        m_rx[mglev].define(amrex::convert(ba, IntVect(1,0,1)), dm, 1, 0);
        m_ry[mglev].define(amrex::convert(ba, IntVect(0,1,1)), dm, 1, 0);
        m_zf[mglev].define(amrex::convert(ba, IntVect(0,0,1)), dm, 3, 0);
    }
    defineLineSolvers();
    m_detJ.clear();
    m_needs_update = true;
}

void
MLTerrainPoisson::setZSplitSolver (ZSplitSolver a_solver)
{
    if (a_solver != m_zsplit) {
        m_zsplit = a_solver;
        if (!m_rx.empty()) { defineLineSolvers(); }
    }
}

void
MLTerrainPoisson::defineLineSolvers ()
{
    const int nmglevs = m_num_mg_levels[0];
    m_column_lu.clear();
    m_column_slab.clear();
    m_spike.clear();
    m_spike_end.clear();
    m_spike_buf.clear();
    m_spike_sum.clear();
    m_zcut.clear();
    m_spike_lu.clear();
    m_col_lu.clear();
    m_col_res.clear();
    m_col_cor.clear();
    m_column_lu.resize(nmglevs);
    m_column_slab.resize(nmglevs);
    m_spike.resize(nmglevs);
    m_spike_end.resize(nmglevs);
    m_spike_buf.resize(nmglevs);
    m_spike_sum.resize(nmglevs);
    m_zcut.resize(nmglevs);
    m_spike_lu.resize(nmglevs);
    m_col_lu.resize(nmglevs);
    m_col_res.resize(nmglevs);
    m_col_cor.resize(nmglevs);
    for (int mglev = 0; mglev < nmglevs; ++mglev) {
        auto const& ba = m_grids[0][mglev];
        auto const& dm = m_dmap[0][mglev];
        // Boxes split in z: whole columns, each on the rank owning most of it.
        Box const& domain = m_geom[0][mglev].Domain();
        bool split = false;
        for (int i = 0, N = static_cast<int>(ba.size()); i < N; ++i) {
            split = split || ba[i].length(2) < domain.length(2);
        }
        if (split) {
            BoxList bl;
            bl.reserve(ba.size());
            for (int i = 0, N = static_cast<int>(ba.size()); i < N; ++i) {
                bl.push_back(Box(ba[i]).setRange(2, domain.smallEnd(2), domain.length(2)));
            }
            BoxArray cba(std::move(bl));
            cba.removeOverlap(false);
            // Ties go to the rank with the fewest cells so far.
            Vector<int> pmap(cba.size());
            Vector<Long> load(ParallelDescriptor::NProcs(), 0);
            for (int ic = 0, N = static_cast<int>(cba.size()); ic < N; ++ic) {
                std::map<int,Long> npts;
                for (auto const& is : ba.intersections(cba[ic])) {
                    npts[dm[is.first]] += is.second.numPts();
                }
                Long nmax = -1;
                for (auto const& [rank, n] : npts) {
                    if (n > nmax || (n == nmax && load[rank] < load[pmap[ic]])) {
                        nmax = n;
                        pmap[ic] = rank;
                    }
                }
                load[pmap[ic]] += cba[ic].numPts();
            }
            DistributionMapping cdm(std::move(pmap));

            // Spikes need the same z cuts in every column; else whole columns.
            std::set<std::pair<int,int>> zr;
            for (int i = 0, N = static_cast<int>(ba.size()); i < N; ++i) {
                zr.emplace(ba[i].smallEnd(2), ba[i].bigEnd(2));
            }
            bool regular = true;
            Vector<int> zcut{domain.smallEnd(2)};
            for (auto const& [zlo, zhi] : zr) {
                regular = regular && zlo == zcut.back();
                zcut.push_back(zhi+1);
            }
            regular = regular && zcut.back() == domain.bigEnd(2)+1;
            if (regular && m_zsplit == ZSplitSolver::Spike) {
                int const np = static_cast<int>(zr.size());
                m_column_lu[mglev].define(ba, dm, 3, 0);
                m_spike[mglev].define(ba, dm, 2, 0);
                // x-y slabs with k = block index: one per box, and per
                // column with all p blocks.
                BoxList sl;
                for (int i = 0, N = static_cast<int>(ba.size()); i < N; ++i) {
                    int const q = static_cast<int>(std::find(zcut.begin(), zcut.end(),
                                                             ba[i].smallEnd(2)) - zcut.begin());
                    sl.push_back(Box(ba[i]).setRange(2, q));
                }
                m_spike_buf[mglev].define(BoxArray(std::move(sl)), dm, 2, 0);
                BoxList cl;
                for (int i = 0, N = static_cast<int>(cba.size()); i < N; ++i) {
                    cl.push_back(Box(cba[i]).setRange(2, 0, np));
                }
                BoxArray const csba(std::move(cl));
                m_spike_end[mglev].define(csba, cdm, 4, 0);
                m_spike_sum[mglev].define(csba, cdm, 2, 0);
                m_spike_lu[mglev].define(csba, cdm, 8, 0);
                m_zcut[mglev] = std::move(zcut);
            } else {
                m_col_lu[mglev].define(cba, cdm, 3, 0);
                m_col_res[mglev].define(cba, cdm, 1, 0);
                m_col_cor[mglev].define(cba, cdm, 1, 0);
            }
        } else {
            m_column_lu[mglev].define(ba, dm, 3, 0);
        }
    }
    m_needs_update = true;
}

void
MLTerrainPoisson::setZPhys (int amrlev, MultiFab const& z_phys_nd)
{
    AMREX_ALWAYS_ASSERT(amrlev == 0 && z_phys_nd.nGrowVect().allGE(1) &&
                        z_phys_nd.ixType().nodeCentered());
    m_zphys.LocalCopy(z_phys_nd, 0, 0, 1, IntVect(1));
    m_needs_update = true;
}

void
MLTerrainPoisson::setAreas (int amrlev, Array<MultiFab const*,AMREX_SPACEDIM> const& area)
{
    AMREX_ALWAYS_ASSERT(amrlev == 0);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT(area[idim]->ixType() == m_area[0][idim].ixType());
        m_area[0][idim].LocalCopy(*area[idim], 0, 0, 1, IntVect(0));
    }
    m_needs_update = true;
}

void
MLTerrainPoisson::setDetJ (int amrlev, MultiFab const& detJ)
{
    AMREX_ALWAYS_ASSERT(amrlev == 0 && detJ.ixType().cellCentered());
    m_detJ.define(m_grids[0][0], m_dmap[0][0], 1, 0);
    m_detJ.LocalCopy(detJ, 0, 0, 1, IntVect(0));
}

mlterrain::BCInfo
MLTerrainPoisson::bcInfo (int mglev) const
{
    mlterrain::BCInfo r;
    Box const& domain = m_geom[0][mglev].Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        r.lo[idim] = domain.smallEnd(idim);
        r.hi[idim] = domain.bigEnd(idim);
        auto f = [] (LinOpBCType t) -> int {
            if (t == LinOpBCType::Periodic) {
                return 0;
            } else if (t == LinOpBCType::Neumann) {
                return 1;
            } else if (t == LinOpBCType::Dirichlet) {
                return -1;
            } else {
                amrex::Abort("MLTerrainPoisson: only Periodic, Neumann and Dirichlet BCs are supported");
                return 0;
            }
        };
        r.bclo[idim] = f(m_lobc[0][idim]);
        r.bchi[idim] = f(m_hibc[0][idim]);
    }
    r.hidden = hiddenDirection();
    return r;
}

IntVect
MLTerrainPoisson::coarsenRatio (int mglev) const
{
    IntVect ratio = mg_coarsen_ratio_vec[mglev-1];
    if (hasHiddenDimension()) { ratio[hiddenDirection()] = 1; }
    return ratio;
}

void
MLTerrainPoisson::fillZPhysGhost (MultiFab& zp, int mglev) const
{
    Geometry const& geom = m_geom[0][mglev];
    // Sync shared nodes, or periodic ghost fills depend on copy order.
    zp.FillBoundaryAndSync(geom.periodicity());

    Box const ndomain = amrex::surroundingNodes(geom.Domain());
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) { continue; }
        int const dlo = ndomain.smallEnd(idim);
        int const dhi = ndomain.bigEnd(idim);
        if (dhi - dlo < 1) { continue; }
        for (MFIter mfi(zp); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            // Directions not yet filled are limited to the domain.
            Box gbx = mfi.fabbox();
            for (int jdim = idim+1; jdim < AMREX_SPACEDIM; ++jdim) {
                if (!geom.isPeriodic(jdim)) {
                    gbx.setSmall(jdim, std::max(gbx.smallEnd(jdim), ndomain.smallEnd(jdim)));
                    gbx.setBig  (jdim, std::min(gbx.bigEnd  (jdim), ndomain.bigEnd  (jdim)));
                }
            }
            auto const& a = zp.array(mfi);
            if (vbx.smallEnd(idim) == dlo) {
                Box b = gbx;
                b.setRange(idim, dlo-1);
                ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect iv(AMREX_D_DECL(i,j,k));
                    IntVect i0 = iv, i1 = iv;
                    i0[idim] = dlo;
                    i1[idim] = dlo+1;
                    a(iv) = Real(2.0)*a(i0) - a(i1);
                });
            }
            if (vbx.bigEnd(idim) == dhi) {
                Box b = gbx;
                b.setRange(idim, dhi+1);
                ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect iv(AMREX_D_DECL(i,j,k));
                    IntVect i0 = iv, i1 = iv;
                    i0[idim] = dhi;
                    i1[idim] = dhi-1;
                    a(iv) = Real(2.0)*a(i0) - a(i1);
                });
            }
        }
    }
}

void
MLTerrainPoisson::prepareForSolve ()
{
    BL_PROFILE("MLTerrainPoisson::prepareForSolve()");

    MLCellLinOpT<MultiFab>::prepareForSolve();

    bool singular = true;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_lobc[0][idim] == LinOpBCType::Dirichlet ||
            m_hibc[0][idim] == LinOpBCType::Dirichlet) {
            singular = false;
        }
    }
    m_is_singular[0] = singular;

    updateCoefs();
}

void
MLTerrainPoisson::update ()
{
    if (MLCellLinOpT<MultiFab>::needsUpdate()) { MLCellLinOpT<MultiFab>::update(); }
    if (m_needs_update) { updateCoefs(); }
}

void
MLTerrainPoisson::updateCoefs ()
{
    BL_PROFILE("MLTerrainPoisson::updateCoefs()");

    const int nmglevs = m_num_mg_levels[0];

    // Coarse heights are only needed for the metric terms.
    MultiFab zp_crse;
    MultiFab const* zp = &m_zphys;

    for (int mglev = 0; mglev < nmglevs; ++mglev) {
        if (mglev > 0) {
            IntVect const ratio = coarsenRatio(mglev);
            MultiFab zp_tmp(amrex::convert(m_grids[0][mglev], IntVect(1)), m_dmap[0][mglev], 1, 1);
            amrex::average_down_nodal(*zp, zp_tmp, ratio);
            fillZPhysGhost(zp_tmp, mglev);
            zp_crse = std::move(zp_tmp);
            zp = &zp_crse;
            amrex::average_down_faces(GetArrOfConstPtrs(m_area[mglev-1]),
                                      GetArrOfPtrs(m_area[mglev]), ratio, 0);
        }

        auto const bci = bcInfo(mglev);
        auto const dxinv = m_geom[0][mglev].InvCellSizeArray();
        auto& lu = m_column_lu[mglev];
        // A single whole column of a singular problem is singular itself.
        Box const& domain = m_geom[0][mglev].Domain();
        bool const singular_column = m_is_singular[0] && bci.bclo[2] != 0
            && domain.length(0) == 1 && domain.length(1) == 1;
        bool const spike = m_spike[mglev].isDefined();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_grids[0][mglev], m_dmap[0][mglev]); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            auto const& zpa = zp->const_array(mfi);
            auto const& axa = m_area[mglev][0].const_array(mfi);
            auto const& aya = m_area[mglev][1].const_array(mfi);
            auto const& aza = m_area[mglev][2].const_array(mfi);
            auto const& rxa = m_rx[mglev].array(mfi);
            auto const& rya = m_ry[mglev].array(mfi);
            auto const& zfa = m_zf[mglev].array(mfi);
            ParallelFor(amrex::convert(vbx,IntVect(1,0,1)), amrex::convert(vbx,IntVect(0,1,1)),
                        amrex::convert(vbx,IntVect(0,0,1)),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rxa(i,j,k) = mlterrain::metric_rx(i, j, k, zpa, dxinv[0]);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rya(i,j,k) = mlterrain::metric_ry(i, j, k, zpa, dxinv[1]);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                mlterrain::metric_zf(i, j, k, zfa, zpa, dxinv[0], dxinv[1]);
            });
            if (!lu.isDefined()) { continue; }
            auto const& lua = lu.array(mfi);
            mlterrain::Metrics const met{.rx = rxa, .ry = rya, .zf = zfa};
            int const klo = vbx.smallEnd(2);
            int const khi = vbx.bigEnd(2);
            // Spikes keep the coupling to the neighboring blocks.
            int const clo = spike ? domain.smallEnd(2) : klo;
            int const chi = spike ? domain.bigEnd(2) : khi;
            bool const pin_top = singular_column && khi == domain.bigEnd(2);
            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                mlterrain::column_coefs(i, j, k, clo, chi, lua, met, axa, aya, aza,
                                        dxinv[0], dxinv[1], dxinv[2], bci);
            });
            if (spike) {
                auto const& spka = m_spike[mglev].array(mfi);
                ParallelFor(amrex::makeSlab(vbx,2,klo), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Real const as = lua(i,j,klo,0);
                    Real const ce = lua(i,j,khi,2);
                    mlterrain::factor_column(i, j, klo, khi, lua, pin_top);
                    // W = D^{-1} (as e_klo) and V = D^{-1} (ce e_khi).
                    Real w = as * lua(i,j,klo,1);
                    spka(i,j,klo,0) = w;
                    for (int k = klo+1; k <= khi; ++k) {
                        w = -(lua(i,j,k,0)*w) * lua(i,j,k,1);
                        spka(i,j,k,0) = w;
                    }
                    for (int k = khi-1; k >= klo; --k) {
                        w = spka(i,j,k,0) - lua(i,j,k+1,2) * w;
                        spka(i,j,k,0) = w;
                    }
                    Real v = ce * lua(i,j,khi,1);
                    spka(i,j,khi,1) = v;
                    for (int k = khi-1; k >= klo; --k) {
                        v = -(lua(i,j,k+1,2) * v);
                        spka(i,j,k,1) = v;
                    }
                });
            } else {
                ParallelFor(amrex::makeSlab(vbx,2,klo), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    mlterrain::factor_column(i, j, klo, khi, lua, pin_top);
                });
            }
        }

        if (spike) {
            // Gather W and V at the ends of all blocks of each column.
            MultiFab emf(m_spike_buf[mglev].boxArray(), m_dmap[0][mglev], 4, 0);
            auto const& ea = emf.arrays();
            auto const& spk = m_spike[mglev].const_arrays();
            ParallelFor(emf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int q)
            {
                auto const& w = spk[b];
                int const klo = w.begin[2];
                int const khi = w.end[2]-1;
                ea[b](i,j,q,0) = w(i,j,klo,0);
                ea[b](i,j,q,1) = w(i,j,khi,0);
                ea[b](i,j,q,2) = w(i,j,klo,1);
                ea[b](i,j,q,3) = w(i,j,khi,1);
            });
            m_spike_end[mglev].ParallelCopy(emf);
        }

        if (m_col_lu[mglev].isDefined()) {
            // Factor whole columns on the column layout.
            auto& clu = m_col_lu[mglev];
            BoxArray const& cba = clu.boxArray();
            DistributionMapping const& cdm = clu.DistributionMap();
            MultiFab crx(amrex::convert(cba,IntVect(1,0,1)), cdm, 1, 0);
            MultiFab cry(amrex::convert(cba,IntVect(0,1,1)), cdm, 1, 0);
            MultiFab czf(amrex::convert(cba,IntVect(0,0,1)), cdm, 3, 0);
            crx.ParallelCopy(m_rx[mglev]);
            cry.ParallelCopy(m_ry[mglev]);
            czf.ParallelCopy(m_zf[mglev]);
            Array<MultiFab,AMREX_SPACEDIM> carea;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                carea[idim].define(amrex::convert(cba, IntVect::TheDimensionVector(idim)), cdm, 1, 0);
                carea[idim].ParallelCopy(m_area[mglev][idim]);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(clu); mfi.isValid(); ++mfi) {
                Box const& vbx = mfi.validbox();
                auto const& lua = clu.array(mfi);
                mlterrain::Metrics const met{.rx = crx.const_array(mfi),
                                             .ry = cry.const_array(mfi),
                                             .zf = czf.const_array(mfi)};
                auto const& axa = carea[0].const_array(mfi);
                auto const& aya = carea[1].const_array(mfi);
                auto const& aza = carea[2].const_array(mfi);
                int const klo = vbx.smallEnd(2);
                int const khi = vbx.bigEnd(2);
                ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    mlterrain::column_coefs(i, j, k, klo, khi, lua, met, axa, aya, aza,
                                            dxinv[0], dxinv[1], dxinv[2], bci);
                });
                ParallelFor(amrex::makeSlab(vbx,2,klo), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    mlterrain::factor_column(i, j, klo, khi, lua, singular_column);
                });
            }
        }
    }

    m_needs_update = false;
}

bool
MLTerrainPoisson::scaleRHS (int amrlev, MultiFab* rhs) const
{
    AMREX_ASSERT(amrlev == 0);
    amrex::ignore_unused(amrlev);
    if (!m_detJ.isDefined()) { return false; }
    if (rhs) {
        MultiFab::Multiply(*rhs, m_detJ, 0, 0, 1, 0);
    }
    return true;
}

void
MLTerrainPoisson::applyBC (int amrlev, int mglev, MultiFab& in, BCMode /*bc_mode*/,
                           StateMode /*s_mode*/, const MLMGBndryT<MultiFab>* /*bndry*/,
                           bool skip_fillboundary) const
{
    BL_PROFILE("MLTerrainPoisson::applyBC()");
    AMREX_ASSERT(amrlev == 0);
    amrex::ignore_unused(amrlev);

    IntVect ng(1);
    if (hasHiddenDimension()) { ng[hiddenDirection()] = 0; }
    AMREX_ASSERT(in.nGrowVect().allGE(ng));

    Geometry const& geom = m_geom[0][mglev];
    if (!skip_fillboundary) {
        in.FillBoundary(0, 1, ng, geom.periodicity());
    }

    auto const bci = bcInfo(mglev);
    auto& tags = m_bc_tags[mglev];
    if (!tags.is_defined()) {
        // Disjoint ghost regions outside non-periodic domain faces.
        Vector<BCTag> vtags;
        for (MFIter mfi(in); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (ng[idim] == 0) { continue; }
                Box gbx = amrex::grow(vbx, ng);
                for (int jdim = idim+1; jdim < AMREX_SPACEDIM; ++jdim) {
                    if (bci.bclo[jdim] != 0) { gbx.setSmall(jdim, std::max(gbx.smallEnd(jdim), bci.lo[jdim])); }
                    if (bci.bchi[jdim] != 0) { gbx.setBig  (jdim, std::min(gbx.bigEnd  (jdim), bci.hi[jdim])); }
                }
                if (bci.bclo[idim] != 0 && vbx.smallEnd(idim) == bci.lo[idim]) {
                    vtags.push_back(BCTag{.bx = Box(gbx).setRange(idim, bci.lo[idim]-1),
                                          .local_index = mfi.LocalIndex()});
                }
                if (bci.bchi[idim] != 0 && vbx.bigEnd(idim) == bci.hi[idim]) {
                    vtags.push_back(BCTag{.bx = Box(gbx).setRange(idim, bci.hi[idim]+1),
                                          .local_index = mfi.LocalIndex()});
                }
            }
        }
        tags.define(vtags);
    }

    if (tags.ntags > 0) {
        auto const& ma = in.arrays();
        // Reflect across every non-periodic domain face the cell is outside of.
        ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, BCTag const& tag)
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            IntVect im = iv;
            int sign = 1;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (iv[d] < bci.lo[d] && bci.bclo[d] != 0) {
                    im[d] = 2*bci.lo[d] - 1 - iv[d];
                    sign *= bci.bclo[d];
                } else if (iv[d] > bci.hi[d] && bci.bchi[d] != 0) {
                    im[d] = 2*bci.hi[d] + 1 - iv[d];
                    sign *= bci.bchi[d];
                }
            }
            auto const& a = ma[tag.local_index];
            a(iv) = Real(sign) * a(im);
        });
    }
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}

void
MLTerrainPoisson::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLTerrainPoisson::Fapply()");
    AMREX_ASSERT(amrlev == 0);
    amrex::ignore_unused(amrlev);

    auto const bci = bcInfo(mglev);
    auto const dxinv = m_geom[0][mglev].InvCellSizeArray();
    auto const& y = out.arrays();
    auto const& x = in.const_arrays();
    auto const& rx = m_rx[mglev].const_arrays();
    auto const& ry = m_ry[mglev].const_arrays();
    auto const& zf = m_zf[mglev].const_arrays();
    auto const& ax = m_area[mglev][0].const_arrays();
    auto const& ay = m_area[mglev][1].const_arrays();
    auto const& az = m_area[mglev][2].const_arrays();
    ParallelFor(out, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        mlterrain::Metrics const met{.rx = rx[b], .ry = ry[b], .zf = zf[b]};
        y[b](i,j,k) = mlterrain::adotx(i, j, k, x[b], met, ax[b], ay[b], az[b],
                                       dxinv[0], dxinv[1], dxinv[2], bci);
    });
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}

void
MLTerrainPoisson::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLTerrainPoisson::interpolation()");
    AMREX_ASSERT(amrlev == 0);
    amrex::ignore_unused(amrlev);

    // Piecewise constant for ratio (2,2), linear in x and y otherwise.
    IntVect const ratio = coarsenRatio(fmglev+1);
    if (ratio[0] == 2 && ratio[1] == 2) {
        MLCellLinOpT<MultiFab>::interpolation(amrlev, fmglev, fine, crse);
        return;
    }
    IntVect ng(0);
    for (int d = 0; d < 2; ++d) { ng[d] = (ratio[d] > 1) ? 1 : 0; }
    MultiFab ct(amrex::coarsen(fine.boxArray(), ratio), fine.DistributionMap(), 1, ng,
                MFInfo().SetArena(The_Async_Arena()));
    ct.ParallelCopy(crse, 0, 0, 1, IntVect(0), ng, m_geom[0][fmglev+1].periodicity());

    auto const bci = bcInfo(fmglev+1);
    auto const& fa = fine.arrays();
    auto const& ca = ct.const_arrays();
    ParallelFor(fine, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect const iv(i,j,k);
        IntVect ic = iv, in = iv;
        Real w[2], s[2];
        for (int d = 0; d < 2; ++d) {
            ic[d] = amrex::coarsen(iv[d], ratio[d]);
            Real const t = (Real(iv[d] - ratio[d]*ic[d]) + Real(0.5)) / Real(ratio[d]) - Real(0.5);
            w[d] = std::abs(t);
            in[d] = ic[d] + ((t < Real(0.0)) ? -1 : 1);
            s[d] = Real(1.0);
            // Reflect at non-periodic domain faces: even (Neumann), odd (Dirichlet).
            if (in[d] < bci.lo[d] && bci.bclo[d] != 0) {
                in[d] = ic[d];
                s[d] = Real(bci.bclo[d]);
            } else if (in[d] > bci.hi[d] && bci.bchi[d] != 0) {
                in[d] = ic[d];
                s[d] = Real(bci.bchi[d]);
            }
        }
        auto const& c = ca[b];
        Real v = (Real(1.0)-w[0])*(Real(1.0)-w[1])*c(ic);
        if (w[0] > Real(0.0)) { v += w[0]*(Real(1.0)-w[1])*s[0]*c(in[0],ic[1],k); }
        if (w[1] > Real(0.0)) { v += (Real(1.0)-w[0])*w[1]*s[1]*c(ic[0],in[1],k); }
        if (w[0] > Real(0.0) && w[1] > Real(0.0)) {
            v += w[0]*w[1]*s[0]*s[1]*c(in[0],in[1],k);
        }
        fa[b](i,j,k) += v;
    });
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}

void
MLTerrainPoisson::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                           int redblack) const
{
    BL_PROFILE("MLTerrainPoisson::Fsmooth()");
    AMREX_ASSERT(amrlev == 0);
    amrex::ignore_unused(amrlev);

    auto const bci = bcInfo(mglev);
    auto const dxinv = m_geom[0][mglev].InvCellSizeArray();

    MFItInfo mfi_info;
    // Tiles must span the box in z.
    mfi_info.EnableTiling(IntVect(AMREX_D_DECL(1024000,8,1024000))).SetDynamic(true);

    // GPUs and boxes split in z: the residual first, then the line solves.
    // Otherwise the fused CPU path below.
    bool const spike = m_spike[mglev].isDefined();
    bool const split = m_col_lu[mglev].isDefined();
    bool const gpu = Gpu::inLaunchRegion();
    if (gpu || spike || split) {
        auto& res = m_smooth_res[mglev];
        if (!res.isDefined()) {
            res.define(m_grids[0][mglev], m_dmap[0][mglev], 1, 0);
        }
        auto const& ta = res.arrays();
        auto const& sa = sol.arrays();

        if (gpu) {
            auto const& ra = rhs.const_arrays();
            auto const& rx = m_rx[mglev].const_arrays();
            auto const& ry = m_ry[mglev].const_arrays();
            auto const& zf = m_zf[mglev].const_arrays();
            auto const& ax = m_area[mglev][0].const_arrays();
            auto const& ay = m_area[mglev][1].const_arrays();
            auto const& az = m_area[mglev][2].const_arrays();
            // Columns with (i+j+redblack) even, as two strided lattices.
            IntVect const stride(2,2,1);
            IntVect const off0(redblack,0,0);
            IntVect const off1(1-redblack,1,0);
            for (auto const& off : {off0, off1}) {
                ParallelForStrided(sol, stride, off,
                [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                {
                    mlterrain::Metrics const met{.rx = rx[b], .ry = ry[b], .zf = zf[b]};
                    ta[b](i,j,k) = ra[b](i,j,k)
                        - mlterrain::adotx(i, j, k, sa[b], met, ax[b], ay[b], az[b],
                                           dxinv[0], dxinv[1], dxinv[2], bci);
                });
            }
        } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(sol, mfi_info); mfi.isValid(); ++mfi) {
                Box const& tbx = mfi.tilebox();
                auto const& tfab = res.array(mfi);
                auto const& sfab = sol.const_array(mfi);
                auto const& rfab = rhs.const_array(mfi);
                mlterrain::Metrics const met{.rx = m_rx[mglev].const_array(mfi),
                                             .ry = m_ry[mglev].const_array(mfi),
                                             .zf = m_zf[mglev].const_array(mfi)};
                auto const& axa = m_area[mglev][0].const_array(mfi);
                auto const& aya = m_area[mglev][1].const_array(mfi);
                auto const& aza = m_area[mglev][2].const_array(mfi);
                auto const lo = amrex::lbound(tbx);
                auto const hi = amrex::ubound(tbx);
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                    AMREX_PRAGMA_SIMD
                    for (int i = ilo; i <= hi.x; i += 2) {
                        tfab(i,j,k) = rfab(i,j,k) - mlterrain::adotx(i, j, k, sfab, met, axa, aya, aza,
                                                                     dxinv[0], dxinv[1], dxinv[2], bci);
                    }
                }}
            }
        }

        if (spike) {
            auto& buf = m_spike_buf[mglev];
            auto& sum = m_spike_sum[mglev];
            int const np = static_cast<int>(m_zcut[mglev].size()) - 1;
            int const nr = 2*(np-1);
            auto const& ya = buf.arrays();
            // Block solves y = D^{-1} r of the active color, in place, with
            // y at the block start and end in the box's x-y slab (k = block).
            if (gpu) {
                auto const& lua = m_column_lu[mglev].const_arrays();
                ParallelFor(buf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int q)
                {
                    ya[b](i,j,q,0) = Real(0.0);
                    ya[b](i,j,q,1) = Real(0.0);
                    if (((i+j+redblack) & 1) == 0) {
                        auto const& lu = lua[b];
                        auto const& t = ta[b];
                        int const klo = lu.begin[2];
                        int const khi = lu.end[2]-1;
                        Real u = Real(0.0);
                        for (int k = klo; k <= khi; ++k) {
                            u = (t(i,j,k) - lu(i,j,k,0)*u) * lu(i,j,k,1);
                            t(i,j,k) = u;
                        }
                        for (int k = khi-1; k >= klo; --k) {
                            u = t(i,j,k) - lu(i,j,k+1,2) * u;
                            t(i,j,k) = u;
                        }
                        ya[b](i,j,q,0) = t(i,j,klo);
                        ya[b](i,j,q,1) = t(i,j,khi);
                    }
                });
            } else {
                buf.setVal(Real(0.0));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(res, mfi_info); mfi.isValid(); ++mfi) {
                    Box const& tbx = mfi.tilebox();
                    auto const& tfab = res.array(mfi);
                    auto const& yfab = buf.array(mfi);
                    auto const& lua = m_column_lu[mglev].const_array(mfi);
                    auto const lo = amrex::lbound(tbx);
                    auto const hi = amrex::ubound(tbx);
                    int const q = yfab.begin[2];
                    for (int k = lo.z; k <= hi.z; ++k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                        if (k == lo.z) {
                            AMREX_PRAGMA_SIMD
                            for (int i = ilo; i <= hi.x; i += 2) {
                                tfab(i,j,k) = tfab(i,j,k) * lua(i,j,k,1);
                            }
                        } else {
                            AMREX_PRAGMA_SIMD
                            for (int i = ilo; i <= hi.x; i += 2) {
                                tfab(i,j,k) = (tfab(i,j,k) - lua(i,j,k,0)*tfab(i,j,k-1)) * lua(i,j,k,1);
                            }
                        }
                    }}
                    for (int k = hi.z-1; k >= lo.z; --k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                        AMREX_PRAGMA_SIMD
                        for (int i = ilo; i <= hi.x; i += 2) {
                            tfab(i,j,k) -= lua(i,j,k+1,2) * tfab(i,j,k+1);
                        }
                    }}
                    for (int j = lo.y; j <= hi.y; ++j) {
                        int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                        for (int i = ilo; i <= hi.x; i += 2) {
                            yfab(i,j,q,0) = tfab(i,j,lo.z);
                            yfab(i,j,q,1) = tfab(i,j,hi.z);
                        }
                    }
                }
            }

            // Interface values from the reduced system of each column:
            // l_m + W_m(e) l_{m-1} + V_m(e) f_{m+1} = y_m(e) and
            // f_{m+1} + W_{m+1}(s) l_m + V_{m+1}(s) f_{m+2} = y_{m+1}(s),
            // with l = x(block end), f = x(block start), unknowns ordered
            // (l_0, f_1, l_1, f_2, ...).  Block q gets (l_{q-1}, f_{q+1}).
            sum.ParallelCopy(buf);
            auto const& ea = m_spike_end[mglev].const_arrays();
            auto const& ua = sum.arrays();
            auto const& la = m_spike_lu[mglev].arrays();
            ParallelFor(sum, [=] AMREX_GPU_DEVICE (int b, int i, int j, int kq)
            {
                if (kq != 0 || ((i+j+redblack) & 1) != 0) { return; }
                auto const& e = ea[b];
                auto const& u = ua[b];
                auto const& lr = la[b];
                // Forward elimination with a window of rows c, c+1, c+2
                // (band columns -2..+2); row c's diagonal, upper band and
                // rhs are kept in lr for the back substitution.
                Real w[3][5] = {};
                Real wr[3] = {};
                for (int c = -2; c < nr; ++c) {
                    int const row = c+2;
                    if (row < nr) {
                        int const m = row/2;
                        if (row % 2 == 0) {
                            w[2][0] = (m > 0) ? e(i,j,m,1) : Real(0.0);
                            w[2][1] = Real(0.0);
                            w[2][2] = Real(1.0);
                            w[2][3] = e(i,j,m,3);
                            w[2][4] = Real(0.0);
                            wr[2] = u(i,j,m,1);
                        } else {
                            w[2][0] = Real(0.0);
                            w[2][1] = e(i,j,m+1,0);
                            w[2][2] = Real(1.0);
                            w[2][3] = Real(0.0);
                            w[2][4] = (row+2 < nr) ? e(i,j,m+1,2) : Real(0.0);
                            wr[2] = u(i,j,m+1,0);
                        }
                    }
                    if (c >= 0) {
                        for (int t = 1; t <= 2 && c+t < nr; ++t) {
                            Real const fac = w[t][2-t] / w[0][2];
                            for (int cc = c; cc <= c+2 && cc < nr; ++cc) {
                                w[t][cc-c-t+2] -= fac * w[0][cc-c+2];
                            }
                            wr[t] -= fac * wr[0];
                        }
                        int const n0 = 4*(c%2);
                        lr(i,j,c/2,n0  ) = w[0][2];
                        lr(i,j,c/2,n0+1) = w[0][3];
                        lr(i,j,c/2,n0+2) = w[0][4];
                        lr(i,j,c/2,n0+3) = wr[0];
                    }
                    for (int n = 0; n < 5; ++n) {
                        w[0][n] = w[1][n];
                        w[1][n] = w[2][n];
                    }
                    wr[0] = wr[1];
                    wr[1] = wr[2];
                }
                Real x1 = Real(0.0); // x[c+1]
                Real x2 = Real(0.0); // x[c+2]
                for (int c = nr-1; c >= 0; --c) {
                    int const n0 = 4*(c%2);
                    Real x = lr(i,j,c/2,n0+3);
                    if (c+1 < nr) { x -= lr(i,j,c/2,n0+1) * x1; }
                    if (c+2 < nr) { x -= lr(i,j,c/2,n0+2) * x2; }
                    x /= lr(i,j,c/2,n0);
                    if (c % 2 == 0) {
                        u(i,j,c/2+1,0) = x; // l_m for block m+1
                    } else {
                        u(i,j,c/2,1) = x;   // f_{m+1} for block m
                    }
                    x2 = x1;
                    x1 = x;
                }
                u(i,j,0,0) = Real(0.0);
                u(i,j,np-1,1) = Real(0.0);
            });
            buf.ParallelCopy(sum);

            // x = y - W l_{q-1} - V f_{q+1}, added to sol.
            auto const& spk = m_spike[mglev].const_arrays();
            if (gpu) {
                ParallelFor(res, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                {
                    if (((i+j+redblack) & 1) == 0) {
                        auto const& y = ya[b];
                        int const q = y.begin[2];
                        sa[b](i,j,k) += ta[b](i,j,k) - spk[b](i,j,k,0)*y(i,j,q,0)
                                                     - spk[b](i,j,k,1)*y(i,j,q,1);
                    }
                });
            } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(res, mfi_info); mfi.isValid(); ++mfi) {
                    Box const& tbx = mfi.tilebox();
                    auto const& tfab = res.const_array(mfi);
                    auto const& sfab = sol.array(mfi);
                    auto const& yfab = buf.const_array(mfi);
                    auto const& wv = m_spike[mglev].const_array(mfi);
                    auto const lo = amrex::lbound(tbx);
                    auto const hi = amrex::ubound(tbx);
                    int const q = yfab.begin[2];
                    for (int k = lo.z; k <= hi.z; ++k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                        AMREX_PRAGMA_SIMD
                        for (int i = ilo; i <= hi.x; i += 2) {
                            sfab(i,j,k) += tfab(i,j,k) - wv(i,j,k,0)*yfab(i,j,q,0)
                                                       - wv(i,j,k,1)*yfab(i,j,q,1);
                        }
                    }}
                }
            }
        } else if (split) {
            auto& cres = m_col_res[mglev];
            auto& ccor = m_col_cor[mglev];
            cres.ParallelCopy(res);
            ccor.setVal(Real(0.0));
            if (gpu) {
                auto& slab = m_column_slab[mglev];
                if (!slab.isDefined()) {
                    BoxList bl;
                    for (int i = 0, N = static_cast<int>(ccor.size()); i < N; ++i) {
                        Box const& cb = ccor.boxArray()[i];
                        bl.push_back(amrex::makeSlab(cb, 2, cb.smallEnd(2)));
                    }
                    slab.define(BoxArray(std::move(bl)), ccor.DistributionMap(), 1, 0,
                                MFInfo().SetAlloc(false));
                }
                auto const& cta = cres.arrays();
                auto const& cca = ccor.arrays();
                auto const& clua = m_col_lu[mglev].const_arrays();
                ParallelFor(slab, [=] AMREX_GPU_DEVICE (int b, int i, int j, int)
                {
                    if (((i+j+redblack) & 1) == 0) {
                        auto const& lu = clua[b];
                        mlterrain::solve_column(i, j, lu.begin[2], lu.end[2]-1, cca[b], cta[b], lu);
                    }
                });
            } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(ccor, mfi_info); mfi.isValid(); ++mfi) {
                    Box const& tbx = mfi.tilebox();
                    auto const& tfab = cres.array(mfi);
                    auto const& cfab = ccor.array(mfi);
                    auto const& lua = m_col_lu[mglev].const_array(mfi);
                    auto const lo = amrex::lbound(tbx);
                    auto const hi = amrex::ubound(tbx);
                    for (int k = lo.z; k <= hi.z; ++k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                        if (k == lo.z) {
                            AMREX_PRAGMA_SIMD
                            for (int i = ilo; i <= hi.x; i += 2) {
                                tfab(i,j,k) = tfab(i,j,k) * lua(i,j,k,1);
                            }
                        } else {
                            AMREX_PRAGMA_SIMD
                            for (int i = ilo; i <= hi.x; i += 2) {
                                tfab(i,j,k) = (tfab(i,j,k) - lua(i,j,k,0)*tfab(i,j,k-1)) * lua(i,j,k,1);
                            }
                        }
                    }}
                    for (int k = hi.z; k >= lo.z; --k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                        if (k < hi.z) {
                            AMREX_PRAGMA_SIMD
                            for (int i = ilo; i <= hi.x; i += 2) {
                                tfab(i,j,k) -= lua(i,j,k+1,2) * tfab(i,j,k+1);
                            }
                        }
                        AMREX_PRAGMA_SIMD
                        for (int i = ilo; i <= hi.x; i += 2) {
                            cfab(i,j,k) = tfab(i,j,k);
                        }
                    }}
                }
            }
            // The other color gets zero.
            sol.ParallelAdd(ccor);
        } else {
            // One launch for the latency-bound column solves.
            auto& slab = m_column_slab[mglev];
            if (!slab.isDefined()) {
                BoxList bl;
                for (int i = 0, N = static_cast<int>(sol.size()); i < N; ++i) {
                    Box const& vb = m_grids[0][mglev][i];
                    bl.push_back(amrex::makeSlab(vb, 2, vb.smallEnd(2)));
                }
                slab.define(BoxArray(std::move(bl)), m_dmap[0][mglev], 1, 0,
                            MFInfo().SetAlloc(false));
            }
            auto const& lua = m_column_lu[mglev].const_arrays();
            ParallelFor(slab, [=] AMREX_GPU_DEVICE (int b, int i, int j, int)
            {
                if (((i+j+redblack) & 1) == 0) {
                    auto const& lu = lua[b];
                    mlterrain::solve_column(i, j, lu.begin[2], lu.end[2]-1, sa[b], ta[b], lu);
                }
            });
        }
        return;
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox tmp;
        for (MFIter mfi(sol, mfi_info); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.tilebox();
            tmp.resize(tbx, 1, The_Async_Arena());
            auto const& ta = tmp.array();
            auto const& sa = sol.array(mfi);
            auto const& ra = rhs.const_array(mfi);
            auto const& lua = m_column_lu[mglev].const_array(mfi);
            mlterrain::Metrics const met{.rx = m_rx[mglev].const_array(mfi),
                                         .ry = m_ry[mglev].const_array(mfi),
                                         .zf = m_zf[mglev].const_array(mfi)};
            auto const& axa = m_area[mglev][0].const_array(mfi);
            auto const& aya = m_area[mglev][1].const_array(mfi);
            auto const& aza = m_area[mglev][2].const_array(mfi);
            // Plane by plane so that the loops over i vectorize.  The
            // first pass computes the residual and eliminates forward,
            // and the second substitutes backward and updates sol.
            auto const lo = amrex::lbound(tbx);
            auto const hi = amrex::ubound(tbx);
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                if (k == lo.z) {
                    AMREX_PRAGMA_SIMD
                    for (int i = ilo; i <= hi.x; i += 2) {
                        Real const r = ra(i,j,k) - mlterrain::adotx(i, j, k, sa, met, axa, aya, aza,
                                                                    dxinv[0], dxinv[1], dxinv[2], bci);
                        ta(i,j,k) = r * lua(i,j,k,1);
                    }
                } else {
                    AMREX_PRAGMA_SIMD
                    for (int i = ilo; i <= hi.x; i += 2) {
                        Real const r = ra(i,j,k) - mlterrain::adotx(i, j, k, sa, met, axa, aya, aza,
                                                                    dxinv[0], dxinv[1], dxinv[2], bci);
                        ta(i,j,k) = (r - lua(i,j,k,0)*ta(i,j,k-1)) * lua(i,j,k,1);
                    }
                }
            }}
            for (int k = hi.z; k >= lo.z; --k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                int const ilo = lo.x + ((lo.x+j+redblack) & 1);
                if (k < hi.z) {
                    AMREX_PRAGMA_SIMD
                    for (int i = ilo; i <= hi.x; i += 2) {
                        ta(i,j,k) -= lua(i,j,k+1,2) * ta(i,j,k+1);
                    }
                }
                AMREX_PRAGMA_SIMD
                for (int i = ilo; i <= hi.x; i += 2) {
                    sa(i,j,k) += ta(i,j,k);
                }
            }}
        }
    }
}

void
MLTerrainPoisson::FFlux (int amrlev, const MFIter& mfi, const Array<FAB*,AMREX_SPACEDIM>& flux,
                         const FAB& sol, Location /*loc*/, int /*face_only*/) const
{
    AMREX_ASSERT(amrlev == 0);
    amrex::ignore_unused(amrlev);

    const int mglev = 0;
    auto const bci = bcInfo(mglev);
    auto const dxinv = m_geom[0][mglev].InvCellSizeArray();
    Box const& tbx = mfi.tilebox();
    auto const& x = sol.const_array();
    mlterrain::Metrics const met{.rx = m_rx[mglev].const_array(mfi),
                                 .ry = m_ry[mglev].const_array(mfi),
                                 .zf = m_zf[mglev].const_array(mfi)};
    auto const& fx = flux[0]->array();
    auto const& fy = flux[1]->array();
    auto const& fz = flux[2]->array();
    ParallelFor(amrex::surroundingNodes(tbx,0), amrex::surroundingNodes(tbx,1),
                amrex::surroundingNodes(tbx,2),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        fx(i,j,k) = mlterrain::flux_x(i, j, k, x, met, dxinv[0], bci);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        fy(i,j,k) = mlterrain::flux_y(i, j, k, x, met, dxinv[1], bci);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        fz(i,j,k) = mlterrain::flux_z(i, j, k, x, met, dxinv[0], dxinv[1], bci);
    });
}

void
MLTerrainPoisson::compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                            MultiFab& sol, Location loc) const
{
    compFlux(amrlev, grad, sol, loc);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        grad[idim]->mult(Real(-1.0));
    }
}

}
