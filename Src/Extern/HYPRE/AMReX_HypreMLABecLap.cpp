#include <AMReX_HypreMLABecLap.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_HypreMLABecLap_K.H>

#include <functional>
#include <numeric>
#include <tuple>

namespace amrex {

HypreMLABecLap::HypreMLABecLap (Vector<Geometry> a_geom,
                                Vector<BoxArray> a_grids,
                                Vector<DistributionMapping> a_dmap,
                                HypreSolverID a_hypre_solver_id,
                                std::string a_parmparse_prefix)
    : m_geom(std::move(a_geom)),
      m_grids(std::move(a_grids)),
      m_dmap(std::move(a_dmap)),
      m_parmparse_prefix(std::move(a_parmparse_prefix)),
      m_nlevels(int(m_grids.size())),
      m_comm(ParallelContext::CommunicatorSub()),
      m_hypre_solver_id(a_hypre_solver_id),
      m_hypre_object_type((a_hypre_solver_id == HypreSolverID::BoomerAMG) ?
                          HYPRE_PARCSR : HYPRE_SSTRUCT)
{
    BL_PROFILE("HypreMLABecLap::HypreMLABecLap");

#ifndef AMREX_FEATURE_HYPRE_SSAMG
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_hypre_solver_id == HypreSolverID::BoomerAMG,
                                     "HypreMLABecLap only supports BoomerAMG ifndef AMREX_FEATURE_HYPRE_SSAMG");
#endif

    m_ref_ratio.resize(m_nlevels-1);
    for (int ilev = 0; ilev < m_nlevels-1; ++ilev) {
        m_ref_ratio[ilev] = m_geom[ilev+1].Domain().length()
            /               m_geom[ilev  ].Domain().length();
        AMREX_ASSERT(m_geom[ilev+1].Domain() == amrex::refine(m_geom[ilev].Domain(),
                                                              m_ref_ratio[ilev])
                     && m_ref_ratio[ilev].allLE(4));
    }

    m_bndry.resize(m_nlevels);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        int const ncomp = 1;
        m_bndry[ilev] = std::make_unique<MLMGBndry>(m_grids[ilev],
                                                    m_dmap[ilev],
                                                    ncomp,
                                                    m_geom[ilev]);
    }

    m_bndry_rhs.resize(m_nlevels);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        int const ncomp = 1;
        m_bndry_rhs[ilev]= std::make_unique<BndryRegister>(m_grids[ilev],
                                                           m_dmap[ilev],
                                                           1, 0, 0, ncomp);
    }

    m_fine_masks.resize(m_nlevels-1);
    m_crse_masks.resize(m_nlevels-1);
    for (int ilev = 0; ilev < m_nlevels-1; ++ilev) {
        m_fine_masks[ilev] = amrex::makeFineMask(m_grids[ilev], m_dmap[ilev], IntVect(1),
                                                 m_grids[ilev+1], m_ref_ratio[ilev],
                                                 m_geom[ilev].periodicity(),
                                                 0, 1);
        m_crse_masks[ilev].define(m_grids[ilev], m_dmap[ilev], 1, 1);
        m_crse_masks[ilev].BuildMask(m_geom[ilev].Domain(),
                                     m_geom[ilev].periodicity(),
                                     1, 0, 0, 1);
    }

    m_c2f_offset_from.resize(m_nlevels-1);
    m_c2f_total_from.resize(m_nlevels-1);
    m_c2f_nentries.resize(m_nlevels-1);
    m_c2f_offset_to.resize(m_nlevels-1);
    m_c2f_total_to.resize(m_nlevels-1);
    for (int ilev = 0; ilev < m_nlevels-1; ++ilev) {
        m_c2f_offset_from[ilev].define(m_grids[ilev], m_dmap[ilev], 1, 0);
        m_c2f_total_from[ilev].define(m_grids[ilev], m_dmap[ilev]);
        m_c2f_nentries[ilev].define(m_grids[ilev], m_dmap[ilev], 1, 0);
        m_c2f_offset_to[ilev].define(m_grids[ilev], m_dmap[ilev], 1, 0);
        m_c2f_total_to[ilev].define(m_grids[ilev], m_dmap[ilev]);
    }

    m_offset_cf_bcoefs.resize(m_nlevels-1);
    m_cf_bcoefs.resize(m_nlevels-1);

    static_assert(std::is_same_v<HYPRE_Real,Real>,
                  "HYPRE_Real and amrex::Real must be the same type");

    HYPRE_SStructGridCreate(m_comm, AMREX_SPACEDIM, m_nlevels, &m_ss_grid);

    constexpr HYPRE_Int nvars = 1;
    constexpr HYPRE_Int ivar = 0;

    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        // Which hypre solver has the limitation of power of 2 restrictions
        // for periodic domains?
        if (m_geom[ilev].isAnyPeriodic()) {
            Array<HYPRE_Int,AMREX_SPACEDIM> periodic;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                periodic[idim] = m_geom[ilev].isPeriodic(idim)
                    ? m_geom[ilev].Domain().length(idim) : 0;
            }
            HYPRE_SStructGridSetPeriodic(m_ss_grid, ilev, periodic.data());
        }

        AMREX_ASSERT(m_grids[ilev].ixType().cellCentered());

        for (MFIter mfi(m_grids[ilev], m_dmap[ilev], MFItInfo().DisableDeviceSync());
             mfi.isValid(); ++mfi)
        {
            Box const& b = mfi.validbox();
            Array<HYPRE_Int,AMREX_SPACEDIM> lo{AMREX_D_DECL(b.smallEnd(0),
                                                            b.smallEnd(1),
                                                            b.smallEnd(2))};
            Array<HYPRE_Int,AMREX_SPACEDIM> hi{AMREX_D_DECL(b.bigEnd(0),
                                                            b.bigEnd(1),
                                                            b.bigEnd(2))};
            HYPRE_SStructGridSetExtents(m_ss_grid, ilev, lo.data(), hi.data());
        }

        auto vartype = HYPRE_SSTRUCT_VARIABLE_CELL;
        HYPRE_SStructVariable vars[nvars] = {vartype};
        HYPRE_SStructGridSetVariables(m_ss_grid, ilev, nvars, vars);
    }

    HYPRE_SStructGridAssemble(m_ss_grid);

#if (AMREX_SPACEDIM == 2)
    HYPRE_Int cross_stencil_offset[5][2] = {{ 0,  0},
                                            {-1,  0},
                                            { 1,  0},
                                            { 0, -1},
                                            { 0,  1}};
#elif (AMREX_SPACEDIM == 3)
    HYPRE_Int cross_stencil_offset[7][3] = {{ 0,  0,  0},
                                            {-1,  0,  0},
                                            { 1,  0,  0},
                                            { 0, -1,  0},
                                            { 0,  1,  0},
                                            { 0,  0, -1},
                                            { 0,  0,  1}};
#endif

    HYPRE_SStructStencilCreate(AMREX_SPACEDIM, 2*AMREX_SPACEDIM+1, &m_ss_stencil);
    for (HYPRE_Int i = 0; i < 2*AMREX_SPACEDIM+1; ++i) {
        HYPRE_SStructStencilSetEntry(m_ss_stencil, i, cross_stencil_offset[i], ivar);
    }

    HYPRE_SStructGraphCreate(m_comm, m_ss_grid, &m_ss_graph);
    HYPRE_SStructGraphSetObjectType(m_ss_graph, m_hypre_object_type);

    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        HYPRE_SStructGraphSetStencil(m_ss_graph, ilev, ivar, m_ss_stencil);
    }

    addNonStencilEntriesToGraph();

    HYPRE_SStructGraphAssemble(m_ss_graph);

    HYPRE_SStructMatrixCreate(m_comm, m_ss_graph, &m_ss_A);
    HYPRE_SStructMatrixSetObjectType(m_ss_A, m_hypre_object_type);
    HYPRE_SStructMatrixInitialize(m_ss_A);
}

HypreMLABecLap::~HypreMLABecLap ()
{
    HYPRE_SStructGridDestroy(m_ss_grid);
    HYPRE_SStructStencilDestroy(m_ss_stencil);
    HYPRE_SStructGraphDestroy(m_ss_graph);
#if 0
    if (m_ss_precond) {
        HYPRE_SStructSolverDestroy(m_ss_precond);
    }
#endif
    HYPRE_SStructMatrixDestroy(m_ss_A);
    if (m_ss_x) {
        HYPRE_SStructVectorDestroy(m_ss_x);
    }
    if (m_ss_b) {
        HYPRE_SStructVectorDestroy(m_ss_b);
    }
    if (m_solver) {
#ifdef AMREX_FEATURE_HYPRE_SSAMG
        if (m_hypre_solver_id == HypreSolverID::SSAMG) {
            if (m_ss_solver) {
                HYPRE_SStructSSAMGDestroy(m_ss_solver);
            }
        } else
#endif
        {
            if (m_solver) {
                HYPRE_BoomerAMGDestroy(m_solver);
            }
        }
    }
}

void HypreMLABecLap::addNonStencilEntriesToGraph ()
{
    BL_PROFILE("HypreMLABecLap::addNonStencilEntriesToGraph");

    Vector<std::tuple<int,int,IntVect,int,IntVect>> entries;

    for (int ilev = 1; ilev < m_nlevels; ++ilev) {
        int const clev = ilev-1;
        int const flev = ilev;

#if (AMREX_SPACEDIM == 3)
        Box const& cgdomain = m_geom[clev].growPeriodicDomain(1);
#endif

        IntVect const& refratio = m_ref_ratio[clev];

        auto const& fine_mask = m_fine_masks[clev];

        for (MFIter mfi(fine_mask, MFItInfo().DisableDeviceSync());
             mfi.isValid(); ++mfi)
        {
            auto const lidx = mfi.LocalIndex();
            int c2f_total_from = 0;
            Long c2f_total_to = 0;
#ifdef AMREX_USE_GPU
            IArrayBox h_mask(fine_mask[mfi].box(), 1, The_Pinned_Arena());
            IArrayBox h_c2f_offset_from(m_c2f_offset_from[clev][mfi].box(), 1, The_Pinned_Arena());
            IArrayBox h_c2f_offset_to(m_c2f_offset_to[clev][mfi].box(), 1, The_Pinned_Arena());
            IArrayBox h_c2f_nentries(m_c2f_nentries[clev][mfi].box(), 1, The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(h_mask.dataPtr(),
                                   fine_mask[mfi].dataPtr(),
                                   h_mask.nBytes());
            Gpu::streamSynchronize();
            auto const& mask = h_mask.const_array();
            auto const& c2f_offset_from = h_c2f_offset_from.array();
            auto const& c2f_offset_to = h_c2f_offset_to.array();
            auto const& c2f_nentries = h_c2f_nentries.array();
#else
            auto const& mask = fine_mask.const_array(mfi);
            auto const& c2f_offset_from = m_c2f_offset_from[clev].array(mfi);
            auto const& c2f_offset_to = m_c2f_offset_to[clev].array(mfi);
            auto const& c2f_nentries = m_c2f_nentries[clev].array(mfi);
#endif
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
            {
                amrex::ignore_unused(k);
                int nc2f = 0;
                if (mask(i,j,k) == 0) { // uncovered coarse cell
                    IntVect const civ(AMREX_D_DECL(i,j,k));
                    for (OrientationIter ori; ori; ++ori) {
                        auto const face = ori();
                        int const idir = face.coordDir();
                        IntVect offset(0);
                        offset[idir] = face.isLow() ? -1 : 1;
                        IntVect const to_civ = civ + offset;
                        if (mask(to_civ) == 1) { // covered by fine cells
                            IntVect lo = to_civ * refratio;
                            IntVect hi = lo + refratio - 1;
                            if (face.isLow()) {
                                lo[idir] = hi[idir] - 1;
                            } else {
                                hi[idir] = lo[idir] + 1;
                            }
                            // [lo,hi]: two layers of adjacent fine cells
                            auto len = hi-lo+1;
                            nc2f += AMREX_D_TERM(len[0], *len[1], *len[2]);
                            amrex::LoopOnCpu(lo.dim3(), hi.dim3(),
                                             [&] (int ii, int jj, int kk)
                            {
                                amrex::ignore_unused(kk);
                                entries.emplace_back(clev, lidx, civ, flev,
                                                     IntVect(AMREX_D_DECL(ii,jj,kk)));
                            });
#if (AMREX_SPACEDIM == 3)
                            int const idir1 = ((idir+1) < AMREX_SPACEDIM)
                                ? idir+1 : idir+1-AMREX_SPACEDIM;
                            int const idir2 = ((idir+2) < AMREX_SPACEDIM)
                                ? idir+2 : idir+2-AMREX_SPACEDIM;
                            IntVect t1(0); t1[idir1] = 1;
                            IntVect t2(0); t2[idir2] = 1;
                            IntVect c1 = civ - t1 - t2;
                            IntVect c2 = civ + t1 - t2;
                            IntVect c3 = civ - t1 + t2;
                            IntVect c4 = civ + t1 + t2;
                            if (mask(c1) == 0 && cgdomain.contains(c1) &&
                                mask(c2) == 0 && cgdomain.contains(c2) &&
                                mask(c3) == 0 && cgdomain.contains(c3) &&
                                mask(c4) == 0 && cgdomain.contains(c4))
                            {
                                entries.emplace_back(clev, lidx, civ, clev, c1);
                                entries.emplace_back(clev, lidx, civ, clev, c2);
                                entries.emplace_back(clev, lidx, civ, clev, c3);
                                entries.emplace_back(clev, lidx, civ, clev, c4);
                                nc2f += 4;
                                // Note a corase will not have fine cells on
                                // both ends of a direction. So we don't
                                // have to worry about duplication.
                            }
#endif
                        }
                    }
                }
                c2f_offset_from(i,j,k) = c2f_total_from;
                if (nc2f > 0) { ++c2f_total_from; }
                c2f_nentries(i,j,k) = nc2f;
                c2f_offset_to(i,j,k) = int(c2f_total_to);
                c2f_total_to += nc2f;
            });
#ifdef AMREX_USE_GPU
            Gpu::htod_memcpy_async(m_c2f_offset_from[clev][mfi].dataPtr(),
                                   h_c2f_offset_from.dataPtr(),
                                   h_c2f_offset_from.nBytes());
            Gpu::htod_memcpy_async(m_c2f_offset_to[clev][mfi].dataPtr(),
                                   h_c2f_offset_to.dataPtr(),
                                   h_c2f_offset_to.nBytes());
            Gpu::htod_memcpy_async(m_c2f_nentries[clev][mfi].dataPtr(),
                                   h_c2f_nentries.dataPtr(),
                                   h_c2f_nentries.nBytes());
            Gpu::streamSynchronize();
#endif
            AMREX_ASSERT(c2f_total_to < Long(std::numeric_limits<int>::max()));
            m_c2f_total_from[clev][mfi] = int(c2f_total_from);
            m_c2f_total_to[clev][mfi] = int(c2f_total_to);
        }

        for (MFIter mfi(m_grids[flev], m_dmap[flev], MFItInfo().DisableDeviceSync());
             mfi.isValid(); ++mfi)
        {
            auto const lidx = mfi.LocalIndex();
            Box const& vbx = mfi.validbox();
            for (OrientationIter ori; ori; ++ori) {
                auto const face = ori();
                int const idir = face.coordDir();
                int const idir1 = ((idir+1) < AMREX_SPACEDIM)
                    ? idir+1 : idir+1-AMREX_SPACEDIM;
#if (AMREX_SPACEDIM == 3)
                int const idir2 = ((idir+2) < AMREX_SPACEDIM)
                    ? idir+2 : idir+2-AMREX_SPACEDIM;
#endif

#ifdef AMREX_USE_GPU
                IArrayBox h_mask(m_bndry[flev]->bndryMasks(face)[mfi].box(), 1, The_Pinned_Arena());
                Gpu::dtoh_memcpy_async(h_mask.dataPtr(),
                                       m_bndry[flev]->bndryMasks(face)[mfi].dataPtr(),
                                       h_mask.nBytes());
                Gpu::streamSynchronize();
                auto const& mask = h_mask.const_array();
#else
                auto const& mask = m_bndry[flev]->bndryMasks(face).const_array(mfi);
#endif

                Box bin = vbx;
                bin.setRange(idir, vbx[face]); // just inside face
                Box bout = amrex::adjCell(vbx, face); // just outside face
                IntVect offset_n = bout.smallEnd() - bin.smallEnd();
                IntVect offset_t1(0), offset_t1_r(0);
                offset_t1  [idir1] = 1;
                offset_t1_r[idir1] = refratio[idir1];
#if (AMREX_SPACEDIM == 3)
                IntVect offset_t2(0), offset_t2_r(0);
                offset_t2  [idir2] = 1;
                offset_t2_r[idir2] = refratio[idir2];
#endif
                amrex::LoopOnCpu(bin, [&] (int i, int j, int k)
                {
#if (AMREX_SPACEDIM == 2)
                    amrex::ignore_unused(k);
#endif
                    IntVect iv_in(AMREX_D_DECL(i,j,k));
                    IntVect iv_out = iv_in + offset_n;
                    if (mask(iv_out) == BndryData::not_covered) {
                        IntVect civ_out = amrex::coarsen(iv_out, refratio);
                        entries.emplace_back(flev, lidx, iv_in, clev, civ_out);

                        if (mask(iv_out+offset_t1_r) == BndryData::not_covered) {
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out+offset_t1);
                        }
                        if (mask(iv_out-offset_t1_r) == BndryData::not_covered) {
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out-offset_t1);
                        }
#if (AMREX_SPACEDIM == 3)
                        if (mask(iv_out+offset_t2_r) == BndryData::not_covered) {
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out+offset_t2);
                        }
                        if (mask(iv_out-offset_t2_r) == BndryData::not_covered) {
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out-offset_t2);
                        }
                        if (mask(iv_out-offset_t1_r-offset_t2_r) == BndryData::not_covered &&
                            mask(iv_out+offset_t1_r-offset_t2_r) == BndryData::not_covered &&
                            mask(iv_out-offset_t1_r+offset_t2_r) == BndryData::not_covered &&
                            mask(iv_out+offset_t1_r+offset_t2_r) == BndryData::not_covered)
                        {
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out-offset_t1-offset_t2);
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out+offset_t1-offset_t2);
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out-offset_t1+offset_t2);
                            entries.emplace_back(flev, lidx, iv_in, clev, civ_out+offset_t1+offset_t2);
                        }
#endif
                    }
                });
            }
        }
    }

    // There are duplicates at corners.
    // After this entries is also sorted.
    amrex::RemoveDuplicates(entries);

    m_f2c_bno.resize(m_nlevels-1);
    m_f2c_cell.resize(m_nlevels-1);
    m_f2c_nentries.resize(m_nlevels-1);
    m_f2c_offset.resize(m_nlevels-1);
    m_f2c_values.resize(m_nlevels-1);

    Vector<IntVect> period(m_nlevels);
    Vector<IntVect> smallend(m_nlevels);
    Vector<IntVect> bigend(m_nlevels);
    for (int ilev = 0; ilev <m_nlevels; ++ilev) {
        period[ilev] = m_geom[ilev].Domain().length();
        smallend[ilev] = m_geom[ilev].Domain().smallEnd();
        bigend[ilev] = m_geom[ilev].Domain().bigEnd();
    }

    for (auto& entry : entries) {
        auto const from_level = std::get<0>(entry);
        auto const   to_level = std::get<3>(entry);

        // HYPRE_Int might be a different type than int
        auto from_iv = std::get<2>(entry);
        auto   to_iv = std::get<4>(entry);
        GpuArray<HYPRE_Int,AMREX_SPACEDIM> from_index{AMREX_D_DECL(from_iv[0],
                                                                   from_iv[1],
                                                                   from_iv[2])};
        GpuArray<HYPRE_Int,AMREX_SPACEDIM> to_index{AMREX_D_DECL(to_iv[0],
                                                                 to_iv[1],
                                                                 to_iv[2])};
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (m_geom[0].isPeriodic(idim)) {
                if (to_index[idim] < smallend[to_level][idim]) {
                    to_index[idim] += period[to_level][idim];
                } else if (to_index[idim] > bigend[to_level][idim]) {
                    to_index[idim] -= period[to_level][idim];
                }
            }
        }
        constexpr int ivar = 0;
        HYPRE_SStructGraphAddEntries(m_ss_graph,
                                     from_level, from_index.data(), ivar,
                                       to_level,   to_index.data(), ivar);

        if (from_level == to_level + 1) {
            auto const bno = std::get<1>(entry);
            if ((! m_f2c_bno[to_level].empty()) &&
                (m_f2c_bno[to_level].back() == bno) &&
                (m_f2c_cell[to_level].back() == from_iv)) {
                ++m_f2c_nentries[to_level].back();
            } else {
                m_f2c_bno[to_level].push_back(bno);
                m_f2c_cell[to_level].push_back(from_iv);
                m_f2c_nentries[to_level].push_back(1);
            }
        }
    }

    for (int clev = 0; clev < int(m_f2c_nentries.size()); ++clev) {
        auto const& nentries = m_f2c_nentries[clev];
        if (!nentries.empty()) {
            auto& offset = m_f2c_offset[clev];
            offset.resize(nentries.size());
#if (__cplusplus >= 201703L) && (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE >= 10)
            // GCC's __cplusplus is not a reliable indication for C++17 support
            std::exclusive_scan(nentries.begin(), nentries.end(), offset.begin(),
                                std::size_t(0), std::plus<std::size_t>{});
#else
            offset[0] = 0;
            auto const offset_size = offset.size();
            if (offset_size > 1) {
                auto const* pin = nentries.data();
                std::partial_sum(pin, pin+offset_size-1, offset.data()+1, std::plus<std::size_t>{});
            }
#endif
            auto nvalues = std::size_t(nentries.back()) + offset.back();
            m_f2c_values[clev].resize(nvalues,Real(0.0));
        }
    }
}

void HypreMLABecLap::setup (Real a_ascalar, Real a_bscalar,
                            Vector<MultiFab const*> const& a_acoefs,
                            Vector<Array<MultiFab const*,AMREX_SPACEDIM>> const& a_bcoefs,
                            Array<LinOpBCType,AMREX_SPACEDIM> const& a_lobc,
                            Array<LinOpBCType,AMREX_SPACEDIM> const& a_hibc,
                            Vector<MultiFab const*> const& a_levelbcdata,
                            std::pair<MultiFab const*, IntVect> const& a_coarse_bc)
{
    BL_PROFILE("HypreMLABecLap::setup");

    constexpr int ncomp = 1;
    constexpr HYPRE_Int ivar = 0;

    m_ascalar = a_ascalar;
    m_bscalar = a_bscalar;
    m_lobc = a_lobc;
    m_hibc = a_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_geom[0].isPeriodic(idim)) {
            AMREX_ALWAYS_ASSERT(a_lobc[idim] == LinOpBCType::Periodic &&
                                a_hibc[idim] == LinOpBCType::Periodic);
        } else {
            AMREX_ALWAYS_ASSERT((a_lobc[idim] == LinOpBCType::Dirichlet ||
                                 a_lobc[idim] == LinOpBCType::Neumann  ) &&
                                (a_hibc[idim] == LinOpBCType::Dirichlet ||
                                 a_hibc[idim] == LinOpBCType::Neumann  ));
        }
    }

    MultiFab empty;

    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        MultiFab const* levelbc;
        if (ilev < a_levelbcdata.size() && a_levelbcdata[ilev]) {
            levelbc = a_levelbcdata[ilev];
        } else {
            levelbc = &empty;
        }

        IntVect br_ref_ratio;

        if (ilev == 0) {
            if (m_grids[0].numPts() != m_geom[0].Domain().numPts()) {
                // Need coarse data for bc
                br_ref_ratio = a_coarse_bc.second;
                AMREX_ALWAYS_ASSERT(br_ref_ratio.allGT(0));
                int const in_rad = 0;
                int const out_rad = 1;
                int const extend_rad = 2;
                BndryRegister crse_br(amrex::coarsen(m_grids[0],br_ref_ratio),
                                      m_dmap[0], in_rad, out_rad, extend_rad, ncomp);
                if (a_coarse_bc.first) {
                    Box const& cbx = amrex::coarsen(m_geom[0].Domain(),br_ref_ratio);
                    crse_br.copyFrom(*a_coarse_bc.first, 0, 0, 0, ncomp,
                                      m_geom[0].periodicity(cbx));
                } else {
                    crse_br.setVal(Real(0.0));
                }
                m_bndry[0]->setBndryValues(crse_br, 0, *levelbc, 0, 0, ncomp, br_ref_ratio,
                                           InterpBndryData::IBD_max_order_DEF, 1);
            } else {
                br_ref_ratio = IntVect(1);
                m_bndry[0]->setPhysBndryValues(*levelbc, 0, 0, ncomp);
            }
        } else {
            br_ref_ratio = m_ref_ratio[ilev-1];
            m_bndry[ilev]->setPhysBndryValues(*levelbc, 0, 0, ncomp);
        }

        RealVect crse_bc_loc;
        m_bndry[ilev]->setLOBndryConds({a_lobc}, {a_hibc}, br_ref_ratio, crse_bc_loc);
    }

    // load matrix

    constexpr HYPRE_Int stencil_size = 2*AMREX_SPACEDIM + 1;
    Array<HYPRE_Int,stencil_size> stencil_indices;
    std::iota(stencil_indices.begin(), stencil_indices.end(), 0);

    // big enough for 3d w/ refratio of 4. worst case: a coarse cell with
    // fine neighbors on 3 faces.
    Vector<HYPRE_Int> nonstencil_indices(4*4*3+3);
    std::iota(nonstencil_indices.begin(), nonstencil_indices.end(), stencil_size);

    BaseFab<GpuArray<Real,stencil_size>> matfab;
    for (int ilev = m_nlevels-1; ilev >= 0; --ilev) {
        auto dx = m_geom[ilev].CellSizeArray();

        for (MFIter mfi(m_grids[ilev], m_dmap[ilev]); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            matfab.resize(vbx);

            Array4<Real const> afab;
            if (ilev < a_acoefs.size() && a_acoefs[ilev]) {
                afab = a_acoefs[ilev]->const_array(mfi);
            }

            GpuArray<Array4<Real const>, AMREX_SPACEDIM> bfabs;
            if (ilev < a_bcoefs.size()) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (a_bcoefs[ilev][idim]) {
                        bfabs[idim] = a_bcoefs[ilev][idim]->const_array(mfi);
                    }
                }
            }

            GpuArray<int,AMREX_SPACEDIM*2> bctype;
            GpuArray<Real,AMREX_SPACEDIM*2> bcl;
            GpuArray<Array4<int const>, AMREX_SPACEDIM*2> bcmsk;
            GpuArray<Array4<Real const>, AMREX_SPACEDIM*2> bcval;
            GpuArray<Array4<Real>, AMREX_SPACEDIM*2> bcrhs;
            for (OrientationIter oit; oit; oit++) {
                Orientation ori = oit();
                int cdir(ori);
                bctype[cdir] = m_bndry[ilev]->bndryConds(mfi)[cdir][0];
                bcl[cdir] = m_bndry[ilev]->bndryLocs(mfi)[cdir];
                bcmsk[cdir] = m_bndry[ilev]->bndryMasks(ori)[mfi].const_array();
                bcval[cdir] = m_bndry[ilev]->bndryValues(ori)[mfi].const_array();
                bcrhs[cdir] = (*m_bndry_rhs[ilev])[ori][mfi].array();
            }

            Real sa = a_ascalar;
            Real sb = a_bscalar;
            const auto boxlo = amrex::lbound(vbx);
            const auto boxhi = amrex::ubound(vbx);
            // Set up stencil part of the matrix
            auto fixed_pt = IntVect::TheMaxVector();
            if (m_is_singular && m_nlevels-1 == ilev) {
                auto const& box0 = m_grids.back()[0];
                fixed_pt = box0.smallEnd() + 1;
                // This cell does not have any non-stencil entries. So it's
                // a good point for fixing singularity.
            }
            amrex::fill(matfab,
            [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,stencil_size>& sten,
                                       int i, int j, int k)
            {
                hypmlabeclap_mat(sten, i, j, k, boxlo, boxhi, sa, afab, sb, dx, bfabs,
                                 bctype, bcl, bcmsk, bcval, bcrhs, ilev, fixed_pt);
            });

            bool need_sync = true;

            // This is always false for the first iteration (ilev == m_nlevels-1).
            if (ilev < m_nlevels-1) {
                // As a coarse level, the coarse cells near the coarse/fine
                // interface have their transverse stencil parts modified
                // because they have participated in the fine flux
                // calculation. For example, suppose there is coarse cell
                // (i,j) in 2d and its left neighbors are two fine cells
                // (assuming the ref_ratio is 2). The computation of the two
                // fine fluxes across the coarse/fine interface involves
                // cells (i,j), (i,j-1) and (i,j+1), if the latter two are
                // uncovered coarse cells. For 3D, this will also involve
                // coarse corners cells in the transverse plane, which are
                // not part of the stencil.
                //
                // So we need to do three things. (1) Add reflux
                // modifications to the stencil entries. (2) Handle
                // non-stencil entries involving fine cells. (3) For 3D, we
                // need handle the non-stencil coarse corner entries in the
                // transverse plane.
                //
                // We have received B coeffs from the fine level in the
                // previous iteration.

                int c2f_total_from = m_c2f_total_from[ilev][mfi];
                int c2f_total_to = m_c2f_total_to[ilev][mfi];
                Gpu::DeviceVector<GpuArray<HYPRE_Int,AMREX_SPACEDIM>> civ(c2f_total_from);
                Gpu::DeviceVector<HYPRE_Int> nentries(c2f_total_from);
                Gpu::DeviceVector<int> entry_offset(c2f_total_from);
                Gpu::DeviceVector<Real> entry_values(c2f_total_to,Real(0.0));
                auto* p_civ = civ.data();
                auto* p_nentries = nentries.data();
                auto* p_entry_offset = entry_offset.data();
                auto* p_entry_values = entry_values.data();
                auto const& c2f_offset_from_a = m_c2f_offset_from[ilev].const_array(mfi);
                auto const& c2f_nentries_a = m_c2f_nentries[ilev].const_array(mfi);
                auto const& c2f_offset_to_a = m_c2f_offset_to[ilev].const_array(mfi);
                auto const& mat_a = matfab.array();
                auto const& fine_mask = m_fine_masks[ilev].const_array(mfi);
                auto const& crse_mask = m_crse_masks[ilev].const_array(mfi);
                AMREX_D_TERM(auto offset_bx_a = m_offset_cf_bcoefs[ilev][0].isDefined()
                                              ? m_offset_cf_bcoefs[ilev][0].const_array(mfi)
                                              : Array4<int const>{};,
                             auto offset_by_a = m_offset_cf_bcoefs[ilev][1].isDefined()
                                              ? m_offset_cf_bcoefs[ilev][1].const_array(mfi)
                                              : Array4<int const>{};,
                             auto offset_bz_a = m_offset_cf_bcoefs[ilev][2].isDefined()
                                              ? m_offset_cf_bcoefs[ilev][2].const_array(mfi)
                                              : Array4<int const>{});
                AMREX_D_TERM(Real const* p_bx = (!m_cf_bcoefs[ilev][0].empty())
                                                ? m_cf_bcoefs[ilev][0][mfi]->data()
                                                : nullptr;,
                             Real const* p_by = (!m_cf_bcoefs[ilev][1].empty())
                                                ? m_cf_bcoefs[ilev][1][mfi]->data()
                                                : nullptr;,
                             Real const* p_bz = (!m_cf_bcoefs[ilev][2].empty())
                                                ? m_cf_bcoefs[ilev][2][mfi]->data()
                                                : nullptr);
                auto rr = m_ref_ratio[ilev];
                amrex::ParallelFor(vbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    hypmlabeclap_c2f(i,j,k,mat_a,
                                     p_civ,p_nentries,p_entry_offset,p_entry_values,
                                     c2f_offset_from_a, c2f_nentries_a,
                                     c2f_offset_to_a, dx, sb,
                                     AMREX_D_DECL(offset_bx_a,offset_by_a,offset_bz_a),
                                     AMREX_D_DECL(p_bx, p_by, p_bz),
                                     fine_mask,rr, crse_mask);
                });
                if (c2f_total_from > 0) {
#ifdef AMREX_USE_GPU
                    Gpu::PinnedVector<Real> h_entry_values(entry_values.size());
                    Gpu::PinnedVector<int> h_entry_offset(entry_offset.size());
                    Gpu::PinnedVector<GpuArray<HYPRE_Int,AMREX_SPACEDIM>> h_civ(civ.size());
                    Gpu::PinnedVector<HYPRE_Int> h_nentries(nentries.size());
                    Gpu::copyAsync(Gpu::deviceToHost,
                                   entry_values.begin(),
                                   entry_values.end(),
                                   h_entry_values.begin());
                    Gpu::copyAsync(Gpu::deviceToHost,
                                   entry_offset.begin(),
                                   entry_offset.end(),
                                   h_entry_offset.begin());
                    Gpu::copyAsync(Gpu::deviceToHost,
                                   civ.begin(),
                                   civ.end(),
                                   h_civ.begin());
                    Gpu::copyAsync(Gpu::deviceToHost,
                                   nentries.begin(),
                                   nentries.end(),
                                   h_nentries.begin());
                    Gpu::streamSynchronize();
                    need_sync = false;
#else
                    auto& h_entry_values = entry_values;
                    auto const& h_entry_offset = entry_offset;
                    auto& h_civ = civ;
                    auto& h_nentries = nentries;
#endif
                    {
                        BL_PROFILE("HYPRE_SStructMatrixSetValues");
                        for (int ientry = 0; ientry < c2f_total_from; ++ientry) {
                            HYPRE_SStructMatrixSetValues(m_ss_A, ilev, h_civ[ientry].data(),
                                                         ivar, h_nentries[ientry],
                                                         nonstencil_indices.data(),
                                                         h_entry_values.data()
                                                         + h_entry_offset[ientry]);
                        }
                    }
                }
            }

            if (need_sync) { Gpu::streamSynchronize(); }

            HYPRE_Int vbxlo[] = {AMREX_D_DECL(vbx.smallEnd(0), vbx.smallEnd(1), vbx.smallEnd(2))};
            HYPRE_Int vbxhi[] = {AMREX_D_DECL(vbx.bigEnd(0), vbx.bigEnd(1), vbx.bigEnd(2))};
            {
                BL_PROFILE("HYPRE_SStructMatrixSetBoxValues");
                HYPRE_SStructMatrixSetBoxValues(m_ss_A, ilev, vbxlo, vbxhi, ivar, stencil_size,
                                                stencil_indices.data(), (Real*)matfab.dataPtr());
                Gpu::hypreSynchronize();
            }
        }

        if (ilev > 0) {
            // As a fine level, we have non-stencil entries involving coarse
            // cells.

            int const flev = ilev;
            int const clev = ilev - 1;

            // Update m_f2c_values
            auto const num_f2c_cell = int(m_f2c_cell[clev].size());
            if (num_f2c_cell > 0) {
#ifdef AMREX_USE_GPU
                Gpu::DeviceVector<int> d_f2c_bno(m_f2c_bno[clev].size());
                Gpu::DeviceVector<IntVect> d_f2c_cell(m_f2c_cell[clev].size());
                Gpu::DeviceVector<std::size_t> d_f2c_offset(m_f2c_offset[clev].size());
                Gpu::DeviceVector<Real> d_f2c_values(m_f2c_values[clev].size());
                Gpu::copyAsync(Gpu::hostToDevice,
                               m_f2c_bno[clev].begin(),
                               m_f2c_bno[clev].end(),
                               d_f2c_bno.begin());
                Gpu::copyAsync(Gpu::hostToDevice,
                               m_f2c_cell[clev].begin(),
                               m_f2c_cell[clev].end(),
                               d_f2c_cell.begin());
                Gpu::copyAsync(Gpu::hostToDevice,
                               m_f2c_offset[clev].begin(),
                               m_f2c_offset[clev].end(),
                               d_f2c_offset.begin());
                auto const* p_f2c_bno = d_f2c_bno.data();
                auto const* p_f2c_cell = d_f2c_cell.data();
                auto const* p_f2c_offset = d_f2c_offset.data();
                auto* p_f2c_values = d_f2c_values.data();
#else
                auto const* p_f2c_bno = m_f2c_bno[clev].data();
                auto const* p_f2c_cell = m_f2c_cell[clev].data();
                auto const* p_f2c_offset = m_f2c_offset[clev].data();
                auto* p_f2c_values = m_f2c_values[clev].data();
#endif
                AMREX_D_TERM(MultiArray4<Real const> bcx;,
                             MultiArray4<Real const> bcy;,
                             MultiArray4<Real const> bcz);
                if (flev < a_bcoefs.size() && (a_bcoefs[flev][0] != nullptr)) {
                    AMREX_D_TERM(bcx = a_bcoefs[flev][0]->const_arrays();,
                                 bcy = a_bcoefs[flev][1]->const_arrays();,
                                 bcz = a_bcoefs[flev][2]->const_arrays());
                }
                GpuArray<MultiArray4<int const>,AMREX_SPACEDIM*2> bmasks;
                for (OrientationIter ori; ori; ++ori) {
                    auto const face = ori();
                    bmasks[face] = m_bndry[flev]->bndryMasks(face).const_arrays();
                }
                auto rr = m_ref_ratio[clev];
                int not_covered = BndryData::not_covered;
                amrex::ParallelFor(num_f2c_cell, [=] AMREX_GPU_DEVICE (int icell)
                {
                    int const bno = p_f2c_bno[icell];
                    // bcoefs may not exist
                    hypmlabeclap_f2c_set_values
                        (p_f2c_cell[icell],
                         p_f2c_values + p_f2c_offset[icell],
                         dx, a_bscalar,
                         {AMREX_D_DECL(bcx ? bcx[bno] : Array4<Real const>{},
                                       bcy ? bcy[bno] : Array4<Real const>{},
                                       bcz ? bcz[bno] : Array4<Real const>{})},
                         {AMREX_D_DECL(bmasks[0][bno],
                                       bmasks[1][bno],
                                       bmasks[2][bno]),
                          AMREX_D_DECL(bmasks[AMREX_SPACEDIM][bno],
                                       bmasks[AMREX_SPACEDIM+1][bno],
                                       bmasks[AMREX_SPACEDIM+2][bno])},
                         rr, not_covered);
                });
#ifdef AMREX_USE_GPU
                Gpu::copyAsync(Gpu::deviceToHost,
                               d_f2c_values.begin(),
                               d_f2c_values.end(),
                               m_f2c_values[clev].begin());
                Gpu::streamSynchronize();
#endif

                // This sets non-stencil part for fine cells adjacent to
                // coarse/fine interface.
                for (int i = 0; i < num_f2c_cell; ++i) {
                    auto const& iv = m_f2c_cell[clev][i];
                    HYPRE_Int index[] = {AMREX_D_DECL(iv[0],iv[1],iv[2])};
                    auto const nentries = m_f2c_nentries[clev][i];
                    auto const offset = m_f2c_offset[clev][i];
                    auto* values = m_f2c_values[clev].data() + offset;
                    HYPRE_SStructMatrixSetValues(m_ss_A, flev, index, ivar, nentries,
                                                 nonstencil_indices.data(), values);
                }
            }

            if (ilev < a_bcoefs.size() && a_bcoefs[ilev][0]) {
                commBCoefs(ilev, a_bcoefs[ilev]);
            }
        }
    }

    {
        BL_PROFILE("HYPRE_SStructMatrixAssemble");
        HYPRE_SStructMatrixAssemble(m_ss_A);
        // HYPRE_SStructMatrixPrint("mat", m_ss_A, 0);
    }

#ifdef AMREX_FEATURE_HYPRE_SSAMG
    if (m_hypre_solver_id == HypreSolverID::SSAMG)
    {
        BL_PROFILE("HYPRE_SSAMG_setup");

        AMREX_ALWAYS_ASSERT(m_solver == nullptr);

        HYPRE_SStructSSAMGCreate(m_comm, &m_ss_solver);

        HYPRE_SStructSSAMGSetNumPreRelax(m_ss_solver, 4);
        HYPRE_SStructSSAMGSetNumPostRelax(m_ss_solver, 4);
        HYPRE_SStructSSAMGSetNumCoarseRelax(m_ss_solver, 4);

        HYPRE_SStructSSAMGSetLogging(m_ss_solver, 1);
        // HYPRE_SStructSSAMGSetPrintLevel(m_ss_solver, 1); /* 0: no, 1: setup, 2: solve, 3:both

//        HYPRE_SStructSSAMGSetup(m_ss_solver, A, b, x);

        HYPRE_SStructSSAMGSetMaxIter(m_ss_solver, m_maxiter);
    } else
#endif
    {
        BL_PROFILE("HYPRE_BoomerAMG_setup");

        AMREX_ALWAYS_ASSERT(m_solver == nullptr);

        HYPRE_BoomerAMGCreate(&m_solver);

        HYPRE_BoomerAMGSetOldDefault(m_solver); // Falgout coarsening with modified classical interpolation
        HYPRE_BoomerAMGSetStrongThreshold(m_solver, (AMREX_SPACEDIM == 3) ? 0.4 : 0.25); // default is 0.25
        HYPRE_BoomerAMGSetRelaxOrder(m_solver, 1);   /* 0: default, natural order, 1: C/F relaxation order */
        HYPRE_BoomerAMGSetNumSweeps(m_solver, 2);   /* Sweeps on fine levels */
        // HYPRE_BoomerAMGSetFCycle(m_solver, 1); // default is 0
        // HYPRE_BoomerAMGSetCoarsenType(m_solver, 6);
        // HYPRE_BoomerAMGSetRelaxType(m_solver, 6);   /* G-S/Jacobi hybrid relaxation */

        HYPRE_BoomerAMGSetLogging(m_solver, 1);
        // HYPRE_BoomerAMGSetPrintLevel(m_solver, 1); /* 0: no, 1: setup, 2: solve, 3:both

        HYPRE_ParCSRMatrix par_A;
        HYPRE_SStructMatrixGetObject(m_ss_A, (void**) &par_A);
        HYPRE_BoomerAMGSetup(m_solver, par_A, nullptr, nullptr);

        HYPRE_BoomerAMGSetMinIter(m_solver, 1);
        HYPRE_BoomerAMGSetMaxIter(m_solver, m_maxiter);
    }
}

void HypreMLABecLap::solve (Vector<MultiFab*> const& a_sol, Vector<MultiFab const*> const& a_rhs,
                            Real a_reltol, Real a_abstol)
{
    BL_PROFILE("HypreMLABecLap::solve()");

    constexpr int ncomp = 1;
    constexpr HYPRE_Int ivar = 0;

    // Load vectors
    {
        BL_PROFILE("HypreMLABecLap::load_vector");

        // Do we still have to do this repeatedly to avoid a hypre bug?
        HYPRE_SStructVectorCreate(m_comm, m_ss_grid, &m_ss_x);
        HYPRE_SStructVectorSetObjectType(m_ss_x, m_hypre_object_type);
        HYPRE_SStructVectorInitialize(m_ss_x);
        //
        HYPRE_SStructVectorCreate(m_comm, m_ss_grid, &m_ss_b);
        HYPRE_SStructVectorSetObjectType(m_ss_b, m_hypre_object_type);
        HYPRE_SStructVectorInitialize(m_ss_b);

        for (int ilev = 0; ilev < m_nlevels; ++ilev) {
            FArrayBox tmp;
            bool has_ghostcells = (a_sol[ilev]->nGrowVect() != 0);
            AMREX_ALWAYS_ASSERT(a_rhs[ilev]->nGrowVect() == 0);
            for (MFIter mfi(*a_rhs[ilev]); mfi.isValid(); ++mfi) {
                Box const& vbx = mfi.validbox();
                HYPRE_Int vbxlo[] = {AMREX_D_DECL(vbx.smallEnd(0), vbx.smallEnd(1), vbx.smallEnd(2))};
                HYPRE_Int vbxhi[] = {AMREX_D_DECL(vbx.bigEnd(0), vbx.bigEnd(1), vbx.bigEnd(2))};

                FArrayBox& solsrc = (*a_sol[ilev])[mfi];
                Real* psol;
                if (has_ghostcells) {
                    tmp.resize(vbx, ncomp);
                    psol = tmp.dataPtr();
                    solsrc.template copyToMem<RunOn::Device>(vbx, 0, ncomp, psol);
                    Gpu::streamSynchronize();
                } else {
                    psol = solsrc.dataPtr();
                }
                HYPRE_SStructVectorSetBoxValues(m_ss_x, ilev, vbxlo, vbxhi, ivar, psol);
                Gpu::hypreSynchronize();

                tmp.resize(vbx, ncomp);
                auto const& rhs1 = tmp.array();
                auto const& rhs0 = a_rhs[ilev]->const_array(mfi);
                auto* prhs = tmp.dataPtr();
                GpuArray<Array4<int const>, AMREX_SPACEDIM*2> bcmsk;
                GpuArray<Array4<Real const>, AMREX_SPACEDIM*2> bcrhs;
                for (OrientationIter oit; oit; oit++) {
                    Orientation ori = oit();
                    int cdir(ori);
                    bcmsk[cdir] = m_bndry[ilev]->bndryMasks(ori)[mfi].const_array();
                    bcrhs[cdir] = (*m_bndry_rhs[ilev])[ori][mfi].const_array();
                }
                const auto boxlo = amrex::lbound(vbx);
                const auto boxhi = amrex::ubound(vbx);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    hypmlabeclap_rhs(i, j, k, boxlo, boxhi, rhs1, rhs0, bcmsk, bcrhs);
                });
                Gpu::streamSynchronize();
                HYPRE_SStructVectorSetBoxValues(m_ss_b, ilev, vbxlo, vbxhi, ivar, prhs);
                Gpu::hypreSynchronize();
            }
        }

        HYPRE_SStructVectorAssemble(m_ss_x);
        HYPRE_SStructVectorAssemble(m_ss_b);

//        HYPRE_SStructVectorPrint("x", m_ss_x, 0);
//        HYPRE_SStructVectorPrint("b", m_ss_b, 0);
    }

    // Solve
    {
        BL_PROFILE("HypreMLABecLap::actual_solve");

        auto reltol = a_reltol;
        if (a_abstol > Real(0.0)) {
            amrex::Abort("HypreMLABecLap::solve: TODO abstol > 0");
        }

        HYPRE_Int num_iterations;
        Real final_res;

#ifdef AMREX_FEATURE_HYPRE_SSAMG
        if (m_hypre_solver_id == HypreSolverID::SSAMG)
        {
            HYPRE_SStructSSAMGSetTol(m_ss_solver, reltol);

            HYPRE_SStructSSAMGSetup(m_ss_solver, m_ss_A, m_ss_b, m_ss_x);

            HYPRE_SStructSSAMGSolve(m_ss_solver, m_ss_A, m_ss_b, m_ss_x);

            HYPRE_SStructSSAMGGetNumIterations(m_ss_solver, &num_iterations);
            HYPRE_SStructSSAMGGetFinalRelativeResidualNorm(m_ss_solver, &final_res);

            if (m_verbose) {
                amrex::Print() << "\n" << num_iterations
                               << " Hypre SSAMG Iterations, Relative Residual "
                               << final_res << '\n';
            }
        } else
#endif
        {
            HYPRE_BoomerAMGSetTol(m_solver, reltol);

            HYPRE_ParCSRMatrix par_A;
            HYPRE_ParVector par_b;
            HYPRE_ParVector par_x;

            HYPRE_SStructMatrixGetObject(m_ss_A, (void**) &par_A);
            HYPRE_SStructVectorGetObject(m_ss_b, (void**) &par_b);
            HYPRE_SStructVectorGetObject(m_ss_x, (void**) &par_x);

            HYPRE_BoomerAMGSolve(m_solver, par_A, par_b, par_x);

            HYPRE_BoomerAMGGetNumIterations(m_solver, &num_iterations);
            HYPRE_BoomerAMGGetFinalRelativeResidualNorm(m_solver, &final_res);

            if (m_verbose) {
                amrex::Print() << "\n" << num_iterations
                               << " Hypre SS BoomerAMG Iterations, Relative Residual "
                               << final_res << '\n';
            }
        }

        if (final_res > reltol) {
            amrex::Abort("Hypre failed to converge after "+std::to_string(num_iterations)+
                         " iterations. Final relative residual is "+std::to_string(final_res));
        }
    }

    // Get solution
    {
        BL_PROFILE("HypreMLABecLap::get_solution");

        HYPRE_SStructVectorGather(m_ss_x);

        for (int ilev = 0; ilev < m_nlevels; ++ilev) {
            FArrayBox sol;
            bool has_ghostcells = (a_sol[ilev]->nGrowVect() != 0);
            for (MFIter mfi(*a_rhs[ilev]); mfi.isValid(); ++mfi) {
                Box const& vbx = mfi.validbox();
                HYPRE_Int vbxlo[] = {AMREX_D_DECL(vbx.smallEnd(0), vbx.smallEnd(1), vbx.smallEnd(2))};
                HYPRE_Int vbxhi[] = {AMREX_D_DECL(vbx.bigEnd(0), vbx.bigEnd(1), vbx.bigEnd(2))};
                FArrayBox& dest = (*a_sol[ilev])[mfi];
                Real* p;
                if (has_ghostcells) {
                    sol.resize(vbx, ncomp);
                    p = sol.dataPtr();
                } else {
                    p = dest.dataPtr();
                }

                HYPRE_SStructVectorGetBoxValues(m_ss_x, ilev, vbxlo, vbxhi, ivar, p);
                Gpu::hypreSynchronize();

                if (has_ghostcells) {
                    dest.template copyFromMem<RunOn::Device>(vbx, 0, ncomp, p);
                    Gpu::streamSynchronize();
                }
            }
        }

        HYPRE_SStructVectorDestroy(m_ss_x);
        HYPRE_SStructVectorDestroy(m_ss_b);
        m_ss_x = nullptr;
        m_ss_b = nullptr;
    }

    for (int ilev = m_nlevels-2; ilev >= 0; --ilev) {
        amrex::average_down(*a_sol[ilev+1], *a_sol[ilev], 0, ncomp, m_ref_ratio[ilev]);
    }
}

#ifdef AMREX_USE_GPU
namespace {
    struct BCCommTag
    {
        Array4<Real const> fsrc;
        Array4<int const> offset;
        Real* pdst;
        Box cbx;
        IntVect d2s;
        int idir;

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box const& box () const noexcept { return cbx; }
    };

    void unpack_bc (Vector<BCCommTag> const& tags, IntVect const& rr)
    {
        amrex::ParallelFor(tags,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, BCCommTag const& tag)
        {
            if (tag.offset(i,j,k) >= 0) {
                IntVect rrface = rr;
                rrface[tag.idir] = 1;
                IntVect fiv(AMREX_D_DECL(i*rr[0]+tag.d2s[0],
                                         j*rr[1]+tag.d2s[1],
                                         k*rr[2]+tag.d2s[2]));
                auto* p = tag.pdst + tag.offset(i,j,k);
#if (AMREX_SPACEDIM == 3)
                for (int irz = 0; irz < rrface[2]; ++irz) {
#endif
                for (int iry = 0; iry < rrface[1]; ++iry) {
                for (int irx = 0; irx < rrface[0]; ++irx) {
                    *p++ = tag.fsrc(fiv+IntVect(AMREX_D_DECL(irx,iry,irz)));
                }}
#if (AMREX_SPACEDIM == 3)
                }
#endif
            }
        });
    }
}
#endif

void HypreMLABecLap::commBCoefs (int flev, Array<MultiFab const*,AMREX_SPACEDIM> const& a_bcoefs)
{
    AMREX_ASSERT(AMREX_D_TERM(a_bcoefs[0], && a_bcoefs[1], && a_bcoefs[2]));

    int const ncomp = 1;
    int const clev = flev-1;
    IntVect const& rr = m_ref_ratio[clev];

    BoxArray const& ba_dst =                amrex::convert(m_grids[clev],IntVect(1));
    BoxArray const& ba_src = amrex::coarsen(amrex::convert(m_grids[flev],IntVect(1)),rr);
    DistributionMapping const& dm_dst = m_dmap[clev];
    DistributionMapping const& dm_src = m_dmap[flev];
    auto const& cperiod = m_geom[clev].periodicity();

    auto & offset_bcoefs = m_offset_cf_bcoefs[clev];
    auto & cf_bcoefs = m_cf_bcoefs[clev];
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        auto const& ba = amrex::convert(ba_dst,IntVect::TheDimensionVector(idim));
        offset_bcoefs[idim].define(ba, dm_dst, 1, 0);
        cf_bcoefs[idim].define(ba, dm_dst);

        IntVect rrface = rr;
        rrface[idim] = 1;
        int const nfaces = AMREX_D_TERM(rrface[0],*rrface[1],*rrface[2]);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(offset_bcoefs[idim]); mfi.isValid(); ++mfi) {
            auto const& offset_a = offset_bcoefs[idim].array(mfi);
            auto const& fmask_a = m_fine_masks[clev].const_array(mfi);
            BoxIndexer box_indexer(mfi.validbox());
#ifdef AMREX_USE_GPU
            auto tot = Scan::PrefixSum<int>
                (std::uint64_t(mfi.validbox().numPts()),
                 [=] AMREX_GPU_DEVICE (std::uint64_t index) -> int
                 {
                     auto ijk = box_indexer(index);
                     IntVect ivm(ijk);
                     ivm[idim] -= 1;
                     return int(fmask_a(ijk) != fmask_a(ivm));
                 },
                 [=] AMREX_GPU_DEVICE (Long index, int psum)
                 {
                     auto ijk = box_indexer(index);
                     IntVect ivm(ijk);
                     ivm[idim] -= 1;
                     if (fmask_a(ijk) == fmask_a(ivm)) {
                         offset_a(ijk) = -1; // not a coarse/fine face
                     } else {
                         offset_a(ijk) = psum*nfaces;
                     }
                 },
                 Scan::Type::exclusive, Scan::retSum);
#else
            int tot = 0;
            amrex::Loop(mfi.validbox(), [&] (int i, int j, int k)
            {
                IntVect ivm(AMREX_D_DECL(i,j,k));
                ivm[idim] -= 1;
                int is_cf = (fmask_a(i,j,k) != fmask_a(ivm));
                int psum = tot;
                tot += int(is_cf);
                offset_a(i,j,k) = is_cf ? psum*nfaces : -1;
            });
#endif
            tot *= nfaces;
            cf_bcoefs[idim][mfi] = std::make_unique<Gpu::DeviceVector<Real>>(tot);
        }
    }

    using Tag = FabArrayBase::CopyComTag;
    Vector<Tag> loc_tags;
    FabArrayBase::MapOfCopyComTagContainers send_tags;
    FabArrayBase::MapOfCopyComTagContainers recv_tags;

    auto const& imap_src = a_bcoefs[0]->IndexArray();
    auto const& imap_dst = m_fine_masks[clev].IndexArray();

    int const myproc = ParallelDescriptor::MyProc();
    AMREX_ALWAYS_ASSERT(ParallelDescriptor::TeamSize() == 1);

    if (!(imap_dst.empty() && imap_src.empty())) {
        auto const nlocal_src = static_cast<int>(imap_src.size());
        auto const nlocal_dst = static_cast<int>(imap_dst.size());

        std::vector<std::pair<int,Box>> isects;
        auto const& cpshifts = cperiod.shiftIntVect();

        if (ParallelContext::NProcsSub() > 1) {
            for (int i = 0; i < nlocal_src; ++i) {
                int const k_src = imap_src[i];
                Box const& bx_src = ba_src[k_src];
                for (auto const& pit : cpshifts) {
                    Box const& bx_src_shifted = bx_src + pit;
                    ba_dst.intersections(bx_src_shifted, isects);
                    for (auto const& is : isects) {
                        int const k_dst = is.first;
                        int const dst_owner = dm_dst[k_dst];
                        if (myproc != dst_owner) { // local copy will be dealt with later
                            Box const& bx_dst = ba_dst[k_dst];
                            for (OrientationIter ori; ori; ++ori) {
                                auto const face = ori();
                                int const idir = face.coordDir();
                                auto const ixtype = IntVect::TheDimensionVector(idir);
                                Box const& face_bx_src_shifted = amrex::convert(bx_src_shifted, ixtype);
                                Box const& face_bx_dst = amrex::convert(bx_dst, ixtype);
                                Box const& b = amrex::bdryNode(face_bx_src_shifted, face) & face_bx_dst;
                                if (b.ok()) {
                                    send_tags[dst_owner].emplace_back
                                        (amrex::refine(b,rr), amrex::refine(b-pit,rr), k_dst, k_src);
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < nlocal_dst; ++i) {
            int const k_dst = imap_dst[i];
            Box const& bx_dst = ba_dst[k_dst];
            for (auto const& pit : cpshifts) {
                Box const& bx_dst_shifted = bx_dst + pit;
                ba_src.intersections(bx_dst_shifted, isects);
                for (auto const& is : isects) {
                    int const k_src = is.first;
                    int const src_owner = dm_src[k_src];
                    Box const& bx_src = ba_src[k_src];
                    auto& tagv = (myproc == src_owner) ? loc_tags : recv_tags[src_owner];
                    for (OrientationIter ori; ori; ++ori) {
                        auto const face = ori();
                        int const idir = face.coordDir();
                        auto const ixtype = IntVect::TheDimensionVector(idir);
                        Box const& face_bx_dst_shifted = amrex::convert(bx_dst_shifted, ixtype);
                        Box const& face_bx_src = amrex::convert(bx_src, ixtype);
                        Box const& b = amrex::bdryNode(face_bx_src, face) & face_bx_dst_shifted;
                        if (b.ok()) {
                            tagv.emplace_back(amrex::refine(b-pit,rr),
                                              amrex::refine(b,rr), k_dst, k_src);
                        }
                    }
                }
            }
        }

        // We need to fix the order so that the send and recv processes match.
        auto f = [] (Tag const& a, Tag const& b) {
                     return (a.sbox.ixType() < b.sbox.ixType() )
                         || ((a.sbox.ixType() == b.sbox.ixType()) && (a < b));
                 };
        for (auto& [k, v] : send_tags) {
            std::sort(v.begin(), v.end(), f);
        }
        for (auto& [k, v] : recv_tags) {
            std::sort(v.begin(), v.end(), f);
        }
    }

    if (ParallelContext::NProcsSub() == 1) {
        commBCoefs_local(flev, a_bcoefs, loc_tags);
        return;
    }

#ifdef AMREX_USE_MPI
    int const mpi_tag = ParallelDescriptor::SeqNum();

    auto const N_snds = int(send_tags.size());
    auto const N_rcvs = int(recv_tags.size());
    auto const N_locs = int(loc_tags.size());

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) { return; }

    using FA = FabArray<FArrayBox>;

    CommHandler handler{};
    handler.mpi_tag = mpi_tag;
    if (N_rcvs > 0) {
        auto& comm_data = handler.recv;
        comm_data.the_data = FA::PostRcvs(recv_tags, comm_data.data, comm_data.size, comm_data.rank,
                                          comm_data.request, ncomp, handler.mpi_tag);
    }

    if (N_snds > 0) {
        auto& comm_data = handler.send;
        comm_data.the_data = FA::PrepareSendBuffers(send_tags, comm_data.data, comm_data.size,
                                                    comm_data.rank, comm_data.request, comm_data.cctc,
                                                    ncomp);
        // pack send buffers
#ifdef AMREX_USE_GPU
        Vector<Array4PairTag<Real>> bc_send_tags;
#endif

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel for
#endif
        for (int isend = 0; isend < N_snds; ++isend) {
            char* dptr = comm_data.data[isend];
            auto const& cctc = *comm_data.cctc[isend];
            for (auto const& t : cctc) {
                Box const& bx = t.sbox;
                auto const type = bx.ixType();
#if (AMREX_SPACEDIM == 2)
                int idir = type.nodeCentered(0) ? 0 : 1;
#else
                int idir = type.nodeCentered(0) ? 0 : (type.nodeCentered(1) ? 1 : 2);
#endif
                auto const& sfab = (*a_bcoefs[idir])[t.srcIndex];
#ifdef AMREX_USE_GPU
                bc_send_tags.emplace_back(Array4PairTag<Real>
                    {makeArray4<Real>((Real*)dptr, bx, ncomp), sfab.const_array(), bx});
                std::size_t nbytes = bx.numPts()*sizeof(Real)*ncomp;
#else
                auto nbytes = sfab.copyToMem(bx, 0, ncomp, dptr);
#endif
                dptr += nbytes;
            }
            AMREX_ASSERT(dptr <= comm_data.data[isend] + comm_data.size[isend] && ncomp == 1);
        }

#ifdef AMREX_USE_GPU
        amrex::ParallelFor(bc_send_tags,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, Array4PairTag<Real> const& tag)
        {
            tag.dfab(i,j,k) = tag.sfab(i,j,k);
        });
#endif

        FA::PostSnds(comm_data.data, comm_data.size, comm_data.rank, comm_data.request,
                     handler.mpi_tag);
    }

    if (N_locs > 0) {
        commBCoefs_local(flev, a_bcoefs, loc_tags);
    }

    if (N_rcvs > 0) {
        auto& comm_data = handler.recv;
        comm_data.stats.resize(comm_data.request.size());
        ParallelDescriptor::Waitall(comm_data.request, comm_data.stats);

        // unpack recv buffers
#ifdef AMREX_USE_GPU
        Vector<BCCommTag> bc_recv_tags;
#endif

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel for
#endif
        for (int irecv = 0; irecv < N_rcvs; ++irecv) {
            char const* dptr = comm_data.data[irecv];
            auto const& cctc = recv_tags.at(comm_data.rank[irecv]);
            for (auto const& t : cctc) {
                Box const& bx = t.dbox;
                Box const& cbx = amrex::coarsen(bx,rr);

                auto const type = cbx.ixType();
#if (AMREX_SPACEDIM == 2)
                int idir = type.nodeCentered(0) ? 0 : 1;
#else
                int idir = type.nodeCentered(0) ? 0 : (type.nodeCentered(1) ? 1 : 2);
#endif

                auto const& fsrc = amrex::makeArray4((Real const*)dptr, bx, ncomp);
                auto const& offset = offset_bcoefs[idir].const_array(t.dstIndex);
                auto* pdst = cf_bcoefs[idir][t.dstIndex]->data();

#ifdef AMREX_USE_GPU
                bc_recv_tags.emplace_back(BCCommTag{fsrc, offset, pdst, cbx, IntVect(0), idir});
#else
                IntVect rrface = rr;
                rrface[idir] = 1;
                Dim3 rrdim3 = rr.dim3();
                amrex::LoopOnCpu(cbx, [&] (int i, int j, int k)
                {
                    if (offset(i,j,k) >= 0) {
                        int ii = i*rrdim3.x;
                        int jj = j*rrdim3.y;
                        int kk = k*rrdim3.z;
                        auto* p = pdst + offset(i,j,k);
#if (AMREX_SPACEDIM == 3)
                        for (int irz = 0; irz < rrface[2]; ++irz) {
#else
                        constexpr int irz = 0;
#endif
                        for (int iry = 0; iry < rrface[1]; ++iry) {
                        for (int irx = 0; irx < rrface[0]; ++irx) {
                            *p++ = fsrc(ii+irx,jj+iry,kk+irz);
                        }}
#if (AMREX_SPACEDIM == 3)
                        }
#endif
                    }
                });
#endif
                dptr += bx.numPts() * ncomp * sizeof(Real);
            }
            AMREX_ASSERT(dptr <= comm_data.data[irecv] + comm_data.size[irecv] && ncomp == 1);
        }

#ifdef AMREX_USE_GPU
        unpack_bc(bc_recv_tags, rr);
#endif
    }

    if (N_snds > 0) {
        auto& comm_data = handler.send;
        comm_data.stats.resize(comm_data.request.size());
        ParallelDescriptor::Waitall(comm_data.request, comm_data.stats);
    }
#endif
}

void HypreMLABecLap::commBCoefs_local (int flev,
                                       Array<MultiFab const*,AMREX_SPACEDIM> const& a_bcoefs,
                                       Vector<FabArrayBase::CopyComTag> const& tags)
{
    if (tags.empty()) { return; }

    int const clev = flev-1;
    auto const& rr = m_ref_ratio[clev];

    auto const& offset_cf_bcoefs = m_offset_cf_bcoefs[clev];
    auto const& cf_bcoefs = m_cf_bcoefs[clev];

#ifdef AMREX_USE_GPU
    Vector<BCCommTag> bc_local_tags;
#endif

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel for
#endif
    for (int itag = 0; itag < int(tags.size()); ++itag) { // NOLINT(modernize-loop-convert)
        auto const& tag = tags[itag];
        Box const& cbx = amrex::coarsen(tag.dbox,rr);
        IntVect d2s = tag.dbox.smallEnd() - tag.sbox.smallEnd();

        auto const type = cbx.ixType();
#if (AMREX_SPACEDIM == 2)
        int idir = type.nodeCentered(0) ? 0 : 1;
#else
        int idir = type.nodeCentered(0) ? 0 : (type.nodeCentered(1) ? 1 : 2);
#endif

        auto const& fsrc = a_bcoefs[idir]->const_array(tag.srcIndex);
        auto const& offset = offset_cf_bcoefs[idir].const_array(tag.dstIndex);
        auto* pdst = cf_bcoefs[idir][tag.dstIndex]->data();

#ifdef AMREX_USE_GPU
        bc_local_tags.emplace_back(BCCommTag{fsrc, offset, pdst, cbx, d2s, idir});
#else
        IntVect rrface = rr;
        rrface[idir] = 1;
        amrex::LoopOnCpu(cbx, [&] (int i, int j, int k)
        {
            if (offset(i,j,k) >= 0) {
                IntVect fiv(AMREX_D_DECL(i*rr[0]+d2s[0],
                                         j*rr[1]+d2s[1],
                                         k*rr[2]+d2s[2]));
                auto* p = pdst + offset(i,j,k);
#if (AMREX_SPACEDIM == 3)
                for (int irz = 0; irz < rrface[2]; ++irz) {
#endif
                for (int iry = 0; iry < rrface[1]; ++iry) {
                for (int irx = 0; irx < rrface[0]; ++irx) {
                   *p++ = fsrc(fiv+IntVect(AMREX_D_DECL(irx,iry,irz)));
                }}
#if (AMREX_SPACEDIM == 3)
                }
#endif
            }
        });
#endif
    }

#ifdef AMREX_USE_GPU
    unpack_bc(bc_local_tags, rr);
#endif
}

}
