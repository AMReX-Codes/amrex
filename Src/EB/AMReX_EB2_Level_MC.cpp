#include <AMReX_EB2_Level.H>
#include <AMReX_MarchingCubes.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelReduce.H>

#include <algorithm>
#include <array>
#include <limits>
#include <string>

/**
 * \file AMReX_EB2_Level_MC.cpp
 *
 * Geometry-source independent driver of the marching-cubes EB builder.  The
 * only geometry-specific work (filling the nodal MC level set and the exact
 * edge crossings) happens in GShopLevel<G>::define_fine_marching_cubes; see
 * AMReX_EB2_Level_MC.H for the adapters.  Everything below consumes just the
 * nodal field in m_sdf and the per-FAB MC::MCFab caches.
 */

namespace amrex::EB2 {

namespace {

// Device kernels live in free functions: CUDA does not allow extended device
// lambdas inside protected member functions.

//! Pre-fill the volume fraction of \p bx from the corner signs: 0 covered,
//! 1 regular, -1 for cut cells that build_cell_fractions overwrites.
void prefill_volume_fractions (Box const& bx, Array4<Real const> const& sdf,
                               Array4<Real> const& vfrac)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int nfluid = 0;
        nfluid += sdf(i, j, k) > 0.0_rt;
        nfluid += sdf(i + 1, j, k) > 0.0_rt;
        nfluid += sdf(i, j + 1, k) > 0.0_rt;
        nfluid += sdf(i + 1, j + 1, k) > 0.0_rt;
        nfluid += sdf(i, j, k + 1) > 0.0_rt;
        nfluid += sdf(i + 1, j, k + 1) > 0.0_rt;
        nfluid += sdf(i, j + 1, k + 1) > 0.0_rt;
        nfluid += sdf(i + 1, j + 1, k + 1) > 0.0_rt;
        vfrac(i, j, k) = nfluid == 0 ? 0.0_rt : (nfluid == 8 ? 1.0_rt : -1.0_rt);
    });
}

//! Flip the sign of the nodal field: EB2's public convention is negative in fluid.
void negate_levelset (Box const& bx, Array4<Real> const& phi)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        phi(i, j, k) = -phi(i, j, k);
    });
}

//! Per-FAB counter blocks (MC::num_fab_counters ints each) for one level.
class FabCounters
{
public:
    explicit FabCounters (int nfabs)
        : m_buffer(static_cast<std::size_t>(nfabs) * MC::num_fab_counters)
    {}

    //! Zero every block on the host and push to the device.
    void reset ()
    {
        std::fill_n(m_buffer.hostData(), m_buffer.size(), 0);
        m_buffer.copyToDeviceAsync();
        // The FAB loops that follow may use other streams.
        Gpu::streamSynchronize();
    }

    //! Device pointer to the block of the FAB with local index \p local_index.
    int* block (int local_index) noexcept
    {
        return m_buffer.data() + static_cast<std::size_t>(local_index) * MC::num_fab_counters;
    }

    //! Copy every block to the host and reduce over FABs and MPI ranks.
    void reduce ()
    {
        Gpu::streamSynchronizeAll();
        m_buffer.copyToHost();
        m_totals.fill(0);
        for (std::size_t n = 0; n < m_buffer.size(); ++n) {
            m_totals[n % MC::num_fab_counters] += m_buffer.hostData()[n];
        }
        ParallelAllReduce::Sum(m_totals.data(), MC::num_fab_counters,
                               ParallelContext::CommunicatorSub());
    }

    //! Level-wide total of \p counter after reduce().
    [[nodiscard]] int total (MC::FabCounter counter) const noexcept { return m_totals[counter]; }

private:
    Gpu::Buffer<int> m_buffer;
    std::array<int, MC::num_fab_counters> m_totals{};
};

} // namespace

void
Level::assert_marching_cubes_supported (Geometry const& geom)
{
    auto const cell_size = geom.CellSizeArray();
    Real const max_cell_size = amrex::max(cell_size[0], amrex::max(cell_size[1], cell_size[2]));
    Real const cubic_tolerance = 16.0_rt * std::numeric_limits<Real>::epsilon() * max_cell_size;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.Coord() == CoordSys::cartesian &&
            std::abs(cell_size[0] - cell_size[1]) <= cubic_tolerance &&
            std::abs(cell_size[0] - cell_size[2]) <= cubic_tolerance,
        "Marching-cubes EB construction requires a 3D Cartesian grid with dx == "
        "dy == dz");
}

LayoutData<MC::MCFab>
Level::define_marching_cubes_caches ()
{
    BL_PROFILE("EB2::Level::define_marching_cubes_caches");

    MC::Initialize();

    LayoutData<MC::MCFab> mc_fabs(m_grids, m_dmap);
    for (MFIter mfi(m_sdf, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
        mc_fabs[mfi].defineEdgeIntersections(m_sdf[mfi].box());
    }
    return mc_fabs;
}

void
Level::build_marching_cubes_level (Geometry const& geom, bool extend_domain_face,
                                   LayoutData<MC::MCFab>& mc_fabs, bool write_stl_output)
{
    BL_PROFILE("EB2::Level::build_marching_cubes_level");

    auto const repair = GetRepairParameters();
    Real const small_volfrac = repair.small_volfrac;
    int const maxiter = repair.maxiter;
    bool const cover_multiple_cuts = repair.cover_multiple_cuts;

    // The geometry source has filled m_sdf (MC convention: > 0 fluid) and the
    // exact edge crossings from the raw geometry.  The domain-face extension
    // is owned by the builder so that every geometry source behaves
    // identically: it is applied to the nodal field and to the crossings.
    // It runs on whole FABs, ghost nodes included, so every node stays a
    // deterministic function of the geometry and no communication is needed.
    Box const domain = geom.Domain();
    GpuArray<int, 3> const is_periodic{geom.isPeriodic(0), geom.isPeriodic(1), geom.isPeriodic(2)};
    // Every MC entry point accumulates its counts into the FAB's device
    // counter block; the blocks are copied to the host once per pass.
    FabCounters counters(m_sdf.local_size());
    counters.reset();
    if (extend_domain_face) {
        for (MFIter mfi(m_sdf, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
            MC::extend_domain_face_levelset(m_sdf[mfi].box(), domain, is_periodic, m_sdf[mfi],
                                            counters.block(mfi.LocalIndex()));
            MC::extend_domain_face_edge_intersections(domain, is_periodic, mc_fabs[mfi]);
        }
    }

    // Every EB field is allocated here and populated from the MC working SDF.
    int const ng = GFab::ng;
    MFInfo mf_info;
    mf_info.SetTag("EB2::Level-MC");
    m_cellflag.define(m_grids, m_dmap, 1, ng, mf_info);
    m_volfrac.define(m_grids, m_dmap, 1, ng, mf_info);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    m_bndryarea.define(m_grids, m_dmap, 1, ng, mf_info);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)), m_dmap,
                                1, ng, mf_info);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)), m_dmap,
                                AMREX_SPACEDIM - 1, ng, mf_info);
        IntVect edge_type{1};
        edge_type[idim] = 0;
        m_edgecent[idim].define(amrex::convert(m_grids, edge_type), m_dmap, 1, ng, mf_info);
    }
    // prefill_volume_fractions rewrites grow(vbx,1) every pass; this only
    // provides the regular default beyond it.
    m_volfrac.setVal(1.0_rt, 0, 1, m_volfrac.nGrowVect());

    // Like the moments, the rejection marks are computed on grow(vbx,1), so
    // the ghost ring of these arrays holds exactly the marks the neighboring
    // FAB computes for its valid cells and faces, and the nodal repair needs
    // no communication.  Degenerate faces in the outer ring are marked and
    // repaired like any other.
    iMultiFab rejected_cells(m_grids, m_dmap, 1, IntVect(1), mf_info);
    Array<iMultiFab, AMREX_SPACEDIM> rejected_faces;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        rejected_faces[idim].define(m_areafrac[idim].boxArray(), m_dmap, 1, IntVect(1), mf_info);
    }
    bool converged = false;
    for (int iter = 0; iter < maxiter; ++iter) {
        rejected_cells.setVal(0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rejected_faces[idim].setVal(0);
        }

        counters.reset();
        MFItInfo info{};
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
        info.SetDynamic(true);
#pragma omp parallel
#endif
        for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
            MC::marching_cubes(geom, m_sdf[mfi], mc_fabs[mfi], counters.block(mfi.LocalIndex()));
        }

        m_centroid.setVal(0.0_rt, 0, AMREX_SPACEDIM, m_centroid.nGrowVect());
        m_bndryarea.setVal(0.0_rt, 0, 1, m_bndryarea.nGrowVect());
        m_bndrycent.setVal(-1.0_rt, 0, AMREX_SPACEDIM, m_bndrycent.nGrowVect());
        m_bndrynorm.setVal(0.0_rt, 0, AMREX_SPACEDIM, m_bndrynorm.nGrowVect());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
            Box const vbx = amrex::enclosedCells(mfi.validbox());
            Box const gbx = amrex::grow(vbx, 1);
            auto& mc_fab = mc_fabs[mfi];
            int* const fab_counters = counters.block(mfi.LocalIndex());

            prefill_volume_fractions(gbx, m_sdf.const_array(mfi), m_volfrac.array(mfi));

            MC::build_face_fractions(
                gbx, mc_fab, m_sdf[mfi], m_areafrac[0][mfi], m_areafrac[1][mfi],
                m_areafrac[2][mfi], m_facecent[0][mfi], m_facecent[1][mfi], m_facecent[2][mfi],
                rejected_faces[0][mfi], rejected_faces[1][mfi], rejected_faces[2][mfi],
                fab_counters);
            MC::build_edge_centroids(gbx, mc_fab, m_sdf[mfi], m_edgecent[0][mfi],
                                     m_edgecent[1][mfi], m_edgecent[2][mfi]);
            // Moments and rejection marks include the ghost ring so that they
            // are available locally for the repair and the final topology.
            MC::build_cell_fractions(
                gbx, geom, mc_fab, m_sdf[mfi], m_areafrac[0][mfi],
                m_areafrac[1][mfi], m_areafrac[2][mfi], m_volfrac[mfi], m_centroid[mfi],
                m_bndryarea[mfi], m_bndrycent[mfi], m_bndrynorm[mfi], fab_counters);

            MC::mark_faces_for_cleanup(
                gbx, mc_fab, m_sdf[mfi], rejected_faces[0][mfi],
                rejected_faces[1][mfi], rejected_faces[2][mfi], fab_counters);
            MC::mark_cells_for_cleanup(gbx, mc_fab, m_sdf[mfi], m_volfrac[mfi], small_volfrac,
                                       rejected_cells[mfi], fab_counters);
        }

        // One device-to-host copy and one reduction per pass.
        counters.reduce();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(counters.total(MC::counter_invalid_triangles) == 0,
                                         "Marching Cubes: invalid triangle");
        int const face_decision_errors = counters.total(MC::counter_face_decision_errors);
        int const degenerate_faces = counters.total(MC::counter_degenerate_faces);
        int const cell_geometry_errors = counters.total(MC::counter_closure_errors)
            + counters.total(MC::counter_volume_errors)
            + counters.total(MC::counter_centroid_errors)
            + counters.total(MC::counter_area_vector_errors);
        int face_rejections = counters.total(MC::counter_face_rejections);
        int const topology_rejections = counters.total(MC::counter_topology_rejections);
        int const small_cell_rejections = counters.total(MC::counter_small_cell_rejections);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(face_decision_errors == 0,
                                         "Marching-cubes EB encountered " +
                                             std::to_string(face_decision_errors) +
                                             " Cartesian faces whose two cells resolved the "
                                             "MC33 ambiguity differently");
        // Degenerate face polygons were marked for nodal repair alongside the
        // multi-aperture faces.
        face_rejections += degenerate_faces;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            cell_geometry_errors == 0 || topology_rejections + face_rejections > 0,
            "Marching-cubes EB could not map an invalid cell-moment record to "
            "the nodal repair set");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            cover_multiple_cuts || face_rejections + topology_rejections == 0,
            "Marching-cubes EB found " + std::to_string(face_rejections) +
                " unsupported faces and " + std::to_string(topology_rejections) +
                " unsupported cells; "
                "set eb2.cover_multiple_cuts=1 to enable legacy-style nodal "
                "repair");

        if (face_rejections + topology_rejections + small_cell_rejections == 0) {
            converged = true;
            break;
        }

        counters.reset();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
            // Ghost nodes are repaired too, like the nodes of the legacy
            // generator's grown level set: ghost nodes beyond every grid have
            // no owner and must be fixed here for the repair to converge.
            MC::zero_nodes_for_cleanup(
                m_sdf[mfi].box(), rejected_cells[mfi], rejected_faces[0][mfi],
                rejected_faces[1][mfi], rejected_faces[2][mfi], m_sdf[mfi],
                counters.block(mfi.LocalIndex()));
        }
        counters.reduce();
        int const changed_nodes = counters.total(MC::counter_changed_nodes);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(changed_nodes > 0,
                                         "Marching-cubes EB repair found rejected "
                                         "geometry but changed no fluid nodes");

        // Shared valid nodes were repaired from identical marks on every FAB
        // and already agree.  Owned ghost nodes also depend on cells outside
        // grow(vbx,1), so refresh them from their owners, as the legacy
        // generator does after each pass.  The domain-face extension follows
        // so that it reads up-to-date reference nodes on every FAB.
        m_sdf.FillBoundary(geom.periodicity());
        if (extend_domain_face) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
                MC::extend_domain_face_levelset(m_sdf[mfi].box(), domain, is_periodic,
                                                m_sdf[mfi], counters.block(mfi.LocalIndex()));
            }
        }
        if (amrex::Verbose() > 0) {
            if (small_cell_rejections > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " fixed "
                               << small_cell_rejections << " small cells\n";
            }
            if (face_rejections + topology_rejections > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " fixed "
                               << face_rejections << " unsupported faces and "
                               << topology_rejections << " unsupported cells\n";
            }
            if (cell_geometry_errors > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " included "
                               << cell_geometry_errors
                               << " invalid cell-moment records in the repair set (closure="
                               << counters.total(MC::counter_closure_errors)
                               << ", volume=" << counters.total(MC::counter_volume_errors)
                               << ", centroid=" << counters.total(MC::counter_centroid_errors)
                               << ", area-vector="
                               << counters.total(MC::counter_area_vector_errors) << ")\n";
            }
        }
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(converged,
                                     "Marching-cubes EB failed to fix small or unsupported cells");

    MFItInfo final_info{};
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    final_info.SetDynamic(true);
#endif
    m_cellflag.setVal(EBCellFlag::TheDefaultCell());
    int geometry_errors = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+ : geometry_errors)
#endif
    for (MFIter mfi(m_sdf, final_info); mfi.isValid(); ++mfi) {
        Box const vbx = amrex::enclosedCells(mfi.validbox());
        geometry_errors += MC::build_cell_topology(
            vbx, mc_fabs[mfi], m_sdf[mfi], m_cellflag[mfi], m_volfrac[mfi],
            m_areafrac[0][mfi], m_areafrac[1][mfi], m_areafrac[2][mfi]);
        m_cellflag[mfi].resetType(GFab::ng);
    }
    ParallelAllReduce::Sum(geometry_errors, ParallelContext::CommunicatorSub());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geometry_errors == 0,
                                     "Final repaired marching-cubes topology is inconsistent");

    // The converged triangulation and the public EB data now describe the same
    // repaired domain.  Only the finest level is written; coarse levels that
    // are rebuilt from the geometry source would otherwise overwrite it.
    std::string stl_output;
    ParmParse("eb2").query("mc_stl_file", stl_output);
    if (write_stl_output && !stl_output.empty()) {
        MC::write_stl(stl_output, mc_fabs);
    }
    mc_fabs.clear();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_sdf, final_info); mfi.isValid(); ++mfi) {
        negate_levelset(m_sdf[mfi].box(), m_sdf.array(mfi));
    }
    m_levelset = std::move(m_sdf);

    m_ok = true;
}

}
