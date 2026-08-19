
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_STL_utils.H>
#include <AMReX_MarchingCubes.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WriteEBSurface.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

/**
 * Host-only sphere implicit function.  It is deliberately not GPU-callable so
 * that the generic marching-cubes adapters exercise their RunOn::Cpu path
 * (pinned staging of fillFab/getIntercept), the same way a user-defined
 * host functor would.
 */
struct HostSphereIF
{
    Real m_radius;
    RealArray m_center;
    bool m_fluid_inside;

    Real operator() (RealArray const& p) const noexcept
    {
        Real d2 = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            d2 += (p[d]-m_center[d])*(p[d]-m_center[d]);
        }
        Real const v = std::sqrt(d2) - m_radius;
        return m_fluid_inside ? v : -v;
    }
};

/**
 * Union of two overlapping spheres (fluid outside) with the implicit function
 * multiplied by an arbitrary scale.  The marching-cubes reconstruction,
 * including its ambiguous-face decisions and the nodal repair, must not depend
 * on the units of the implicit function.
 */
struct ScaledTwoSphereIF : amrex::GPUable
{
    Real m_scale;

    AMREX_GPU_HOST_DEVICE Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        constexpr Real radius = 0.6_rt;
        constexpr Real offset = 0.35_rt;
        Real const d1 = std::sqrt((x-offset)*(x-offset) + y*y + z*z);
        Real const d2 = std::sqrt((x+offset)*(x+offset) + y*y + z*z);
        // Positive inside the body (EB2 convention), fluid outside.
        return m_scale*amrex::max(radius-d1, radius-d2);
    }

    Real operator() (RealArray const& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
    }
};

/**
 * Coarsely resolved high-frequency gyroid (fluid outside).  Its saddles produce
 * MC33 case-4/10 cells, so the interior (tunnel) test and the nodal repair of
 * tunnel cells are exercised; both must be independent of the scale.
 */
struct ScaledGyroidIF : amrex::GPUable
{
    Real m_scale;

    AMREX_GPU_HOST_DEVICE Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        constexpr Real k = 21.0_rt;
        return m_scale*(std::sin(k*x)*std::cos(k*y) + std::sin(k*y)*std::cos(k*z)
                        + std::sin(k*z)*std::cos(k*x) - 0.1_rt);
    }

    Real operator() (RealArray const& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
    }
};

//! Implicit-function factories for check_scale_invariance().  These are named
//! types rather than function-local lambdas because nvcc rejects an extended
//! __device__ lambda inside a function template instantiated with a type local
//! to a function.
struct MakeTwoSphereShop
{
    auto operator() (Real scale) const { return EB2::makeShop(ScaledTwoSphereIF{{}, scale}); }
};

struct MakeGyroidShop
{
    auto operator() (Real scale) const { return EB2::makeShop(ScaledGyroidIF{{}, scale}); }
};

//! Build \p make_shop(scale) with the implicit function scaled by 1e-6, 1 and
//! 1e6 (and once with a different EB2 box size) and require the same fluid
//! volume and repaired-node count.
template <class MakeShop>
void check_scale_invariance (std::string const& name, int ncell, MakeShop const& make_shop)
{
    Box const domain(IntVect(0), IntVect(ncell - 1));
    Geometry const geom(domain, RealBox({-1.2_rt, -1.2_rt, -1.2_rt}, {1.2_rt, 1.2_rt, 1.2_rt}),
                        0, {0, 0, 0});
    BoxArray ba(domain);
    ba.maxSize(16);
    DistributionMapping const dm(ba);

    // {scale of the implicit function, EB2 box size used for the build}: the
    // reconstruction must not depend on either.
    struct Variant { Real scale; int eb_max_grid_size; };
    Variant const variants[] = {{.scale = 1.e-6_rt, .eb_max_grid_size = 16},
                                {.scale = 1.0_rt, .eb_max_grid_size = 16},
                                {.scale = 1.e6_rt, .eb_max_grid_size = 16},
                                {.scale = 1.0_rt, .eb_max_grid_size = 8}};
    int const saved_max_grid_size = EB2::max_grid_size;
    Real reference_volume = -1.0_rt;
    Long reference_zero_nodes = -1;
    for (auto const& [scale, eb_max_grid_size] : variants) {
        EB2::max_grid_size = eb_max_grid_size;
        EB2::Build(make_shop(scale), geom, 0, 0, 4);
        auto factory = makeEBFabFactory(geom, ba, dm, {1, 1, 1}, EBSupport::full);
        Real const volume = factory->getVolFrac().sum() * geom.CellSize(0) * geom.CellSize(1)
                            * geom.CellSize(2);
        auto const& levelset = factory->getLevelSet();
        Gpu::DeviceScalar<Long> zero_node_count(static_cast<Long>(0));
        Long* const zero_nodes = zero_node_count.dataPtr();
        for (MFIter mfi(levelset); mfi.isValid(); ++mfi) {
            auto const phi = levelset.const_array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (phi(i, j, k) == 0.0_rt) {
                    Gpu::Atomic::AddNoRet(zero_nodes, static_cast<Long>(1));
                }
            });
        }
        Long global_zero_nodes = zero_node_count.dataValue();
        ParallelAllReduce::Sum(global_zero_nodes, ParallelContext::CommunicatorSub());
        amrex::Print() << name << ": scale " << scale << ", eb2.max_grid_size "
                       << eb_max_grid_size << ": fluid volume " << volume
                       << ", zero level-set nodes " << global_zero_nodes << "\n";
        // The variants share every sign decision, so only the root finding
        // (scale) and the summation order (decomposition) may differ at
        // roundoff level.
        if (reference_volume < 0.0_rt) {
            reference_volume = volume;
            reference_zero_nodes = global_zero_nodes;
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::abs(volume - reference_volume) <= 1.e-8_rt * reference_volume,
                "Marching-cubes fluid volume depends on the scale of the implicit "
                "function or on the box decomposition");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                global_zero_nodes == reference_zero_nodes,
                "Marching-cubes nodal repair depends on the scale of the implicit "
                "function or on the box decomposition");
        }
        EB2::IndexSpace::pop();
    }
    EB2::max_grid_size = saved_max_grid_size;
}

void validate_scale_invariance ()
{
    check_scale_invariance("two spheres", 48, MakeTwoSphereShop{});
    // Tunnel-rich, badly under-resolved geometry: needs the nodal repair and
    // takes a long cascade of repair passes.
    ParmParse ppeb2("eb2");
    int cover_multiple_cuts = 0;
    int maxiter = 32;
    ppeb2.queryAdd("cover_multiple_cuts", cover_multiple_cuts);
    ppeb2.queryAdd("maxiter", maxiter);
    ppeb2.add("cover_multiple_cuts", 1);
    ppeb2.add("maxiter", 200);
    check_scale_invariance("gyroid", 24, MakeGyroidShop{});
    ppeb2.add("cover_multiple_cuts", cover_multiple_cuts);
    ppeb2.add("maxiter", maxiter);
}

//! A zeroed marching-cubes counter block for one FAB.
Gpu::Buffer<int> make_fab_counters ()
{
    Gpu::Buffer<int> counters(MC::num_fab_counters);
    std::fill_n(counters.hostData(), counters.size(), 0);
    counters.copyToDeviceAsync();
    return counters;
}

//! Level-wide sum of one counter after copying the block to the host.
int fab_counter (Gpu::Buffer<int>& counters, MC::FabCounter which)
{
    counters.copyToHost();
    return counters.hostData()[which];
}

bool resolved_bottom_face_is_fluid_connected (GpuArray<Real, 8> const& values)
{
    Box const cell_box(IntVect(0), IntVect(0));
    Box const node_box = amrex::surroundingNodes(cell_box);
    FArrayBox sdf(node_box, 1);
    sdf.setVal<RunOn::Device>(-1.0_rt);
    auto const phi = sdf.array();
    Box const cube_nodes = amrex::surroundingNodes(cell_box);
    ParallelFor(cube_nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int const n = 4 * k + (j == 0 ? i : 3 - i);
        phi(i, j, k) = values[n];
    });

    Geometry geom(cell_box, RealBox({0.0_rt, 0.0_rt, 0.0_rt},
                                    {1.0_rt, 1.0_rt, 1.0_rt}), 0, {0, 0, 0});
    MC::MCFab result;
    auto counters = make_fab_counters();
    MC::marching_cubes(geom, sdf, result, counters.data());
    AMREX_ALWAYS_ASSERT(fab_counter(counters, MC::counter_invalid_triangles) == 0);

    GpuArray<int, MC::num_cell_data_components> host{};
    Gpu::dtoh_memcpy(host.data(), result.m_cell_data.dataPtr(), sizeof(host));
    Gpu::streamSynchronize();
    int const bit = 1 << 4; // MC33 face 5: the low-z face.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((host[MC::face_decision_valid_mask] & bit) != 0,
                                     "MC33 did not retain its resolved low-z face decision");
    return (host[MC::face_fluid_connected_mask] & bit) != 0;
}

void validate_mc33_face_decisions ()
{
    MC::Initialize();
    AMREX_ALWAYS_ASSERT(resolved_bottom_face_is_fluid_connected(
        {2.0_rt, -1.0_rt, 2.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt}));
    AMREX_ALWAYS_ASSERT(!resolved_bottom_face_is_fluid_connected(
        {1.0_rt, -2.0_rt, 1.0_rt, -2.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt}));
    AMREX_ALWAYS_ASSERT(resolved_bottom_face_is_fluid_connected(
        {1.0_rt, -1.0_rt, 1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt}));
}

void validate_mc33_vertex_indices ()
{
    Box const cell_box(IntVect(0), IntVect(0));
    Box const node_box = amrex::surroundingNodes(cell_box);
    Box const cube_nodes = amrex::surroundingNodes(cell_box);
    Geometry geom(cell_box, RealBox({0.0_rt, 0.0_rt, 0.0_rt},
                                    {1.0_rt, 1.0_rt, 1.0_rt}), 0, {0, 0, 0});
    for (int mask = 0; mask < 256; ++mask) {
        FArrayBox sdf(node_box, 1);
        sdf.setVal<RunOn::Device>(-1.0_rt);
        auto const phi = sdf.array();
        ParallelFor(cube_nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int const n = 4 * k + (j == 0 ? i : 3 - i);
            Real const magnitude = 1.0_rt + 0.125_rt * static_cast<Real>(n);
            phi(i, j, k) = (mask & (1 << n)) != 0 ? magnitude : -magnitude;
        });

        MC::MCFab result;
        auto counters = make_fab_counters();
        MC::marching_cubes(geom, sdf, result, counters.data());
        AMREX_ALWAYS_ASSERT(fab_counter(counters, MC::counter_invalid_triangles) == 0);

        GpuArray<int, MC::num_cell_data_components> cell_data{};
        Gpu::dtoh_memcpy(cell_data.data(), result.m_cell_data.dataPtr(), sizeof(cell_data));

        int const ntri = cell_data[MC::triangle_count];
        int const interior_count = cell_data[MC::interior_vertex_count];
        int const interior_offset = cell_data[MC::interior_vertex_offset];
        int const nvertices = static_cast<int>(result.m_vertices.x.size());
        int const edge_vertex_count = nvertices - interior_count;
        AMREX_ALWAYS_ASSERT(std::cmp_equal(ntri, result.m_triangles.v1.size()));
        AMREX_ALWAYS_ASSERT(interior_count == 0 || interior_count == 1);
        if (interior_count != 0) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                interior_offset == edge_vertex_count,
                "MC33 interior vertex overlaps the Cartesian-edge vertex block");
            // The interior vertex is the average of the crossings on the 12
            // cube edges, so it must be their centroid and lie inside the cell.
            Gpu::PinnedVector<Real> vx(nvertices), vy(nvertices), vz(nvertices);
            Gpu::dtoh_memcpy(vx.data(), result.m_vertices.x.data(), nvertices * sizeof(Real));
            Gpu::dtoh_memcpy(vy.data(), result.m_vertices.y.data(), nvertices * sizeof(Real));
            Gpu::dtoh_memcpy(vz.data(), result.m_vertices.z.data(), nvertices * sizeof(Real));
            Gpu::streamSynchronize();
            Real cx = 0.0_rt, cy = 0.0_rt, cz = 0.0_rt;
            for (int n = 0; n < edge_vertex_count; ++n) {
                cx += vx[n]; cy += vy[n]; cz += vz[n];
            }
            cx /= static_cast<Real>(edge_vertex_count);
            cy /= static_cast<Real>(edge_vertex_count);
            cz /= static_cast<Real>(edge_vertex_count);
            Real const tol = 16.0_rt * std::numeric_limits<Real>::epsilon();
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::abs(vx[interior_offset] - cx) <= tol &&
                std::abs(vy[interior_offset] - cy) <= tol &&
                std::abs(vz[interior_offset] - cz) <= tol,
                "MC33 interior vertex is not the centroid of the edge crossings");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                vx[interior_offset] >= 0.0_rt && vx[interior_offset] <= 1.0_rt &&
                vy[interior_offset] >= 0.0_rt && vy[interior_offset] <= 1.0_rt &&
                vz[interior_offset] >= 0.0_rt && vz[interior_offset] <= 1.0_rt,
                "MC33 interior vertex lies outside the cell");
        }

        Gpu::PinnedVector<int> v1(ntri);
        Gpu::PinnedVector<int> v2(ntri);
        Gpu::PinnedVector<int> v3(ntri);
        if (ntri != 0) {
            Gpu::dtoh_memcpy(v1.data(), result.m_triangles.v1.data(), ntri * sizeof(int));
            Gpu::dtoh_memcpy(v2.data(), result.m_triangles.v2.data(), ntri * sizeof(int));
            Gpu::dtoh_memcpy(v3.data(), result.m_triangles.v3.data(), ntri * sizeof(int));
        }
        Gpu::streamSynchronize();
        for (int n = 0; n < ntri; ++n) {
            AMREX_ALWAYS_ASSERT(v1[n] >= 0 && v1[n] < nvertices);
            AMREX_ALWAYS_ASSERT(v2[n] >= 0 && v2[n] < nvertices);
            AMREX_ALWAYS_ASSERT(v3[n] >= 0 && v3[n] < nvertices);
        }
    }
}

/**
 * A single-valued cell may hold exactly one face-connected fluid region.  For
 * every corner-sign mask and several magnitude patterns (chosen so that the
 * MC33 face and interior tests take both branches, including the 4.1.2 /
 * 10.1.2 tunnel tilings), run the extraction and the cell rejection rule and
 * check that a cell is rejected exactly when its fluid corners form more than
 * one group under cube-edge adjacency plus MC33-resolved connected ambiguous
 * faces.  In particular every tunnel (one connected surface patch joining two
 * corner groups through the interior) must be rejected.
 */
void validate_cell_topology_rejection ()
{
    Box const cell_box(IntVect(0), IntVect(0));
    Box const node_box = amrex::surroundingNodes(cell_box);
    Geometry const geom(cell_box, RealBox({0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt}),
                        0, {0, 0, 0});
    // MC33 face ids 1..6 in cube-corner order (test_face's A,B,C,D).
    constexpr int face_corner[6][4] = {{0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3},
                                       {3, 7, 4, 0}, {0, 3, 2, 1}, {4, 7, 6, 5}};
    constexpr int edge_lo[12] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3};
    constexpr int edge_hi[12] = {1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7};

    unsigned int seed = 12345U;
    auto next_magnitude = [&seed] () {
        seed = seed * 1664525U + 1013904223U;              // LCG, deterministic
        return 0.05_rt + 3.0_rt * static_cast<Real>(seed >> 8) / 16777216.0_rt;
    };

    // Every mask a few times, plus many draws of the MC33 case-4 masks (two
    // body-diagonal fluid corners), whose interior test selects the 4.1.2
    // tunnel for roughly one draw in a hundred.
    constexpr int case4_masks[8] = {65, 130, 20, 40, 190, 125, 235, 215};
    auto configurations = [&] (int trial, int slot) {
        return trial < 24 ? slot + 1 : case4_masks[slot % 8];
    };
    Long checked = 0, rejected_cells = 0, tunnels = 0, mismatches = 0;
    for (int trial = 0; trial < 24 + 40; ++trial) {
        int const nslots = trial < 24 ? 254 : 8 * 8;
        for (int slot = 0; slot < nslots; ++slot) {
            int const mask = configurations(trial, slot);
            GpuArray<Real, 8> values{};
            for (int n = 0; n < 8; ++n) {
                Real const magnitude = next_magnitude();
                values[n] = (mask & (1 << n)) != 0 ? magnitude : -magnitude;
            }
            FArrayBox sdf(node_box, 1);
            auto const phi = sdf.array();
            ParallelFor(node_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int const n = 4 * k + (j == 0 ? i : 3 - i);
                phi(i, j, k) = values[n];
            });

            MC::MCFab result;
            auto counters = make_fab_counters();
            MC::marching_cubes(geom, sdf, result, counters.data());
            FArrayBox vfrac(cell_box, 1);
            vfrac.setVal<RunOn::Device>(0.5_rt);      // no small-cell or sentinel rejection
            IArrayBox rejected(cell_box, 1);
            rejected.setVal<RunOn::Device>(0);
            MC::mark_cells_for_cleanup(cell_box, result, sdf, vfrac, 0.0_rt, rejected,
                                       counters.data());
            Gpu::streamSynchronize();

            GpuArray<int, MC::num_cell_data_components> cell_data{};
            Gpu::dtoh_memcpy(cell_data.data(), result.m_cell_data.dataPtr(), sizeof(cell_data));
            int rejected_flag = 0;
            Gpu::dtoh_memcpy(&rejected_flag, rejected.dataPtr(), sizeof(int));
            AMREX_ALWAYS_ASSERT(fab_counter(counters, MC::counter_invalid_triangles) == 0);

            // Reference: fluid corner groups joined along edges and across
            // ambiguous faces whose stored MC33 decision is "connected".
            bool fluid[8];
            for (int n = 0; n < 8; ++n) { fluid[n] = values[n] > 0.0_rt; }
            int parent[8];
            for (int n = 0; n < 8; ++n) { parent[n] = n; }
            auto root = [&] (int n) { while (parent[n] != n) { n = parent[n]; } return n; };
            auto join = [&] (int a, int b) { int ra = root(a), rb = root(b); if (ra != rb) { parent[rb] = ra; } };
            for (int e = 0; e < 12; ++e) {
                if (fluid[edge_lo[e]] && fluid[edge_hi[e]]) { join(edge_lo[e], edge_hi[e]); }
            }
            for (int f = 0; f < 6; ++f) {
                int const c0 = face_corner[f][0], c1 = face_corner[f][1];
                int const c2 = face_corner[f][2], c3 = face_corner[f][3];
                if (!(fluid[c0] == fluid[c2] && fluid[c1] == fluid[c3] && fluid[c0] != fluid[c1])) {
                    continue;
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    (cell_data[MC::face_decision_valid_mask] & (1 << f)) != 0,
                    "An ambiguous face has no stored MC33 decision");
                if ((cell_data[MC::face_fluid_connected_mask] & (1 << f)) != 0) {
                    if (fluid[c0]) { join(c0, c2); } else { join(c1, c3); }
                }
            }
            int groups = 0;
            for (int n = 0; n < 8; ++n) { groups += fluid[n] && root(n) == n; }

            // Tunnel: several corner groups but a single connected surface patch.
            int const ntri = cell_data[MC::triangle_count];
            bool single_patch = ntri > 0;
            if (groups > 1 && ntri > 0) {
                Gpu::PinnedVector<int> v1(ntri), v2(ntri), v3(ntri);
                Gpu::dtoh_memcpy(v1.data(), result.m_triangles.v1.data(), ntri * sizeof(int));
                Gpu::dtoh_memcpy(v2.data(), result.m_triangles.v2.data(), ntri * sizeof(int));
                Gpu::dtoh_memcpy(v3.data(), result.m_triangles.v3.data(), ntri * sizeof(int));
                Gpu::streamSynchronize();
                std::vector<int> tparent(ntri);
                for (int t = 0; t < ntri; ++t) { tparent[t] = t; }
                auto troot = [&] (int t) { while (tparent[t] != t) { t = tparent[t]; } return t; };
                for (int a = 0; a < ntri; ++a) {
                    for (int b = a + 1; b < ntri; ++b) {
                        int const va[3] = {v1[a], v2[a], v3[a]};
                        int const vb[3] = {v1[b], v2[b], v3[b]};
                        bool share = false;
                        for (int const x : va) { for (int const y : vb) { share = share || x == y; } }
                        if (share) { int ra = troot(a), rb = troot(b); if (ra != rb) { tparent[rb] = ra; } }
                    }
                }
                int patches = 0;
                for (int t = 0; t < ntri; ++t) { patches += troot(t) == t; }
                single_patch = patches == 1;
                if (single_patch) { ++tunnels; }
            }

            bool const expect_rejected = groups > 1;
            ++checked;
            rejected_cells += rejected_flag != 0;
            if ((rejected_flag != 0) != expect_rejected) {
                ++mismatches;
                amrex::Print() << "  mask " << mask << " trial " << trial << ": groups=" << groups
                               << " rejected=" << rejected_flag << " triangles=" << ntri << "\n";
            }
            if (groups > 1 && single_patch) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rejected_flag != 0,
                                                 "An MC33 tunnel cell was not rejected");
            }
        }
    }
    amrex::Print() << "Cell topology rule: " << checked << " configurations, " << rejected_cells
                   << " rejected, " << tunnels << " tunnel tilings (all rejected)\n";
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mismatches == 0,
                                     "Cell rejection disagrees with the face-connected corner groups");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(tunnels > 0,
                                     "The topology test did not exercise any tunnel tiling");
}

void validate_exact_ambiguous_face_fraction ()
{
    Box const cell_box(IntVect(0), IntVect(0));
    Box const node_box = amrex::surroundingNodes(cell_box);
    FArrayBox sdf(node_box, 1);
    sdf.setVal<RunOn::Device>(1.0_rt);
    auto const phi = sdf.array();
    ParallelFor(amrex::surroundingNodes(cell_box),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (k == 0) {
                        phi(i, j, k) = (i == j) ? 1.0_rt : -1.0_rt;
                    }
                });

    MC::MCFab mc_fab;
    mc_fab.m_cell_data.resize(cell_box, MC::num_cell_data_components);
    mc_fab.m_cell_data.setVal<RunOn::Device>(0);
    auto const cell_data = mc_fab.m_cell_data.array();
    ParallelFor(cell_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cell_data(i, j, k, MC::face_decision_valid_mask) = 1 << 4;
    });

    mc_fab.defineEdgeIntersections(sdf.box());
    auto const exact_x = mc_fab.m_edge_intersections[0].array();
    auto const exact_y = mc_fab.m_edge_intersections[1].array();
    ParallelFor(cell_box, [=] AMREX_GPU_DEVICE(int, int, int) noexcept {
        exact_x(0, 0, 0) = 0.2_rt;
        exact_y(1, 0, 0) = 0.3_rt;
        exact_x(0, 1, 0) = 0.6_rt;
        exact_y(0, 0, 0) = 0.9_rt;
    });

    FArrayBox apx(amrex::surroundingNodes(cell_box, 0), 1);
    FArrayBox apy(amrex::surroundingNodes(cell_box, 1), 1);
    FArrayBox apz(amrex::surroundingNodes(cell_box, 2), 1);
    FArrayBox fcx(amrex::surroundingNodes(cell_box, 0), 2);
    FArrayBox fcy(amrex::surroundingNodes(cell_box, 1), 2);
    FArrayBox fcz(amrex::surroundingNodes(cell_box, 2), 2);
    IArrayBox rejected_x(amrex::surroundingNodes(cell_box, 0), 1);
    IArrayBox rejected_y(amrex::surroundingNodes(cell_box, 1), 1);
    IArrayBox rejected_z(amrex::surroundingNodes(cell_box, 2), 1);
    rejected_x.setVal<RunOn::Device>(0);
    rejected_y.setVal<RunOn::Device>(0);
    rejected_z.setVal<RunOn::Device>(0);
    auto counters = make_fab_counters();
    MC::build_face_fractions(cell_box, mc_fab, sdf, apx, apy, apz, fcx, fcy, fcz, rejected_x,
                             rejected_y, rejected_z, counters.data());
    int const errors = fab_counter(counters, MC::counter_face_decision_errors)
                       + fab_counter(counters, MC::counter_degenerate_faces);

    Gpu::Buffer<Real> area_buffer(1);
    Box const sample_box = amrex::makeSingleCellBox(IntVect(0), apz.box().ixType());
    apz.copyToMem<RunOn::Device>(sample_box, 0, 1, area_buffer.data());
    area_buffer.copyToHost();
    Real const area = area_buffer.hostData()[0];
    AMREX_ALWAYS_ASSERT(errors == 0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(area - 0.23_rt) < 64.0_rt * std::numeric_limits<Real>::epsilon(),
        "Ambiguous MC face fraction did not use exact STL edge crossings");
}

void validate_nodal_stl_face_levelset ()
{
    Box const domain(IntVect(0), IntVect(7));
    Geometry const geom(
        domain,
        RealBox({-0.4875_rt, -0.3625_rt, -0.075_rt},
                {-0.3375_rt, -0.2125_rt, 0.075_rt}),
        0, {0, 0, 0});
    BoxArray const nodal_grid(amrex::surroundingNodes(domain));
    DistributionMapping const dm(nodal_grid);
    MultiFab levelset(nodal_grid, dm, 1, 1);

    STLtools stl;
    stl.setUseMarchingCubes(true);
    stl.read_stl_file("cube.stl", 0.05_rt, {-0.4_rt, -0.275_rt, 0.0_rt}, 0);
    stl.fillMarchingCubesLevelSet(levelset, levelset.nGrowVect(), geom);

    Gpu::DeviceScalar<Long> error_count(static_cast<Long>(0));
    Long* const errors = error_count.dataPtr();
    for (MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        auto const phi = levelset.const_array(mfi);
        Box const face_nodes(IntVect(2, 2, 2), IntVect(2, 7, 6));
        ParallelFor(face_nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (phi(i, j, k) != 0.0_rt) {
                Gpu::Atomic::AddNoRet(errors, static_cast<Long>(1));
            }
        });
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        error_count.dataValue() == 0,
        "STL samples on a grid-aligned face were not consistently classified as on-surface");
}

void validate_narrow_band_levelset ()
{
    Box const domain(IntVect(0), IntVect(15));
    Geometry const geom(domain, RealBox({-1.0_rt, -1.0_rt, -1.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt}), 0,
                        {0, 0, 0});
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    MultiFab exact(nba, dm, 1, 1);
    MultiFab narrow(nba, dm, 1, 1);

    STLtools stl;
    stl.setUseMarchingCubes(true);
    stl.read_stl_file("cube.stl", 0.43_rt, {0.037_rt, -0.021_rt, 0.015_rt}, 0);
    stl.fillSignedDistance(exact, exact.nGrowVect(), geom);
    exact.OverrideSync(geom.periodicity());
    exact.FillBoundary(geom.periodicity());
    stl.fillMarchingCubesLevelSet(narrow, narrow.nGrowVect(), geom);

    Gpu::Buffer<Long> errors(5);
    std::fill_n(errors.hostData(), errors.size(), static_cast<Long>(0));
    errors.copyToDeviceAsync();
    Long* const error = errors.data();
    for (MFIter mfi(exact); mfi.isValid(); ++mfi) {
        auto const full = exact.const_array(mfi);
        auto const band = narrow.const_array(mfi);
        MC::MCFab mc_fab;
        mc_fab.defineEdgeIntersections(narrow[mfi].box());
        stl.fillMarchingCubesEdgeIntersections(narrow[mfi], geom, mc_fab);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            auto const crossing = mc_fab.m_edge_intersections[idim].const_array();
            ParallelFor(Box{crossing}, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int const di = idim == 0;
                int const dj = idim == 1;
                int const dk = idim == 2;
                Real const lo = band(i, j, k);
                Real const hi = band(i + di, j + dj, k + dk);
                Real const refined = crossing(i, j, k);
                // Edges without an exact STL crossing keep the sentinel; a
                // sign-changing edge inside the band must always have one.
                if ((lo > 0.0_rt) != (hi > 0.0_rt)) {
                    if (refined == MC::invalid_edge_intersection) {
                        Gpu::Atomic::AddNoRet(error + 4, static_cast<Long>(1));
                        return;
                    }
                    Gpu::Atomic::AddNoRet(error + 2, static_cast<Long>(1));
                    Real const linear = lo / (lo - hi);
                    if (std::abs(refined - linear) >
                        64.0_rt * std::numeric_limits<Real>::epsilon()) {
                        Gpu::Atomic::AddNoRet(error + 3, static_cast<Long>(1));
                    }
                }
            });
        }
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if ((full(i, j, k) > 0.0_rt) != (band(i, j, k) > 0.0_rt)) {
                Gpu::Atomic::AddNoRet(error, static_cast<Long>(1));
            }
        });

        Box const cbx = amrex::enclosedCells(mfi.validbox());
        ParallelFor(cbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            bool has_fluid = false;
            bool has_covered = false;
            for (int kk = 0; kk <= 1; ++kk) {
                for (int jj = 0; jj <= 1; ++jj) {
                    for (int ii = 0; ii <= 1; ++ii) {
                        bool const fluid = full(i + ii, j + jj, k + kk) > 0.0_rt;
                        has_fluid = has_fluid || fluid;
                        has_covered = has_covered || !fluid;
                    }
                }
            }
            if (has_fluid && has_covered) {
                for (int kk = 0; kk <= 1; ++kk) {
                    for (int jj = 0; jj <= 1; ++jj) {
                        for (int ii = 0; ii <= 1; ++ii) {
                            Real const a = full(i + ii, j + jj, k + kk);
                            Real const b = band(i + ii, j + jj, k + kk);
                            Real const tolerance = 16.0_rt * std::numeric_limits<Real>::epsilon() *
                                                   amrex::max(1.0_rt, std::abs(a));
                            if (std::abs(a - b) > tolerance) {
                                Gpu::Atomic::AddNoRet(error + 1, static_cast<Long>(1));
                            }
                        }
                    }
                }
            }
        });
    }
    errors.copyToHost();
    GpuArray<Long, 5> global_errors{errors.hostData()[0], errors.hostData()[1],
                                    errors.hostData()[2], errors.hostData()[3],
                                    errors.hostData()[4]};
    ParallelAllReduce::Sum(global_errors.data(), static_cast<int>(global_errors.size()),
                           ParallelContext::CommunicatorSub());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_errors[0] == 0, "Narrow-band MC field changed the global STL sign classification");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(global_errors[1] == 0,
                                     "Narrow-band MC field is not exact on a "
                                     "cut-cell interpolation stencil");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(global_errors[2] > 0,
                                     "Exact STL refinement found no sign-changing Cartesian edges");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_errors[3] > 0, "Exact STL refinement did not move any sampled-SDF edge crossing");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_errors[4] == 0, "A sign-changing edge in the band has no exact STL crossing");
}

char const* execution_backend () noexcept
{
#if defined(AMREX_USE_CUDA)
    return "CUDA";
#elif defined(AMREX_USE_HIP)
    return "HIP";
#elif defined(AMREX_USE_SYCL)
    return "SYCL";
#else
    return "CPU";
#endif
}

} // namespace

void main_main ()
{
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        Gpu::inLaunchRegion(), "The marching-cubes GPU test must execute in a GPU launch region");
#endif
    bool algorithm_tests = true;
    ParmParse pp;
    pp.query("algorithm_tests", algorithm_tests);
    if (algorithm_tests) {
        validate_mc33_face_decisions();
        validate_mc33_vertex_indices();
        validate_cell_topology_rejection();
        validate_exact_ambiguous_face_fraction();
        validate_nodal_stl_face_levelset();
    }
    bool scale_invariance_test = false;
    pp.query("scale_invariance_test", scale_invariance_test);

    int nx = 64;
    int ny = 64;
    int nz = 64;
    int max_grid_size = 32;
    int required_coarsening_level = 0;
    int max_coarsening_level = 0;
    bool build_coarse_level_by_coarsening = true;
    Real xmin = -1.2;
    Real xmax = 1.2;
    Real ymin = -1.2;
    Real ymax = 1.2;
    Real zmin = -1.2;
    Real zmax = 1.2;
    bool custom_stl_test = false;
    bool narrow_band_test = false;
    std::string eb_method("marching_cubes");
    std::string geom_type("stl");
    std::string api_build;
    bool api_all_levels = false;
    bool allow_full_volume_cut_cells = false;
    std::vector<int> is_periodic{0, 0, 0};
    std::string stl_file("cube.stl");
    Real stl_scale = 1.0;
    std::vector<Real> stl_center{0.0, 0.0, 0.0};
    Real body_volume = -1.0_rt; // analytic solid volume, < 0 when unknown
    bool fluid_inside = false;
    {
        pp.query("nx", nx);
        pp.query("ny", ny);
        pp.query("nz", nz);
        pp.query("max_grid_size", max_grid_size);
        pp.query("xmin", xmin);
        pp.query("xmax", xmax);
        pp.query("ymin", ymin);
        pp.query("ymax", ymax);
        pp.query("zmin", zmin);
        pp.query("zmax", zmax);
        pp.query("required_coarsening_level", required_coarsening_level);
        pp.query("max_coarsening_level", max_coarsening_level);
        pp.query("build_coarse_level_by_coarsening", build_coarse_level_by_coarsening);
        pp.query("custom_stl_test", custom_stl_test);
        pp.query("narrow_band_test", narrow_band_test);
        // Build through the C++ API with an explicit GeometryShop instead of
        // the ParmParse-driven EB2::Build: "sphere" (GPU-callable SphereIF) or
        // "host_sphere" (host-only functor).
        pp.query("api_build", api_build);
        // With api_build, also build every coarse level directly from the
        // GeometryShop (EB2::Build overload taking a Vector<Geometry>) instead
        // of coarsening the fine level.
        pp.query("api_all_levels", api_all_levels);
        std::vector<Real> prob_lo{xmin, ymin, zmin};
        std::vector<Real> prob_hi{xmax, ymax, zmax};
        ParmParse ppgeom("geometry");
        if (ppgeom.queryarr("prob_lo", prob_lo)) {
            xmin = prob_lo[0];
            ymin = prob_lo[1];
            zmin = prob_lo[2];
        }
        if (ppgeom.queryarr("prob_hi", prob_hi)) {
            xmax = prob_hi[0];
            ymax = prob_hi[1];
            zmax = prob_hi[2];
        }
        ppgeom.queryarr("is_periodic", is_periodic);
        pp.query("allow_full_volume_cut_cells", allow_full_volume_cut_cells);

        ParmParse ppeb2("eb2");
        ppeb2.queryAdd("geom_type", geom_type);
        ppeb2.queryAdd("geometry_method", eb_method);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            eb_method == "legacy" || eb_method == "marching_cubes",
            "eb2.geometry_method must be legacy or marching_cubes");

        // Analytic solid volume for the geometry types the test knows.
        if (!api_build.empty() || geom_type == "sphere") {
            Real radius = 0.5_rt;
            ppeb2.queryAdd("sphere_radius", radius);
            std::vector<Real> center{0.0_rt, 0.0_rt, 0.0_rt};
            ppeb2.queryAdd("sphere_center", center);
            int has_fluid_inside = 0;
            ppeb2.queryAdd("sphere_has_fluid_inside", has_fluid_inside);
            fluid_inside = has_fluid_inside != 0;
            body_volume = 4.0_rt/3.0_rt*std::acos(-1.0_rt)*radius*radius*radius;
        } else if (geom_type == "stl") {
            ppeb2.queryAdd("stl_file", stl_file);
            ppeb2.queryAdd("stl_scale", stl_scale);
            ppeb2.queryAdd("stl_center", stl_center);
            // cube.stl is the [-1,1]^3 cube; other STL files use custom_stl_test.
            body_volume = 8.0_rt*stl_scale*stl_scale*stl_scale;
        } else if (geom_type == "box") {
            std::vector<Real> lo, hi;
            ppeb2.getarr("box_lo", lo);
            ppeb2.getarr("box_hi", hi);
            int has_fluid_inside = 0;
            ppeb2.queryAdd("box_has_fluid_inside", has_fluid_inside);
            fluid_inside = has_fluid_inside != 0;
            // Only the part of the box inside the domain contributes.
            Real const clipped[3] = {
                amrex::max(0.0_rt, amrex::min(hi[0], xmax) - amrex::max(lo[0], xmin)),
                amrex::max(0.0_rt, amrex::min(hi[1], ymax) - amrex::max(lo[1], ymin)),
                amrex::max(0.0_rt, amrex::min(hi[2], zmax) - amrex::max(lo[2], zmin))};
            body_volume = clipped[0]*clipped[1]*clipped[2];
        } else if (geom_type == "cylinder") {
            Real radius = 0.5_rt;
            Real height = -1.0_rt;
            ppeb2.get("cylinder_radius", radius);
            ppeb2.queryAdd("cylinder_height", height);
            int has_fluid_inside = 0;
            ppeb2.queryAdd("cylinder_has_fluid_inside", has_fluid_inside);
            fluid_inside = has_fluid_inside != 0;
            if (height > 0.0_rt) {
                body_volume = std::acos(-1.0_rt)*radius*radius*height;
            }
        }
        if (body_volume < 0.0_rt) {
            // parser, torus, plane, ...: only range and consistency checks.
            custom_stl_test = true;
        }
    }

    if (narrow_band_test) {
        validate_narrow_band_levelset();
    }
    if (scale_invariance_test) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(eb_method == "marching_cubes",
                                         "scale_invariance_test requires marching_cubes");
        validate_scale_invariance();
    }

    Geometry geom(Box(IntVect(0), IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1))),
                  RealBox({AMREX_D_DECL(xmin, ymin, zmin)}, {AMREX_D_DECL(xmax, ymax, zmax)}), 0,
                  {AMREX_D_DECL(is_periodic[0], is_periodic[1], is_periodic[2])});
    BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    double t0 = amrex::second();
    if (api_build.empty()) {
        EB2::Build(geom, required_coarsening_level, max_coarsening_level, 4,
                   build_coarse_level_by_coarsening);
    } else {
        Real radius = 0.5_rt;
        std::vector<Real> center{0.0_rt, 0.0_rt, 0.0_rt};
        ParmParse ppeb2("eb2");
        ppeb2.queryAdd("sphere_radius", radius);
        ppeb2.queryAdd("sphere_center", center);
        RealArray const c{AMREX_D_DECL(center[0], center[1], center[2])};
        Vector<Geometry> all_geoms;
        if (api_all_levels) {
            all_geoms.push_back(geom);
            for (int ilev = 1; ilev <= max_coarsening_level; ++ilev) {
                all_geoms.push_back(amrex::coarsen(all_geoms.back(), 2));
            }
        }
        auto build_with = [&] (auto const& gshop) {
            if (api_all_levels) {
                EB2::Build(gshop, all_geoms, 4);
            } else {
                EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level, 4,
                           build_coarse_level_by_coarsening);
            }
        };
        if (api_build == "sphere") {
            build_with(EB2::makeShop(EB2::SphereIF(radius, c, fluid_inside)));
        } else if (api_build == "host_sphere") {
            static_assert(!EB2::IsGPUable<HostSphereIF>::value,
                          "HostSphereIF must exercise the host-only adapter path");
            build_with(EB2::makeShop(HostSphereIF{.m_radius = radius, .m_center = c,
                                                  .m_fluid_inside = fluid_inside}));
        } else {
            amrex::Abort("api_build must be sphere or host_sphere");
        }
    }
    double t1 = amrex::second();
    amrex::Print() << "EB method: " << eb_method << ", build time: " << t1 - t0 << "\n";
    auto factory = makeEBFabFactory(geom, ba, dm, {1, 1, 1}, EBSupport::full);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(factory->maxCoarseningLevel() >= required_coarsening_level,
                                     "Marching-cubes EB did not build every required coarse level");
    std::string eb_surface_stl_file;
    if (ParmParse("eb2").query("eb_surface_stl_file", eb_surface_stl_file)) {
        WriteEBSurfaceSTL(ba, dm, geom, factory.get(), eb_surface_stl_file);
    }
    std::string mc_stl_file;
    if (ParmParse("eb2").query("mc_stl_file", mc_stl_file) && eb_method == "marching_cubes") {
        // The file holds the converged triangulation of the finest level.
        Long facets = 0;
        if (ParallelDescriptor::IOProcessor()) {
            std::ifstream ifs(mc_stl_file);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ifs.good(), "eb2.mc_stl_file was not written");
            std::string line;
            while (std::getline(ifs, line)) {
                if (line.find("facet normal") != std::string::npos) { ++facets; }
            }
        }
        ParallelDescriptor::Bcast(&facets, 1);
        amrex::Print() << "Marching-cubes STL facets: " << facets << "\n";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(facets > 0, "eb2.mc_stl_file contains no facets");
        Long expected_facets = -1;
        pp.query("expected_mc_stl_facets", expected_facets);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(expected_facets < 0 || facets == expected_facets,
                                         "eb2.mc_stl_file facet count differs from the "
                                         "expected finest-level count");
    }
    auto const& volfrac = factory->getVolFrac();
    Real const fluid_volume =
        volfrac.sum() * geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    Real const domain_volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    Real const expected_volume = fluid_inside ? body_volume : domain_volume - body_volume;
    Real const error = std::abs(fluid_volume - expected_volume);
    Real const max_dx =
        amrex::max(geom.CellSize(0), amrex::max(geom.CellSize(1), geom.CellSize(2)));
    Real const roundoff_tolerance = 64.0_rt * std::numeric_limits<Real>::epsilon();
    if (custom_stl_test) {
        amrex::Print() << "Fluid volume: " << fluid_volume << "\n";
    } else {
        amrex::Print() << "Fluid volume: " << fluid_volume << ", expected: " << expected_volume
                       << ", error: " << error << "\n";
    }
    AMREX_ALWAYS_ASSERT(volfrac.min(0) >= -roundoff_tolerance);
    AMREX_ALWAYS_ASSERT(volfrac.max(0) <= 1.0_rt + roundoff_tolerance);
    if (custom_stl_test) {
        AMREX_ALWAYS_ASSERT(fluid_volume > 0.0_rt);
        AMREX_ALWAYS_ASSERT(fluid_volume < domain_volume);
    }

    auto volume_centroid = factory->getCentroid().ToMultiFab(0.0_rt, 0.0_rt);
    auto boundary_centroid = factory->getBndryCent().ToMultiFab(0.0_rt, 0.0_rt);
    auto boundary_normal = factory->getBndryNormal().ToMultiFab(0.0_rt, 0.0_rt);
    auto bndry_area = factory->getBndryArea().ToMultiFab(0.0_rt, 0.0_rt);
    AMREX_ALWAYS_ASSERT(!volume_centroid.contains_nan());
    AMREX_ALWAYS_ASSERT(!boundary_centroid.contains_nan());
    AMREX_ALWAYS_ASSERT(!boundary_normal.contains_nan());
    AMREX_ALWAYS_ASSERT(!bndry_area.contains_nan());
    AMREX_ALWAYS_ASSERT(!volume_centroid.contains_inf());
    AMREX_ALWAYS_ASSERT(!boundary_centroid.contains_inf());
    AMREX_ALWAYS_ASSERT(!boundary_normal.contains_inf());
    AMREX_ALWAYS_ASSERT(!bndry_area.contains_inf());
    AMREX_ALWAYS_ASSERT(bndry_area.min(0) >= -roundoff_tolerance);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        AMREX_ALWAYS_ASSERT(volume_centroid.min(d) >= -0.5_rt - roundoff_tolerance);
        AMREX_ALWAYS_ASSERT(volume_centroid.max(d) <= 0.5_rt + roundoff_tolerance);
        // Non-cut entries in a mixed CutFab retain the EB sentinel -1.
        AMREX_ALWAYS_ASSERT(boundary_centroid.min(d) >= -1.0_rt - roundoff_tolerance);
        AMREX_ALWAYS_ASSERT(boundary_centroid.max(d) <= 0.5_rt + roundoff_tolerance);
        AMREX_ALWAYS_ASSERT(boundary_normal.min(d) >= -1.0_rt - roundoff_tolerance);
        AMREX_ALWAYS_ASSERT(boundary_normal.max(d) <= 1.0_rt + roundoff_tolerance);
    }

    auto const area_fraction = factory->getAreaFrac();
    auto const face_centroid = factory->getFaceCent();
    Array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> dense_area;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dense_area[d] = std::make_unique<MultiFab>(area_fraction[d]->ToMultiFab(1.0_rt, 0.0_rt));
        auto face_cent = face_centroid[d]->ToMultiFab(0.0_rt, 0.0_rt);
        AMREX_ALWAYS_ASSERT(!dense_area[d]->contains_nan());
        AMREX_ALWAYS_ASSERT(!face_cent.contains_nan());
        AMREX_ALWAYS_ASSERT(!dense_area[d]->contains_inf());
        AMREX_ALWAYS_ASSERT(!face_cent.contains_inf());
        AMREX_ALWAYS_ASSERT(dense_area[d]->min(0) >= -roundoff_tolerance);
        AMREX_ALWAYS_ASSERT(dense_area[d]->max(0) <= 1.0_rt + roundoff_tolerance);
        for (int n = 0; n < AMREX_SPACEDIM - 1; ++n) {
            AMREX_ALWAYS_ASSERT(face_cent.min(n) >= -0.5_rt - roundoff_tolerance);
            AMREX_ALWAYS_ASSERT(face_cent.max(n) <= 0.5_rt + roundoff_tolerance);
        }
    }

    auto const& flags = factory->getMultiEBCellFlagFab();
    auto const& levelset = factory->getLevelSet();
    AMREX_ALWAYS_ASSERT(!levelset.contains_nan());
    AMREX_ALWAYS_ASSERT(!levelset.contains_inf());

    auto const edge_centroid = factory->getEdgeCent();
    Array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> dense_edge;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dense_edge[d] = std::make_unique<MultiFab>(edge_centroid[d]->ToMultiFab(1.0_rt, -1.0_rt));
        AMREX_ALWAYS_ASSERT(!dense_edge[d]->contains_nan());
        AMREX_ALWAYS_ASSERT(!dense_edge[d]->contains_inf());
        AMREX_ALWAYS_ASSERT(dense_edge[d]->min(0) >= -1.0_rt - roundoff_tolerance);
        AMREX_ALWAYS_ASSERT(dense_edge[d]->max(0) <= 1.0_rt + roundoff_tolerance);
    }

    Gpu::Buffer<int> levelset_errors(2);
    std::fill_n(levelset_errors.hostData(), levelset_errors.size(), 0);
    levelset_errors.copyToDeviceAsync();
    int* const level_errors = levelset_errors.data();
    for (MFIter mfi(flags); mfi.isValid(); ++mfi) {
        auto const flag = flags.const_array(mfi);
        auto const phi = levelset.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (flag(i, j, k).isCovered()) {
                for (int kk = 0; kk <= 1; ++kk) {
                    for (int jj = 0; jj <= 1; ++jj) {
                        for (int ii = 0; ii <= 1; ++ii) {
                            if (phi(i + ii, j + jj, k + kk) < 0.0_rt) {
                                Gpu::Atomic::AddNoRet(level_errors, 1);
                            }
                        }
                    }
                }
            }
        });
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (MFIter mfi(*dense_edge[d]); mfi.isValid(); ++mfi) {
            auto const edge = dense_edge[d]->const_array(mfi);
            auto const phi = levelset.const_array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int const ih = i + (d == 0);
                int const jh = j + (d == 1);
                int const kh = k + (d == 2);
                bool const lo_fluid = phi(i, j, k) < 0.0_rt;
                bool const hi_fluid = phi(ih, jh, kh) < 0.0_rt;
                bool mismatch = false;
                if (lo_fluid && hi_fluid) {
                    mismatch = edge(i, j, k) != 1.0_rt;
                } else if (!lo_fluid && !hi_fluid) {
                    mismatch = edge(i, j, k) != -1.0_rt;
                } else {
                    mismatch = edge(i, j, k) < -0.5_rt || edge(i, j, k) > 0.5_rt;
                }
                if (mismatch) {
                    Gpu::Atomic::AddNoRet(level_errors + 1, 1);
                }
            });
        }
    }
    levelset_errors.copyToHost();
    GpuArray<int, 2> global_levelset_errors{levelset_errors.hostData()[0],
                                            levelset_errors.hostData()[1]};
    ParallelAllReduce::Sum(global_levelset_errors.data(),
                           static_cast<int>(global_levelset_errors.size()),
                           ParallelContext::CommunicatorSub());
    // The legacy generator repairs the level set independently in every FAB
    // and the exported level set takes an arbitrary owner on shared nodes.
    // Without deterministic floating point (-ffast-math) two FABs can round
    // the same borderline cell differently, so cross-FAB agreement of the
    // exported signs is only guaranteed for the marching-cubes generator.
#ifdef __FAST_MATH__
    bool const check_levelset_signs = (eb_method == "marching_cubes");
#else
    bool const check_levelset_signs = true;
#endif
    if (check_levelset_signs) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            global_levelset_errors[0] == 0,
            "A covered cell retains a public negative-in-fluid level-set node");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            global_levelset_errors[1] == 0,
            "Edge centroids disagree with repaired public level-set signs");
    } else {
        amrex::Print() << "Skipping level-set sign consistency checks for the "
                       << "legacy generator in a fast-math build ("
                       << global_levelset_errors[0] << " covered-cell, "
                       << global_levelset_errors[1] << " edge mismatches)\n";
    }

    Gpu::DeviceScalar<Long> zero_node_count(static_cast<Long>(0));
    Long* const zero_nodes = zero_node_count.dataPtr();
    for (MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        auto const phi = levelset.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (phi(i, j, k) == 0.0_rt) {
                Gpu::Atomic::AddNoRet(zero_nodes, static_cast<Long>(1));
            }
        });
    }
    Long repaired_nodes = zero_node_count.dataValue();
    ParallelAllReduce::Sum(repaired_nodes, ParallelContext::CommunicatorSub());

    Real small_volfrac = 0.0_rt;
    ParmParse("eb2").query("small_volfrac", small_volfrac);
    Gpu::Buffer<int> topology_counts(6);
    std::fill_n(topology_counts.hostData(), topology_counts.size(), 0);
    topology_counts.copyToDeviceAsync();
    int* const topology = topology_counts.data();
    Box const domain = geom.Domain();
    for (MFIter mfi(flags); mfi.isValid(); ++mfi) {
        Box const bx = mfi.validbox();
        auto const flag = flags.const_array(mfi);
        auto const apx = dense_area[0]->const_array(mfi);
        auto const apy = dense_area[1]->const_array(mfi);
        auto const apz = dense_area[2]->const_array(mfi);
        auto const vf = volfrac.const_array(mfi);
        auto const bndry_area_arr = bndry_area.const_array(mfi);
        auto const bn = boundary_normal.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (flag(i, j, k).isSingleValued()) {
                Gpu::Atomic::AddNoRet(topology + 1, 1);
                if (vf(i, j, k) + roundoff_tolerance < small_volfrac) {
                    Gpu::Atomic::AddNoRet(topology + 2, 1);
                }
                if (vf(i, j, k) == 1.0_rt && bndry_area_arr(i, j, k) > 0.0_rt) {
                    Gpu::Atomic::AddNoRet(topology + 4, 1);
                }
                Real const area_x = apx(i, j, k) - apx(i + 1, j, k);
                Real const area_y = apy(i, j, k) - apy(i, j + 1, k);
                Real const area_z = apz(i, j, k) - apz(i, j, k + 1);
                Real const norm =
                    std::sqrt(area_x * area_x + area_y * area_y + area_z * area_z);
                Real const normal_tolerance =
                    amrex::max(5.e-12_rt, roundoff_tolerance * amrex::max(1.0_rt, norm));
                if (norm <= 1.e-12_rt ||
                    std::abs(norm * bn(i, j, k, 0) - area_x) > normal_tolerance ||
                    std::abs(norm * bn(i, j, k, 1) - area_y) > normal_tolerance ||
                    std::abs(norm * bn(i, j, k, 2) - area_z) > normal_tolerance) {
                    Gpu::Atomic::AddNoRet(topology + 5, 1);
                }
            }
            if (flag(i, j, k).isMultiValued()) {
                Gpu::Atomic::AddNoRet(topology + 3, 1);
            }

            int mismatches = 0;
            if (i > domain.smallEnd(0)) {
                mismatches += flag(i, j, k).isConnected(-1, 0, 0) != (apx(i, j, k) > 0.0_rt);
            }
            if (i < domain.bigEnd(0)) {
                mismatches += flag(i, j, k).isConnected(1, 0, 0) != (apx(i + 1, j, k) > 0.0_rt);
            }
            if (j > domain.smallEnd(1)) {
                mismatches += flag(i, j, k).isConnected(0, -1, 0) != (apy(i, j, k) > 0.0_rt);
            }
            if (j < domain.bigEnd(1)) {
                mismatches += flag(i, j, k).isConnected(0, 1, 0) != (apy(i, j + 1, k) > 0.0_rt);
            }
            if (k > domain.smallEnd(2)) {
                mismatches += flag(i, j, k).isConnected(0, 0, -1) != (apz(i, j, k) > 0.0_rt);
            }
            if (k < domain.bigEnd(2)) {
                mismatches += flag(i, j, k).isConnected(0, 0, 1) != (apz(i, j, k + 1) > 0.0_rt);
            }
            if (mismatches != 0) {
                Gpu::Atomic::AddNoRet(topology, mismatches);
            }
        });
    }
    topology_counts.copyToHost();
    GpuArray<int, 6> global_topology_counts;
    std::copy_n(topology_counts.hostData(), global_topology_counts.size(),
                global_topology_counts.begin());
    ParallelAllReduce::Sum(global_topology_counts.data(),
                           static_cast<int>(global_topology_counts.size()),
                           ParallelContext::CommunicatorSub());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[0] == 0,
        "EBCellFlag coordinate connectivity disagrees with MC face apertures");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(global_topology_counts[1] > 0,
                                     "Marching-cubes test did not produce any single-valued cells");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[2] == 0,
        "Nodal repair left a single-valued cell below eb2.small_volfrac");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(global_topology_counts[3] == 0,
                                     "Single-valued EB regression retained a multi-valued cell");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[5] == 0,
        "EB boundary normals disagree with repaired face-aperture closure");
    // A single-valued cell with unit volume and a nonzero boundary area owns
    // an EB patch coincident with one of its faces or nodes.  This state only
    // arises when the surface passes exactly through grid nodes, either
    // because the geometry does (allow_full_volume_cut_cells = 1) or because
    // the nodal repair moved nodes onto the surface.
    if (!allow_full_volume_cut_cells && repaired_nodes == 0) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            global_topology_counts[4] == 0,
            "A single-valued cell has unit volume and a nonzero boundary area");
    }
    amrex::Print() << "Single-valued cells: " << global_topology_counts[1]
                   << ", retained multi-valued cells: " << global_topology_counts[3]
                   << ", full-volume cut cells: " << global_topology_counts[4] << "\n";
    amrex::Print() << std::setprecision(std::numeric_limits<Real>::max_digits10)
                   << "MC_TEST_SIGNATURE backend=" << execution_backend()
                   << " precision=" << 8 * sizeof(Real) << " method=" << eb_method
                   << " fluid_volume=" << fluid_volume
                   << " single_valued_cells=" << global_topology_counts[1]
                   << " repaired_nodes=" << repaired_nodes << '\n';
    if (!custom_stl_test) {
        Real const domain_scale = amrex::max(xmax - xmin, amrex::max(ymax - ymin, zmax - zmin));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            error < 6.0_rt * max_dx * max_dx * domain_scale,
            "Marching-cubes cut-cell volume error exceeds the expected O(dx^2) bound");
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    main_main();
    amrex::Finalize();
}
