
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_STL_utils.H>
#include <AMReX_MarchingCubes.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WriteEBSurface.H>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

using namespace amrex;

namespace {

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
    MC::marching_cubes(geom, sdf, result);

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
        MC::marching_cubes(geom, sdf, result);

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
    int const errors =
        MC::build_face_fractions(cell_box, mc_fab, sdf, apx, apy, apz, fcx, fcy, fcz);

    Real area = 0.0_rt;
    Gpu::dtoh_memcpy(&area, apz.dataPtr(), sizeof(Real));
    Gpu::streamSynchronize();
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

void validate_domain_extension_is_feature_local ()
{
    Box const cells(IntVect(0), IntVect(2));
    Box const extended_nodes(IntVect(-2, -1, -1), IntVect(5, 4, 4));
    FArrayBox levelset(extended_nodes, 1);
    levelset.setVal<RunOn::Device>(1.0_rt);
    auto const phi = levelset.array();
    ParallelFor(Box(IntVect(0, 1, 1), IntVect(0, 2, 2)),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    phi(i, j, k) = -1.0_rt;
                });

    MC::extend_domain_face_levelset(
        extended_nodes, cells, GpuArray<int, 3>{0, 0, 0}, levelset);
    auto const extended_phi = levelset.const_array();
    Gpu::DeviceScalar<Long> error_count(static_cast<Long>(0));
    Long* const errors = error_count.dataPtr();
    ParallelFor(extended_nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bool const normal_extrusion = i < 0 && j >= 1 && j <= 2 && k >= 1 && k <= 2;
        bool const source_feature = i == 0 && j >= 1 && j <= 2 && k >= 1 && k <= 2;
        bool const covered = extended_phi(i, j, k) < 0.0_rt;
        if (covered != (normal_extrusion || source_feature)) {
            Gpu::Atomic::AddNoRet(errors, static_cast<Long>(1));
        }
    });
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        error_count.dataValue() == 0,
        "MC domain extension covered nodes outside the selected exterior feature");

    MC::MCFab mc_fab;
    mc_fab.defineEdgeIntersections(extended_nodes);
    auto const transverse_crossing = mc_fab.m_edge_intersections[1].array();
    ParallelFor(Box(IntVect(0, 1, 1), IntVect(0, 1, 1)),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    transverse_crossing(i, j, k) = 0.25_rt;
                });
    MC::extend_domain_face_edge_intersections(
        cells, GpuArray<int, 3>{0, 0, 0}, mc_fab);
    auto const normal_crossing = mc_fab.m_edge_intersections[0].const_array();
    Gpu::DeviceScalar<Long> cache_error_count(static_cast<Long>(0));
    Long* const cache_errors = cache_error_count.dataPtr();
    ParallelFor(Box(IntVect(-2, 1, 1), IntVect(-1, 1, 1)),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (transverse_crossing(i, j, k) != 0.25_rt) {
                        Gpu::Atomic::AddNoRet(cache_errors, static_cast<Long>(1));
                    }
                });
    ParallelFor(Box(IntVect(-1, 1, 0), IntVect(-1, 1, 0)),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (!amrex::isnan(transverse_crossing(i, j, k))
                        || !amrex::isnan(normal_crossing(i, j, k + 1))) {
                        Gpu::Atomic::AddNoRet(cache_errors, static_cast<Long>(1));
                    }
                });
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        cache_error_count.dataValue() == 0,
        "MC exact crossings were not extended only along the selected face normal");
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

    Gpu::Buffer<Long> errors(4);
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
                if ((lo > 0.0_rt) != (hi > 0.0_rt) && !amrex::isnan(refined)) {
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
    GpuArray<Long, 4> global_errors{errors.hostData()[0], errors.hostData()[1],
                                    errors.hostData()[2], errors.hostData()[3]};
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

Long validate_ascii_stl (std::string const& filename)
{
    if (filename.empty()) {
        return 0;
    }
    ParallelDescriptor::Barrier();
    Long facets = 0;
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream input(filename);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(input.good(), "Could not open generated STL " + filename);
        std::string line;
        Long headers = 0;
        Long trailers = 0;
        Long vertex_lines = 0;
        GpuArray<GpuArray<Real, 3>, 3> vertices;
        int vertex_count = 0;
        while (std::getline(input, line)) {
            facets += line.starts_with("facet normal ");
            headers += line.starts_with("solid Created by AMReX");
            trailers += line.starts_with("endsolid Created by AMReX");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                line.find("nan") == std::string::npos && line.find("inf") == std::string::npos,
                "Generated STL contains a non-finite value: " + filename);
            if (line.find("vertex") != std::string::npos) {
                std::istringstream parser(line);
                std::string tag;
                parser >> tag >> vertices[vertex_count][0] >> vertices[vertex_count][1] >>
                    vertices[vertex_count][2];
                ++vertex_lines;
                ++vertex_count;
                if (vertex_count == 3) {
                    Real const ax = vertices[1][0] - vertices[0][0];
                    Real const ay = vertices[1][1] - vertices[0][1];
                    Real const az = vertices[1][2] - vertices[0][2];
                    Real const bx = vertices[2][0] - vertices[0][0];
                    Real const by = vertices[2][1] - vertices[0][1];
                    Real const bz = vertices[2][2] - vertices[0][2];
                    Real const cx = ay * bz - az * by;
                    Real const cy = az * bx - ax * bz;
                    Real const cz = ax * by - ay * bx;
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cx * cx + cy * cy + cz * cz > 0.0_rt,
                                                     "Generated STL contains a degenerate facet: " +
                                                         filename);
                    vertex_count = 0;
                }
            }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            facets > 0 && headers == 1 && trailers == 1 && vertex_count == 0 &&
                vertex_lines == 3 * facets,
            "Generated STL has incomplete facets or invalid solid records: " + filename);
    }
    ParallelDescriptor::Barrier();
    ParallelDescriptor::Bcast(&facets, 1, ParallelDescriptor::IOProcessorNumber());
    return facets;
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
        validate_exact_ambiguous_face_fraction();
        validate_nodal_stl_face_levelset();
        validate_domain_extension_is_feature_local();
    }

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
    std::vector<int> is_periodic{0, 0, 0};
    bool cleanup_test = false;
    bool custom_stl_test = false;
    bool narrow_band_test = false;
    bool extend_domain_test = false;
    std::string eb_method("marching_cubes");
    std::string stl_file("cube.stl");
    Real stl_scale = 1.0;
    std::vector<Real> stl_center{0.0, 0.0, 0.0};
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
        pp.queryarr("is_periodic", is_periodic);
        pp.query("required_coarsening_level", required_coarsening_level);
        pp.query("max_coarsening_level", max_coarsening_level);
        pp.query("build_coarse_level_by_coarsening", build_coarse_level_by_coarsening);
        pp.query("cleanup_test", cleanup_test);
        pp.query("custom_stl_test", custom_stl_test);
        pp.query("narrow_band_test", narrow_band_test);
        pp.query("extend_domain_test", extend_domain_test);
        pp.query("test_stl_file", stl_file);
        pp.query("test_stl_scale", stl_scale);
        pp.queryarr("test_stl_center", stl_center);

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

        ParmParse ppeb2("eb2");
        std::string geom_type("stl");
        ppeb2.queryAdd("geom_type", geom_type);
        ppeb2.queryAdd("stl_file", stl_file);
        ppeb2.queryAdd("stl_scale", stl_scale);
        ppeb2.queryAdd("stl_center", stl_center);
        ppeb2.queryAdd("stl_geometry_method", eb_method);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            eb_method == "legacy" || eb_method == "marching_cubes",
            "eb2.stl_geometry_method must be legacy or marching_cubes");
    }

    if (narrow_band_test) {
        validate_narrow_band_levelset();
    }

    Geometry geom(Box(IntVect(0), IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1))),
                  RealBox({AMREX_D_DECL(xmin, ymin, zmin)}, {AMREX_D_DECL(xmax, ymax, zmax)}), 0,
                  {AMREX_D_DECL(is_periodic[0], is_periodic[1], is_periodic[2])});
    BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    double t0 = amrex::second();
    EB2::Build(geom, required_coarsening_level, max_coarsening_level, 4,
               build_coarse_level_by_coarsening);
    double t1 = amrex::second();
    amrex::Print() << "EB method: " << eb_method << ", build time: " << t1 - t0 << "\n";
    auto factory = makeEBFabFactory(geom, ba, dm, {1, 1, 1}, EBSupport::full);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(factory->maxCoarseningLevel() >= required_coarsening_level,
                                     "Marching-cubes EB did not build every required coarse level");
    std::string eb_surface_stl_file;
    ParmParse("eb2").query("eb_surface_stl_file", eb_surface_stl_file);
    if (!eb_surface_stl_file.empty()) {
        WriteEBSurfaceSTL(ba, dm, geom, factory.get(), eb_surface_stl_file);
        validate_ascii_stl(eb_surface_stl_file);
    }
    std::string mc_stl_file;
    ParmParse("eb2").query("mc_stl_file", mc_stl_file);
    Long const mc_facets = eb_method == "marching_cubes" ? validate_ascii_stl(mc_stl_file) : 0;
    auto const& volfrac = factory->getVolFrac();
    Real const fluid_volume =
        volfrac.sum() * geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    Real const domain_volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    Real const expected_volume = domain_volume - 8.0_rt * stl_scale * stl_scale * stl_scale;
    Real const error = std::abs(fluid_volume - expected_volume);
    Real const max_dx =
        amrex::max(geom.CellSize(0), amrex::max(geom.CellSize(1), geom.CellSize(2)));
    Real const roundoff_tolerance = 64.0_rt * std::numeric_limits<Real>::epsilon();
    if (custom_stl_test || extend_domain_test) {
        amrex::Print() << "Fluid volume: " << fluid_volume << "\n";
    } else {
        amrex::Print() << "Fluid volume: " << fluid_volume << ", expected: " << expected_volume
                       << ", error: " << error << "\n";
    }
    AMREX_ALWAYS_ASSERT(volfrac.min(0) >= -roundoff_tolerance);
    AMREX_ALWAYS_ASSERT(volfrac.max(0) <= 1.0_rt + roundoff_tolerance);
    if (custom_stl_test || extend_domain_test) {
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

    if (extend_domain_test) {
        IntVect const covered_probe_iv(
            AMREX_D_DECL(geom.Domain().smallEnd(0) - 1,
                         (geom.Domain().smallEnd(1) + geom.Domain().bigEnd(1) + 1) / 2,
                         (geom.Domain().smallEnd(2) + geom.Domain().bigEnd(2) + 1) / 2));
        IntVect const fluid_probe_iv(
            AMREX_D_DECL(geom.Domain().smallEnd(0) - 1,
                         geom.Domain().smallEnd(1) + 25,
                         (geom.Domain().smallEnd(2) + geom.Domain().bigEnd(2) + 1) / 2));
        Gpu::Buffer<Long> extension_counts(4);
        std::fill_n(extension_counts.hostData(), extension_counts.size(), static_cast<Long>(0));
        extension_counts.copyToDeviceAsync();
        Long* const extension = extension_counts.data();
        for (MFIter mfi(levelset); mfi.isValid(); ++mfi) {
            Box const covered_probe = Box(covered_probe_iv, covered_probe_iv) & levelset[mfi].box();
            Box const fluid_probe = Box(fluid_probe_iv, fluid_probe_iv) & levelset[mfi].box();
            auto const phi = levelset.const_array(mfi);
            ParallelFor(covered_probe, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Gpu::Atomic::AddNoRet(extension, static_cast<Long>(1));
                if (phi(i, j, k) >= 0.0_rt) {
                    Gpu::Atomic::AddNoRet(extension + 1, static_cast<Long>(1));
                }
            });
            ParallelFor(fluid_probe, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Gpu::Atomic::AddNoRet(extension + 2, static_cast<Long>(1));
                if (phi(i, j, k) < 0.0_rt) {
                    Gpu::Atomic::AddNoRet(extension + 3, static_cast<Long>(1));
                }
            });
        }
        extension_counts.copyToHost();
        Long extension_total = extension_counts.hostData()[0];
        Long extension_covered = extension_counts.hostData()[1];
        Long fluid_total = extension_counts.hostData()[2];
        Long fluid_unchanged = extension_counts.hostData()[3];
        ParallelAllReduce::Sum<Long>({extension_total, extension_covered, fluid_total, fluid_unchanged},
                                     ParallelContext::CommunicatorSub());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            extension_total > 0 && extension_covered == extension_total,
            "MC domain-face extension did not keep the exterior probe covered");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fluid_total > 0 && fluid_unchanged == fluid_total,
            "MC domain-face extension covered an exterior location unrelated to the boundary feature");
    }

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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_levelset_errors[0] == 0,
        "A covered cell retains a public negative-in-fluid level-set node");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_levelset_errors[1] == 0,
        "Edge centroids disagree with repaired public level-set signs");

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
    if (cleanup_test) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            repaired_nodes > 0, "Legacy-style MC repair did not move any fluid node to zero");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[5] == 0,
        "EB boundary normals disagree with repaired face-aperture closure");
    amrex::Print() << "Single-valued cells: " << global_topology_counts[1]
                   << ", retained multi-valued cells: " << global_topology_counts[3] << "\n";
    amrex::Print() << std::setprecision(std::numeric_limits<Real>::max_digits10)
                   << "MC_TEST_SIGNATURE backend=" << execution_backend()
                   << " precision=" << 8 * sizeof(Real) << " method=" << eb_method
                   << " fluid_volume=" << fluid_volume
                   << " single_valued_cells=" << global_topology_counts[1]
                   << " repaired_nodes=" << repaired_nodes << " mc_facets=" << mc_facets << '\n';
    if (!cleanup_test && !custom_stl_test && !extend_domain_test) {
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
