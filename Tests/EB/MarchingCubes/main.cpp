
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_STL_utils.H>
#include <AMReX_MarchingCubes.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WriteEBSurface.H>

#include <fstream>
#include <sstream>

using namespace amrex;

namespace {

bool resolved_bottom_face_is_fluid_connected(GpuArray<Real, 8> const &values) {
  Box const cell_box(IntVect(0), IntVect(0));
  Box const node_box = amrex::surroundingNodes(amrex::grow(cell_box, 1));
  FArrayBox sdf(node_box, 1);
  sdf.setVal(-1.0_rt);
  auto const phi = sdf.array();
  Box const cube_nodes = amrex::surroundingNodes(cell_box);
  ParallelFor(cube_nodes, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    int const n = 4 * k + (j == 0 ? i : 3 - i);
    phi(i, j, k) = values[n];
  });

  Geometry geom(cell_box,
                RealBox({0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt}), 0,
                {0, 0, 0});
  MC::MCFab result;
  MC::marching_cubes(geom, sdf, result);

  Gpu::PinnedVector<int> host(result.m_cell_data.nComp());
  Gpu::dtoh_memcpy(host.data(), result.m_cell_data.dataPtr(),
                   host.size() * sizeof(int));
  Gpu::streamSynchronize();
  int const bit = 1 << 4; // MC33 face 5: the low-z face.
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      (host[MC::face_decision_valid_mask] & bit) != 0,
      "MC33 did not retain its resolved low-z face decision");
  AMREX_ALWAYS_ASSERT(host[MC::case_id] >= 0);
  return (host[MC::face_fluid_connected_mask] & bit) != 0;
}

void validate_mc33_face_decisions() {
  MC::Initialize();
  AMREX_ALWAYS_ASSERT(resolved_bottom_face_is_fluid_connected(
      {2.0_rt, -1.0_rt, 2.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt}));
  AMREX_ALWAYS_ASSERT(!resolved_bottom_face_is_fluid_connected(
      {1.0_rt, -2.0_rt, 1.0_rt, -2.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt}));
  AMREX_ALWAYS_ASSERT(resolved_bottom_face_is_fluid_connected(
      {1.0_rt, -1.0_rt, 1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt, -1.0_rt}));
}

void validate_narrow_band_levelset() {
  Box const domain(IntVect(0), IntVect(15));
  Geometry const geom(
      domain, RealBox({-1.0_rt, -1.0_rt, -1.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt}), 0,
      {0, 0, 0});
  BoxArray ba(domain);
  ba.maxSize(8);
  DistributionMapping dm(ba);
  BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
  MultiFab exact(nba, dm, 1, 2);
  MultiFab narrow(nba, dm, 1, 2);

  STLtools stl;
  stl.setUseMarchingCubes(true);
  stl.read_stl_file("cube.stl", 0.43_rt, {0.037_rt, -0.021_rt, 0.015_rt}, 0);
  stl.fillSignedDistance(exact, exact.nGrowVect(), geom);
  exact.OverrideSync(geom.periodicity());
  exact.FillBoundary(geom.periodicity());
  stl.fillMarchingCubesLevelSet(narrow, narrow.nGrowVect(), geom);

  Gpu::Buffer<Long> errors({0L, 0L});
  Long *const error = errors.data();
  for (MFIter mfi(exact); mfi.isValid(); ++mfi) {
    auto const full = exact.const_array(mfi);
    auto const band = narrow.const_array(mfi);
    ParallelFor(mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  if ((full(i, j, k) > 0.0_rt) != (band(i, j, k) > 0.0_rt)) {
                    Gpu::Atomic::AddNoRet(error, 1L);
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
        for (int kk = -1; kk <= 2; ++kk) {
          for (int jj = -1; jj <= 2; ++jj) {
            for (int ii = -1; ii <= 2; ++ii) {
              Real const a = full(i + ii, j + jj, k + kk);
              Real const b = band(i + ii, j + jj, k + kk);
              Real const tolerance = 16.0_rt *
                                     std::numeric_limits<Real>::epsilon() *
                                     amrex::max(1.0_rt, std::abs(a));
              if (std::abs(a - b) > tolerance) {
                Gpu::Atomic::AddNoRet(error + 1, 1L);
              }
            }
          }
        }
      }
    });
  }
  errors.copyToHost();
  GpuArray<Long, 2> global_errors{
      errors.hostData()[0], errors.hostData()[1]};
  ParallelAllReduce::Sum(global_errors.data(), int(global_errors.size()),
                         ParallelContext::CommunicatorSub());
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      global_errors[0] == 0,
      "Narrow-band MC field changed the global STL sign classification");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(global_errors[1] == 0,
                                   "Narrow-band MC field is not exact on a "
                                   "cut-cell interpolation/normal stencil");
}

void validate_ascii_stl(std::string const &filename) {
  if (filename.empty()) {
    return;
  }
  ParallelDescriptor::Barrier();
  if (ParallelDescriptor::IOProcessor()) {
    std::ifstream input(filename);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        input.good(), "Could not open generated STL " + filename);
    std::string line;
    Long facets = 0;
    Long headers = 0;
    Long trailers = 0;
    Long vertex_lines = 0;
    GpuArray<GpuArray<Real, 3>, 3> vertices;
    int vertex_count = 0;
    while (std::getline(input, line)) {
      facets += line.rfind("facet normal ", 0) == 0;
      headers += line.rfind("solid Created by AMReX", 0) == 0;
      trailers += line.rfind("endsolid Created by AMReX", 0) == 0;
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          line.find("nan") == std::string::npos &&
              line.find("inf") == std::string::npos,
          "Generated STL contains a non-finite value: " + filename);
      if (line.find("vertex") != std::string::npos) {
        std::istringstream parser(line);
        std::string tag;
        parser >> tag >> vertices[vertex_count][0] >>
            vertices[vertex_count][1] >> vertices[vertex_count][2];
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
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
              cx * cx + cy * cy + cz * cz > 0.0_rt,
              "Generated STL contains a degenerate facet: " + filename);
          vertex_count = 0;
        }
      }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        facets > 0 && headers == 1 && trailers == 1 && vertex_count == 0 &&
            vertex_lines == 3 * facets,
        "Generated STL has incomplete facets or invalid solid records: " +
            filename);
  }
  ParallelDescriptor::Barrier();
}

} // namespace

void main_main ()
{
    validate_mc33_face_decisions();

    int nx = 64;
    int ny = 64;
    int nz = 64;
    int max_grid_size = 32;
    int required_coarsening_level = 0;
    int max_coarsening_level = 0;
    bool build_coarse_level_by_coarsening = true;
    Real xmin = -1.2;
    Real xmax =  1.2;
    Real ymin = -1.2;
    Real ymax =  1.2;
    Real zmin = -1.2;
    Real zmax =  1.2;
    std::vector<int> is_periodic{0,0,0};
    bool cleanup_test = false;
    bool bunny_test = false;
    bool narrow_band_test = false;
    bool extend_domain_test = false;
    std::string eb_method("marching_cubes");
    std::string stl_file("cube.stl");
    Real stl_scale = 1.0;
    std::vector<Real> stl_center{0.0, 0.0, 0.0};
    {
        ParmParse pp;
        pp.query("nx", nx);
        pp.query("ny", ny);
        pp.query("nz", nz);
        pp.query("max_grid_size", max_grid_size);
        pp.query("xmin",xmin);
        pp.query("xmax",xmax);
        pp.query("ymin",ymin);
        pp.query("ymax",ymax);
        pp.query("zmin",zmin);
        pp.query("zmax",zmax);
        pp.queryarr("is_periodic",is_periodic);
        pp.query("required_coarsening_level",required_coarsening_level);
        pp.query("max_coarsening_level",max_coarsening_level);
        pp.query("build_coarse_level_by_coarsening",
                 build_coarse_level_by_coarsening);
        pp.query("cleanup_test", cleanup_test);
        pp.query("bunny_test", bunny_test);
        pp.query("narrow_band_test", narrow_band_test);
        pp.query("extend_domain_test", extend_domain_test);
        pp.query("test_stl_file", stl_file);
        pp.query("test_stl_scale", stl_scale);
        pp.queryarr("test_stl_center", stl_center);

        ParmParse ppeb2("eb2");
        std::string geom_type("stl");
        ppeb2.queryAdd("geom_type", geom_type);
        ppeb2.queryAdd("stl_file", stl_file);
        ppeb2.queryAdd("stl_scale", stl_scale);
        ppeb2.queryAdd("stl_center", stl_center);
        ppeb2.queryAdd("stl_geometry_method",eb_method);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            eb_method == "legacy" || eb_method == "marching_cubes",
            "eb2.stl_geometry_method must be legacy or marching_cubes");
    }

    if (narrow_band_test) {
        validate_narrow_band_levelset();
    }

    Geometry geom(Box(IntVect(0),IntVect(AMREX_D_DECL(nx-1,ny-1,nz-1))),
                  RealBox({AMREX_D_DECL(xmin,ymin,zmin)},
                          {AMREX_D_DECL(xmax,ymax,zmax)}),
                  0, {AMREX_D_DECL(is_periodic[0],is_periodic[1],is_periodic[2])});
    BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    double t0 = amrex::second();
    EB2::Build(geom,required_coarsening_level,max_coarsening_level,4,
               build_coarse_level_by_coarsening);
    double t1 = amrex::second();
    amrex::Print() << "EB method: " << eb_method
                   << ", build time: " << t1-t0 << "\n";
    auto factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        factory->maxCoarseningLevel() >= required_coarsening_level,
        "Marching-cubes EB did not build every required coarse level");
    std::string eb_surface_stl_file;
    ParmParse("eb2").query("eb_surface_stl_file", eb_surface_stl_file);
    if (!eb_surface_stl_file.empty()) {
      WriteEBSurfaceSTL(ba, dm, geom, factory.get(), eb_surface_stl_file);
      validate_ascii_stl(eb_surface_stl_file);
    }
    std::string mc_stl_file;
    ParmParse("eb2").query("mc_stl_file", mc_stl_file);
    validate_ascii_stl(mc_stl_file);
    auto const &volfrac = factory->getVolFrac();
    Real const fluid_volume =
        volfrac.sum() * geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
    Real const domain_volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    Real const expected_volume =
        domain_volume - 8.0_rt * stl_scale * stl_scale * stl_scale;
    Real const error = std::abs(fluid_volume - expected_volume);
    Real const max_dx = amrex::max(
        geom.CellSize(0), amrex::max(geom.CellSize(1), geom.CellSize(2)));
    if (bunny_test || extend_domain_test) {
      amrex::Print() << "Fluid volume: " << fluid_volume << "\n";
    } else {
      amrex::Print() << "Fluid volume: " << fluid_volume
                     << ", expected: " << expected_volume
                     << ", error: " << error << "\n";
    }
    AMREX_ALWAYS_ASSERT(volfrac.min(0) >= 0.0_rt);
    AMREX_ALWAYS_ASSERT(volfrac.max(0) <= 1.0_rt);
    if (bunny_test || extend_domain_test) {
      AMREX_ALWAYS_ASSERT(fluid_volume > 0.0_rt);
      AMREX_ALWAYS_ASSERT(fluid_volume < domain_volume);
    }

    auto volume_centroid = factory->getCentroid().ToMultiFab(0.0_rt, 0.0_rt);
    auto boundary_centroid = factory->getBndryCent().ToMultiFab(0.0_rt, 0.0_rt);
    auto boundary_normal = factory->getBndryNormal().ToMultiFab(0.0_rt, 0.0_rt);
    auto boundary_area = factory->getBndryArea().ToMultiFab(0.0_rt, 0.0_rt);
    AMREX_ALWAYS_ASSERT(!volume_centroid.contains_nan());
    AMREX_ALWAYS_ASSERT(!boundary_centroid.contains_nan());
    AMREX_ALWAYS_ASSERT(!boundary_normal.contains_nan());
    AMREX_ALWAYS_ASSERT(!boundary_area.contains_nan());
    AMREX_ALWAYS_ASSERT(!volume_centroid.contains_inf());
    AMREX_ALWAYS_ASSERT(!boundary_centroid.contains_inf());
    AMREX_ALWAYS_ASSERT(!boundary_normal.contains_inf());
    AMREX_ALWAYS_ASSERT(!boundary_area.contains_inf());
    AMREX_ALWAYS_ASSERT(boundary_area.min(0) >= 0.0_rt);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      AMREX_ALWAYS_ASSERT(volume_centroid.min(d) >= -0.5_rt);
      AMREX_ALWAYS_ASSERT(volume_centroid.max(d) <= 0.5_rt);
      // Non-cut entries in a mixed CutFab retain the EB sentinel -1.
      AMREX_ALWAYS_ASSERT(boundary_centroid.min(d) >= -1.0_rt);
      AMREX_ALWAYS_ASSERT(boundary_centroid.max(d) <= 0.5_rt);
      AMREX_ALWAYS_ASSERT(boundary_normal.min(d) >= -1.0_rt);
      AMREX_ALWAYS_ASSERT(boundary_normal.max(d) <= 1.0_rt);
    }

    auto const area_fraction = factory->getAreaFrac();
    auto const face_centroid = factory->getFaceCent();
    Array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> dense_area;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      dense_area[d] = std::make_unique<MultiFab>(
          area_fraction[d]->ToMultiFab(1.0_rt, 0.0_rt));
      auto face_cent = face_centroid[d]->ToMultiFab(0.0_rt, 0.0_rt);
      AMREX_ALWAYS_ASSERT(!dense_area[d]->contains_nan());
      AMREX_ALWAYS_ASSERT(!face_cent.contains_nan());
      AMREX_ALWAYS_ASSERT(!dense_area[d]->contains_inf());
      AMREX_ALWAYS_ASSERT(!face_cent.contains_inf());
      AMREX_ALWAYS_ASSERT(dense_area[d]->min(0) >= 0.0_rt);
      AMREX_ALWAYS_ASSERT(dense_area[d]->max(0) <= 1.0_rt);
      for (int n = 0; n < AMREX_SPACEDIM - 1; ++n) {
        AMREX_ALWAYS_ASSERT(face_cent.min(n) >= -0.5_rt);
        AMREX_ALWAYS_ASSERT(face_cent.max(n) <= 0.5_rt);
      }
    }

    auto const &flags = factory->getMultiEBCellFlagFab();
    auto const &levelset = factory->getLevelSet();
    AMREX_ALWAYS_ASSERT(!levelset.contains_nan());
    AMREX_ALWAYS_ASSERT(!levelset.contains_inf());

    if (extend_domain_test) {
      IntVect const probe_iv(
          AMREX_D_DECL(geom.Domain().smallEnd(0) - 1,
                       (geom.Domain().smallEnd(1) + geom.Domain().bigEnd(1) + 1) / 2,
                       (geom.Domain().smallEnd(2) + geom.Domain().bigEnd(2) + 1) / 2));
      Gpu::Buffer<Long> extension_counts({0L, 0L});
      Long *const extension = extension_counts.data();
      for (MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        Box const probe = Box(probe_iv, probe_iv) & levelset[mfi].box();
        auto const phi = levelset.const_array(mfi);
        ParallelFor(probe, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          Gpu::Atomic::AddNoRet(extension, 1L);
          if (phi(i, j, k) >= 0.0_rt) {
            Gpu::Atomic::AddNoRet(extension + 1, 1L);
          }
        });
      }
      extension_counts.copyToHost();
      Long extension_total = extension_counts.hostData()[0];
      Long extension_covered = extension_counts.hostData()[1];
      ParallelAllReduce::Sum<Long>({extension_total, extension_covered},
                                   ParallelContext::CommunicatorSub());
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          extension_total > 0 && extension_covered == extension_total,
          "MC domain-face extension did not keep the exterior probe covered");
    }

    auto const edge_centroid = factory->getEdgeCent();
    Array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> dense_edge;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      dense_edge[d] = std::make_unique<MultiFab>(
          edge_centroid[d]->ToMultiFab(1.0_rt, -1.0_rt));
      AMREX_ALWAYS_ASSERT(!dense_edge[d]->contains_nan());
      AMREX_ALWAYS_ASSERT(!dense_edge[d]->contains_inf());
      AMREX_ALWAYS_ASSERT(dense_edge[d]->min(0) >= -1.0_rt);
      AMREX_ALWAYS_ASSERT(dense_edge[d]->max(0) <= 1.0_rt);
    }

    Gpu::Buffer<int> levelset_errors({0, 0});
    int *const level_errors = levelset_errors.data();
    for (MFIter mfi(flags); mfi.isValid(); ++mfi) {
      auto const flag = flags.const_array(mfi);
      auto const phi = levelset.const_array(mfi);
      ParallelFor(mfi.validbox(),
                  [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
        ParallelFor(
            mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
    GpuArray<int, 2> global_levelset_errors{
        levelset_errors.hostData()[0], levelset_errors.hostData()[1]};
    ParallelAllReduce::Sum(global_levelset_errors.data(),
                           int(global_levelset_errors.size()),
                           ParallelContext::CommunicatorSub());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_levelset_errors[0] == 0,
        "A covered cell retains a public negative-in-fluid level-set node");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_levelset_errors[1] == 0,
        "Edge centroids disagree with repaired public level-set signs");

    Gpu::Buffer<Long> zero_node_count({0L});
    Long *const zero_nodes = zero_node_count.data();
    for (MFIter mfi(levelset); mfi.isValid(); ++mfi) {
      auto const phi = levelset.const_array(mfi);
      ParallelFor(mfi.validbox(),
                  [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (phi(i, j, k) == 0.0_rt) {
                      Gpu::Atomic::AddNoRet(zero_nodes, 1L);
                    }
                  });
    }
    zero_node_count.copyToHost();
    Long repaired_nodes = zero_node_count.hostData()[0];
    ParallelAllReduce::Sum(repaired_nodes, ParallelContext::CommunicatorSub());

    Real small_volfrac = 0.0_rt;
    ParmParse("eb2").query("small_volfrac", small_volfrac);
    Gpu::Buffer<int> topology_counts({0, 0, 0, 0, 0, 0});
    int *const topology = topology_counts.data();
    Box const domain = geom.Domain();
    for (MFIter mfi(flags); mfi.isValid(); ++mfi) {
      Box const bx = mfi.validbox();
      auto const flag = flags.const_array(mfi);
      auto const apx = dense_area[0]->const_array(mfi);
      auto const apy = dense_area[1]->const_array(mfi);
      auto const apz = dense_area[2]->const_array(mfi);
      auto const vf = volfrac.const_array(mfi);
      auto const ba = boundary_area.const_array(mfi);
      auto const bn = boundary_normal.const_array(mfi);
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (flag(i, j, k).isSingleValued()) {
          Gpu::Atomic::AddNoRet(topology + 1, 1);
          if (vf(i, j, k) < small_volfrac) {
            Gpu::Atomic::AddNoRet(topology + 2, 1);
          }
          if (vf(i, j, k) == 1.0_rt && ba(i, j, k) > 0.0_rt) {
            Gpu::Atomic::AddNoRet(topology + 4, 1);
          }
          Real const nx = apx(i, j, k) - apx(i + 1, j, k);
          Real const ny = apy(i, j, k) - apy(i, j + 1, k);
          Real const nz = apz(i, j, k) - apz(i, j, k + 1);
          Real const norm = std::sqrt(nx * nx + ny * ny + nz * nz);
          if (norm <= 1.e-12_rt ||
              std::abs(norm * bn(i, j, k, 0) - nx) > 5.e-12_rt ||
              std::abs(norm * bn(i, j, k, 1) - ny) > 5.e-12_rt ||
              std::abs(norm * bn(i, j, k, 2) - nz) > 5.e-12_rt) {
            Gpu::Atomic::AddNoRet(topology + 5, 1);
          }
        }
        if (flag(i, j, k).isMultiValued()) {
          Gpu::Atomic::AddNoRet(topology + 3, 1);
        }

        int mismatches = 0;
        if (i > domain.smallEnd(0)) {
          mismatches +=
              flag(i, j, k).isConnected(-1, 0, 0) != (apx(i, j, k) > 0.0_rt);
        }
        if (i < domain.bigEnd(0)) {
          mismatches +=
              flag(i, j, k).isConnected(1, 0, 0) != (apx(i + 1, j, k) > 0.0_rt);
        }
        if (j > domain.smallEnd(1)) {
          mismatches +=
              flag(i, j, k).isConnected(0, -1, 0) != (apy(i, j, k) > 0.0_rt);
        }
        if (j < domain.bigEnd(1)) {
          mismatches +=
              flag(i, j, k).isConnected(0, 1, 0) != (apy(i, j + 1, k) > 0.0_rt);
        }
        if (k > domain.smallEnd(2)) {
          mismatches +=
              flag(i, j, k).isConnected(0, 0, -1) != (apz(i, j, k) > 0.0_rt);
        }
        if (k < domain.bigEnd(2)) {
          mismatches +=
              flag(i, j, k).isConnected(0, 0, 1) != (apz(i, j, k + 1) > 0.0_rt);
        }
        if (mismatches != 0) {
          Gpu::Atomic::AddNoRet(topology, mismatches);
        }
      });
    }
    topology_counts.copyToHost();
    GpuArray<int, 6> global_topology_counts;
    for (int n = 0; n < int(global_topology_counts.size()); ++n) {
      global_topology_counts[n] = topology_counts.hostData()[n];
    }
    ParallelAllReduce::Sum(global_topology_counts.data(),
                           int(global_topology_counts.size()),
                           ParallelContext::CommunicatorSub());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[0] == 0,
        "EBCellFlag coordinate connectivity disagrees with MC face apertures");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[1] > 0,
        "Marching-cubes test did not produce any single-valued cells");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[2] == 0,
        "Nodal repair left a single-valued cell below eb2.small_volfrac");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[3] == 0,
        "Single-valued EB regression retained a multi-valued cell");
    if (cleanup_test) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            repaired_nodes > 0,
            "Legacy-style MC repair did not move any fluid node to zero");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        global_topology_counts[5] == 0,
        "EB boundary normals disagree with repaired face-aperture closure");
    amrex::Print() << "Single-valued cells: "
                   << global_topology_counts[1]
                   << ", retained multi-valued cells: "
                   << global_topology_counts[3] << "\n";
    if (!cleanup_test && !bunny_test && !extend_domain_test) {
        Real const domain_scale = amrex::max(
            xmax-xmin,amrex::max(ymax-ymin,zmax-zmin));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            error < 6.0_rt*max_dx*max_dx*domain_scale,
            "Marching-cubes cut-cell volume error exceeds the expected O(dx^2) bound");
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    main_main();
    amrex::Finalize();
}
