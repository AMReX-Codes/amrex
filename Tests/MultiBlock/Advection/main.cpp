
#include "AMReX_NonLocalBC.H"

#include "AMReX.H"
#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"

#include "AMReX_PlotFileUtil.H"

void MyMain();

int main(int argc, char** argv) {
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif
    // Let me throw exceptions for triggering my debugger
    amrex::Initialize(MPI_COMM_WORLD, std::cout, std::cerr, [](const char* msg) { throw std::runtime_error(msg); });
    MyMain();
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}

using namespace amrex;

enum idirs { ix, iy };
enum num_components { three_components = 3 };

static constexpr IntVect e_x = IntVect::TheDimensionVector(ix);
static constexpr IntVect e_y = IntVect::TheDimensionVector(iy);
class AdvectionAmrCore : public AmrCore {

  public:
    AdvectionAmrCore(Direction vel, Geometry const& level_0_geom,
                     AmrInfo const& amr_info = AmrInfo())
        : AmrCore(level_0_geom, amr_info), velocity{vel} {
        AmrCore::InitFromScratch(0.0);
        InitData();  // CUDA does not allow extended lambdas in ctors.
    }

    void InitData() {
        const auto problo = Geom(0).ProbLoArray();
        const auto dx = Geom(0).CellSizeArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mass,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real> m = mass.array(mfi);
            Array4<Real> vx = mass.array(mfi, 1);
            Array4<Real> vy = mass.array(mfi, 2);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x[] = {AMREX_D_DECL(problo[0] + (0.5+i)*dx[0],
                                         problo[1] + (0.5+j)*dx[1],
                                         problo[2] + (0.5+k)*dx[2])};
                const double r2 = AMREX_D_TERM(x[0] * x[0], + x[1] * x[1], + x[2] * x[2]);
                constexpr double R = 0.1 * 0.1;
                m(i, j, k) = r2 < R ? 1.0 : 0.0;
                vx(i, j, k) = r2 < R ? 1.0 : 0.0;
                vy(i, j, k) = 0.0;
            });
        }
    }

    void AdvanceInTime(double dt) {
        // Do first order accurate godunov splitting
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            DoOperatorSplitStep(dt, static_cast<Direction>(d));
        }
    }

    void DoOperatorSplitStep(double dt, Direction dir) {
        // Perform first order accurate upwinding with velocity 1 in the stored direction.
        const double dx = Geom(0).CellSize(0);
        const double a_dt_over_dx = dt / dx * (velocity == dir);
        if (dir == Direction::x) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
                Array4<Real> m = mass.array(mfi);
                Array4<Real> next = mass_next.array(mfi);
                ParallelFor(mfi.growntilebox(e_y), int(three_components),
                            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                                next(i, j, k, n) = m(i, j, k, n) - a_dt_over_dx * (m(i, j, k, n) - m(i - 1, j, k, n));
                            });
            }
        } else if (dir == Direction::y) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
                Array4<Real> m = mass.array(mfi);
                Array4<Real> next = mass_next.array(mfi);
                ParallelFor(mfi.growntilebox(ix), int(three_components),
                            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                                next(i, j, k, n) = m(i, j, k, n) - a_dt_over_dx * (m(i, j, k, n) - m(i, j - 1, k, n));
                            });
            }
        }
        std::swap(mass, mass_next);
    }

    MultiFab mass{};
    MultiFab mass_next{};
    Direction velocity{};

  private:
    void ErrorEst(int /*level*/, ::amrex::TagBoxArray& /*tags*/, double /*time_point*/,
                  int /* ngrow */) override {
        throw std::runtime_error("For simplicity, this example supports only one level.");
    }

    void
    MakeNewLevelFromScratch(int level, double, const ::amrex::BoxArray& box_array,
                            const ::amrex::DistributionMapping& distribution_mapping) override {
        if (level > 0) {
            throw std::runtime_error("For simplicity, this example supports only one level.");
        }
        const IntVect ngrow{AMREX_D_DECL(1, 1, 0)};
        mass.define(box_array, distribution_mapping, three_components, ngrow);
        mass_next.define(box_array, distribution_mapping, three_components, ngrow);
    }

    void MakeNewLevelFromCoarse(int /*level*/, double /*time_point*/, const ::amrex::BoxArray&,
                                const ::amrex::DistributionMapping&) override {
        throw std::runtime_error("For simplicity, this example supports only one level.");
    }

    void RemakeLevel(int /*level*/, double /*time_point*/, const ::amrex::BoxArray&,
                     const ::amrex::DistributionMapping&) override {
        throw std::runtime_error("For simplicity, this example supports only one level.");
    }

    void ClearLevel(int level) override {
        if (level > 0) {
            throw std::runtime_error("For simplicity, this example supports only one level.");
        }
        mass.clear();
    }
};

using namespace NonLocalBC;
struct OnesidedMultiBlockBoundaryFn {
  AdvectionAmrCore* dest;
  const AdvectionAmrCore* src;
  MultiBlockIndexMapping dtos;
  Box boundary_to_fill;
  std::unique_ptr<MultiBlockCommMetaData> cmd{};
  FabArrayBase::BDKey cached_dest_bd_key{};
  FabArrayBase::BDKey cached_src_bd_key{};
  ApplyDtosAndProjectionOnReciever<MultiBlockIndexMapping,
                                   MapComponents<Identity, SwapComponents<1, 2>>>
      packing{PackComponents{0, 0, three_components}, dtos};

  AMREX_NODISCARD CommHandler FillBoundary_nowait() {
    if (!cmd || cached_dest_bd_key != dest->mass.getBDKey() || cached_src_bd_key != src->mass.getBDKey()) {
        cmd = std::make_unique<MultiBlockCommMetaData>(dest->mass, boundary_to_fill, src->mass, dest->mass.nGrowVect(), dtos);
        cached_dest_bd_key = dest->mass.getBDKey();
        cached_src_bd_key = src->mass.getBDKey();
    }

    return ParallelCopy_nowait(no_local_copy, dest->mass, src->mass, *cmd, packing);
  }

  void FillBoundary_do_local_copy() const {
    AMREX_ASSERT(cmd && cached_dest_bd_key == dest->mass.getBDKey() && cached_src_bd_key == src->mass.getBDKey());
    if (cmd->m_LocTags && cmd->m_LocTags->size() > 0) {
        LocalCopy(packing, dest->mass, src->mass, *cmd->m_LocTags);
    }
  }

  void FillBoundary_finish(CommHandler handler) const {
    ParallelCopy_finish(dest->mass, std::move(handler), *cmd, packing);
  }
};

struct FillBoundaryFn {
    std::vector<OnesidedMultiBlockBoundaryFn> boundaries;

    void operator()(AdvectionAmrCore& core_x, AdvectionAmrCore& core_y) {
        enum { coarsest_level = 0 };
        core_x.mass.FillBoundary(core_x.Geom(coarsest_level).periodicity());
        core_y.mass.FillBoundary(core_y.Geom(coarsest_level).periodicity());
        std::vector<CommHandler> comms;
        for (auto& boundary : boundaries) {
            comms.push_back(boundary.FillBoundary_nowait());
        }
        for (auto& boundary : boundaries) {
            boundary.FillBoundary_do_local_copy();
        }
        const std::size_t n_boundaries = boundaries.size();
        for (std::size_t i = 0; i < n_boundaries; ++i) {
            boundaries[i].FillBoundary_finish(std::move(comms[i]));
        }
        AMREX_ASSERT(!core_x.mass.contains_nan());
    }
};

void WritePlotfiles(const AdvectionAmrCore& core_x, const AdvectionAmrCore& core_y, double time_point, int step)
{
    static const Vector<std::string> varnames{"Mass", "Vector_X", "Vector_Y"};
    int nlevels = 1;
    std::array<char, 256> x_pbuffer{};
    std::array<char, 256> y_pbuffer{};
    snprintf(x_pbuffer.data(), x_pbuffer.size(), "MultiBlock/core_x/plt%04d", step);
    snprintf(y_pbuffer.data(), y_pbuffer.size(), "MultiBlock/core_y/plt%04d", step);

    {
        Vector<const MultiFab*> mf{&core_x.mass};
        Vector<Geometry> geoms{core_x.Geom(0)};
        Vector<int> level_steps{step};
        Vector<IntVect> ref_ratio{};
        std::string plotfilename{x_pbuffer.data()};
        WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames, geoms, time_point, level_steps, ref_ratio);
    }
    {
        Vector<const MultiFab*> mf{&core_y.mass};
        Vector<Geometry> geoms{core_y.Geom(0)};
        Vector<int> level_steps{step};
        Vector<IntVect> ref_ratio{};
        std::string plotfilename{y_pbuffer.data()};
        WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames, geoms, time_point, level_steps, ref_ratio);
    }
}

void MyMain() {
    Box domain(IntVect{}, IntVect{AMREX_D_DECL(63, 63, 0)});
    RealBox real_box1{{AMREX_D_DECL(-0.5, -0.3, 0.0)}, {AMREX_D_DECL(0.5, 0.7, 1.0)}};
    RealBox real_box2{{AMREX_D_DECL(+0.55, -0.3, 0.0)}, {AMREX_D_DECL(1.55, 0.7, 1.0)}};

    Array<int, AMREX_SPACEDIM> is_periodic1{0, 1};
    Geometry geom1{domain, real_box1, CoordSys::cartesian, is_periodic1};

    Array<int, AMREX_SPACEDIM> is_periodic2{1, 0};
    Geometry geom2{domain, real_box2, CoordSys::cartesian, is_periodic2};

    AmrInfo amr_info{};
#if AMREX_SPACEDIM > 2
    amr_info.blocking_factor[0][2] = 1;
#endif

    AdvectionAmrCore core_x(Direction::x, geom1, amr_info);
    AdvectionAmrCore core_y(Direction::y, geom2, amr_info);

    std::vector<OnesidedMultiBlockBoundaryFn> multi_block_boundaries{};
    {   // Fill right boundary of core_x with lower mirror data of core_y
        NonLocalBC::MultiBlockIndexMapping dtos{};
        dtos.permutation = IntVect{AMREX_D_DECL(1, 0, 2)};
        dtos.offset = (domain.bigEnd(iy) + 1) * e_y + domain.bigEnd(ix) * e_x;
        dtos.sign = IntVect{AMREX_D_DECL(-1, 1, 1)};
        Box right_boundary_to_fill_in_x = grow(shift(Box{domain.bigEnd(ix) * e_x, domain.bigEnd()}, e_x), e_y);
        multi_block_boundaries.push_back({&core_x, &core_y, dtos, right_boundary_to_fill_in_x});
    } { // Fill lower boundary of core_y with right mirror data of core_x
        NonLocalBC::MultiBlockIndexMapping dtos{};
        dtos.permutation = IntVect{AMREX_D_DECL(1, 0, 2)};
        dtos.offset = domain.bigEnd(iy) * e_y - (domain.bigEnd(ix) + 1) * e_x;
        dtos.sign = IntVect{AMREX_D_DECL(1, -1, 1)};
        Box lower_boundary_to_fill_in_y = grow(shift(Box{domain.smallEnd(), domain.bigEnd() - domain.bigEnd(iy) * e_y}, -e_y), e_x);
        multi_block_boundaries.push_back({&core_y, &core_x, dtos, lower_boundary_to_fill_in_y});
    } { // Fill left boundary of core_x with upper mirror data of core_y
        NonLocalBC::MultiBlockIndexMapping dtos{};
        dtos.permutation = IntVect{AMREX_D_DECL(1, 0, 2)};
        dtos.offset = -(domain.bigEnd(iy) + 1) * e_y;
        Box left_boundary_to_fill_in_x = grow(shift(Box{domain.smallEnd(), domain.bigEnd() - domain.bigEnd(ix) * e_x}, -e_x), e_y);
        multi_block_boundaries.push_back({&core_x, &core_y, dtos, left_boundary_to_fill_in_x});
    } { // Fill upper boundary of core_y with left mirror data of core_x
        NonLocalBC::MultiBlockIndexMapping dtos{};
        dtos.permutation = IntVect{AMREX_D_DECL(1, 0, 2)};
        dtos.offset = (domain.bigEnd(ix) + 1) * e_x;
        Box upper_boundary_to_fill_in_y = grow(shift(Box{domain.bigEnd(iy) * e_y, domain.bigEnd()}, e_y), e_x);
        multi_block_boundaries.push_back({&core_y, &core_x, dtos, upper_boundary_to_fill_in_y});
    }

    FillBoundaryFn FillBoundary{std::move(multi_block_boundaries)};

    int step = 1;
    const double min_dx1_dy2 = std::min(geom1.CellSize(0), geom2.CellSize(1));
    const double cfl = 1.0;
    const double dt = cfl * min_dx1_dy2;
    const double final_time = 4.0;
    double time_point = 0.0;

    WritePlotfiles(core_x, core_y, time_point, step);
    while (time_point < final_time) {
        FillBoundary(core_x, core_y);

        core_x.AdvanceInTime(dt);
        core_y.AdvanceInTime(dt);

        amrex::Print() << "Step #" << step << ", Time Point = " << time_point << '\n';

        time_point += dt;
        step += 1;
        WritePlotfiles(core_x, core_y, time_point, step);
    }
}
