
#include "AMReX_NonLocalBC.H"

#include "AMReX.H"
#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"

#include "AMReX_PlotFileUtil.H"


void MyMain();

int main(int argc, char** argv) try {
    MPI_Init(nullptr, nullptr);
    amrex::Initialize(MPI_COMM_WORLD, std::cout, std::cerr, [](const char* msg) { throw std::runtime_error(msg); });
    MyMain();
    amrex::Finalize();
    MPI_Finalize();
} catch (const std::exception& e) {
    std::cerr << "An exception occured: " << e.what() << "\nFinalizing AMReX.\n";
    amrex::Finalize();
    MPI_Finalize();
} catch (...) {
    std::cerr << "An unusual exception occured. Finalizing AMReX.\n";
    amrex::Finalize();
    MPI_Finalize();
}

using namespace amrex;

class AdvectionAmrCore : public AmrCore {
  enum { one_component = 1 };
  
  public:
    AdvectionAmrCore(Direction vel, Geometry const& level_0_geom,
                     AmrInfo const& amr_info = AmrInfo())
        : AmrCore(level_0_geom, amr_info), velocity{vel} {
        AmrCore::InitFromScratch(0.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
            Array4<Real> m = mass.array(mfi);
            LoopConcurrentOnCpu(mfi.tilebox(), [level_0_geom, m](int i, int j, int k) {
                Real x[AMREX_SPACEDIM] = {};
                level_0_geom.CellCenter(IntVect{AMREX_D_DECL(i, j, k)}, x);
                const double r2 = AMREX_D_TERM(x[0] * x[0], + x[1] * x[1], + x[2] * x[2]);
                m(i, j, k) = r2 < 0.4 * 0.4 ? 1.0 : 0.0;
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
#pragma omp parallel
#endif
        for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
            Array4<Real> m = mass.array(mfi, 0);
            Array4<Real> next = mass_next.array(mfi, 0);
            LoopConcurrentOnCpu(mfi.tilebox(), [=](int i, int j, int k) {
                next(i, j, k) = m(i, j, k) + a_dt_over_dx * m(i - 1, j, k);
            });
        }
        } else if (dir == Direction::y) {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
                Array4<Real> m = mass.array(mfi, 0);
                Array4<Real> next = mass_next.array(mfi, 0);
                LoopConcurrentOnCpu(mfi.tilebox(), [=](int i, int j, int k) {
                    next(i, j, k) = m(i, j, k) + a_dt_over_dx * m(i, j - 1, k);
                });
            }
        }
        std::swap(mass, mass_next);
    }

    MultiFab mass{};
    MultiFab mass_next{};
    Direction velocity{};

  private:
    void ErrorEst(int level, ::amrex::TagBoxArray& tags, double time_point,
                  int /* ngrow */) override {
        throw std::runtime_error("For simplicity, this example supports only one level.");
    }

    void
    MakeNewLevelFromScratch(int level, double, const ::amrex::BoxArray& box_array,
                            const ::amrex::DistributionMapping& distribution_mapping) override {
        if (level > 0) {
            throw std::runtime_error("For simplicity, this example supports only one level.");
        }
        mass.define(box_array, distribution_mapping, one_component, ::amrex::IntVect{AMREX_D_DECL(1, 1, 0)});
        mass_next.define(box_array, distribution_mapping, one_component, ::amrex::IntVect{AMREX_D_DECL(1, 1, 0)});
    }

    void MakeNewLevelFromCoarse(int level, double time_point, const ::amrex::BoxArray& box_array,
                                const ::amrex::DistributionMapping& distribution_mapping) override {
        throw std::runtime_error("For simplicity, this example supports only one level.");
    }

    void RemakeLevel(int level, double time_point, const ::amrex::BoxArray& box_array,
                     const ::amrex::DistributionMapping& distribution_mapping) override {
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
struct OnsidedMultiBlockBoundaryFn {

  MultiBlockIndexMapping dtos;
  Box boundary_to_fill;
  std::unique_ptr<MultiBlockCommMetaData> cmd{};
  FabArrayBase::BDKey cached_dest_bd_key{};
  FabArrayBase::BDKey cached_src_bd_key{};
  ApplyDtosAndProjectionOnReciever<MultiBlockIndexMapping> packing{PackComponents{0, 0, 1}, dtos};

  AMREX_NODISCARD CommHandler FillBoundary_nowait(AdvectionAmrCore& dest, const AdvectionAmrCore& src) {
    if (!cmd || cached_dest_bd_key != dest.mass.getBDKey() || cached_src_bd_key != src.mass.getBDKey()) {
        cmd = std::make_unique<MultiBlockCommMetaData>(dest.mass, boundary_to_fill, src.mass, dest.mass.nGrowVect(), dtos);
        cached_dest_bd_key = dest.mass.getBDKey();
        cached_src_bd_key = src.mass.getBDKey();
    }
    
    return ParallelCopy_nowait(dest.mass, src.mass, *cmd, packing);
  }

  void FillBoundary_local_copy(AdvectionAmrCore& dest, const AdvectionAmrCore& src) const {
    AMREX_ASSERT(cmd && cached_dest_bd_key == dest.mass.getBDKey() && cached_src_bd_key == src.mass.getBDKey());
    if (cmd->m_LocTags && cmd->m_LocTags->size() > 0) {
        LocalCopy(packing, dest.mass, src.mass, *cmd->m_LocTags);
    }
  }

  void FillBoundary_finish(AdvectionAmrCore& dest, const AdvectionAmrCore& src, CommHandler handler) const {
    AMREX_ASSERT(cmd && cached_dest_bd_key == dest.mass.getBDKey() && cached_src_bd_key == src.mass.getBDKey());
    ParallelCopy_finish(no_local_copy, dest.mass, src.mass, std::move(handler), *cmd, packing);
  }
};

template <typename T> std::add_const_t<T> as_const(T&& x) noexcept { return x; }

struct FillBoundaryFn {
    enum { coarsest_level = 0 };

    OnsidedMultiBlockBoundaryFn x_to_y;
    OnsidedMultiBlockBoundaryFn y_to_x;

    void operator()(AdvectionAmrCore& core_x, AdvectionAmrCore& core_y) {
        core_x.mass.FillBoundary(core_x.Geom(coarsest_level).periodicity());
        core_y.mass.FillBoundary(core_y.Geom(coarsest_level).periodicity());
        // Initiate communications
        NonLocalBC::CommHandler x_to_y_comm = x_to_y.FillBoundary_nowait(core_y, as_const(core_x));
        NonLocalBC::CommHandler y_to_x_comm = y_to_x.FillBoundary_nowait(core_x, as_const(core_y));
        // Do local work and hope for some overlap with communication (throw in some MPI_Tests?)
        x_to_y.FillBoundary_local_copy(core_y, as_const(core_x));
        y_to_x.FillBoundary_local_copy(core_x, as_const(core_y));
        // Wait for all communication to be done
        x_to_y.FillBoundary_finish(core_y, as_const(core_x), std::move(x_to_y_comm));
        y_to_x.FillBoundary_finish(core_x, as_const(core_y), std::move(y_to_x_comm));
    }
};

enum idirs { ix, iy };

void MyMain() {
    Box domain(IntVect{}, IntVect{AMREX_D_DECL(63, 63, 0)});
    RealBox real_box{{AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)}};
    static constexpr auto cartesian = static_cast<int>(CoordSys::CoordType::cartesian);

    Array<int, AMREX_SPACEDIM> is_periodic1{0, 1};
    Geometry geom1{domain, real_box, cartesian, is_periodic1};

    Array<int, AMREX_SPACEDIM> is_periodic2{1, 0};
    Geometry geom2{domain, real_box, cartesian, is_periodic2};

    AmrInfo amr_info{};
    amr_info.blocking_factor[0][2] = 1;

    AdvectionAmrCore core_x(Direction::x, geom1, amr_info);
    AdvectionAmrCore core_y(Direction::y, geom2, amr_info);

    Vector<std::string> var_names{"Mass"};

    WriteSingleLevelPlotfile("MultiBlock/core_x/plt0000", core_x.mass, var_names, core_x.Geom(0), 0.0, 0);
    WriteSingleLevelPlotfile("MultiBlock/core_y/plt0000", core_y.mass, var_names, core_y.Geom(0), 0.0, 0);
    std::array<char, 256> pbuffer{};

    const IntVect e_x = IntVect::TheDimensionVector(ix);
    const IntVect e_y = IntVect::TheDimensionVector(iy);

    NonLocalBC::MultiBlockDestToSrc dtos_x{};
    dtos_x.permutation = IntVect{AMREX_D_DECL(1, 0, 2)};
    dtos_x.offset = (domain.bigEnd(iy) + 1) * e_y;
    NonLocalBC::MultiBlockIndexMapping y_to_x{dtos_x};
    Box right_boundary_to_fill_in_x = grow(shift(Box{domain.bigEnd(ix) * e_x, domain.bigEnd()}, e_x), e_y);

    NonLocalBC::MultiBlockDestToSrc dtos_y{};
    dtos_y.permutation = IntVect{AMREX_D_DECL(1, 0, 2)};
    dtos_y.offset = -(domain.bigEnd(ix) + 1) * e_x;
    NonLocalBC::MultiBlockIndexMapping x_to_y{dtos_y};
    Box lower_boundary_to_fill_in_y = grow(shift(Box{domain.smallEnd(), domain.bigEnd() - domain.bigEnd(iy) * e_y}, -e_y), e_x);

    OnsidedMultiBlockBoundaryFn x_to_y_FB{x_to_y, lower_boundary_to_fill_in_y};
    OnsidedMultiBlockBoundaryFn y_to_x_FB{y_to_x, right_boundary_to_fill_in_x};

    FillBoundaryFn FillBoundary{std::move(x_to_y_FB), std::move(y_to_x_FB)};

    int step = 1;
    const double min_dx1_dy2 = std::min(geom1.CellSize(0), geom2.CellSize(1));
    const double cfl = 0.9;
    const double dt = cfl * min_dx1_dy2;
    const double final_time = 2.0;
    double time = 0.0;

    while (time < final_time) {
        FillBoundary(core_x, core_y);

        core_x.AdvanceInTime(dt);
        core_y.AdvanceInTime(dt);

        amrex::Print() << "Step #" << step << ", Time = " << time << '\n';

        snprintf(pbuffer.data(), pbuffer.size(), "MultiBlock/core_x/plt%04d", step);
        WriteSingleLevelPlotfile(std::string{pbuffer.data()}, core_x.mass, var_names, core_x.Geom(0), 0.0, 0);
        snprintf(pbuffer.data(), pbuffer.size(), "MultiBlock/core_y/plt%04d", step);
        WriteSingleLevelPlotfile(std::string{pbuffer.data()}, core_y.mass, var_names, core_y.Geom(0), 0.0, 0);

        time += dt;
        step += 1;
    }
}