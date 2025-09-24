#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_IntSuperAccumulator.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <array>
#include <cstdint>
#include <limits>
#include <random>

namespace
{

struct Range
{
    float min_value;
    float max_value;
    int   count;
    int   seed;
};

} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;

    {
        using amrex::IntSuperAccumulatorFab;

        const amrex::Box bx(amrex::IntVect::TheZeroVector(), amrex::IntVect::TheZeroVector());
        IntSuperAccumulatorFab accumulator(bx, 1);

        const std::array<Range, 4> ranges = {{
            {1.0e-12f, 1.0e-3f,   32768, 11},
            {1.0e-4f,  1.0f,      65536, 23},
            {1.0f,     1.0e6f,    65536, 37},
            {1.0e-2f,  1.0e8f,   131072, 41}
        }};

        auto compute_reference = [] (amrex::Gpu::HostVector<float> const& values) -> float {
            volatile double sum = 0.0;
            volatile double corr = 0.0;
            for (float v : values) {
                const volatile double y = static_cast<double>(v) - corr;
                const volatile double t = sum + y;
                corr = (t - sum) - y;
                sum = t;
            }
            return static_cast<float>(sum);
        };

        auto accumulate_values = [&] (amrex::Gpu::HostVector<float> const& values) -> float {
            accumulator.reset(bx, 1);
            auto digits = accumulator.array();
            constexpr int i = 0;
            constexpr int j = 0;
            constexpr int k = 0;
            constexpr int entry = 0;

            const int count = static_cast<int>(values.size());

            if (amrex::Gpu::inLaunchRegion()) {
                amrex::Gpu::DeviceVector<float> d_values(count);
                amrex::Gpu::copy(amrex::Gpu::hostToDevice, values.begin(), values.end(), d_values.begin());

                auto ptr = d_values.data();
                auto arr = digits;

                amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) {
                    for (int n = 0; n < count; ++n) {
                        IntSuperAccumulatorFab::accumulate(arr, i, j, k, entry, ptr[n]);
                    }
                });

                amrex::Gpu::streamSynchronize();

                amrex::Gpu::DeviceScalar<float> d_result;
                auto res_ptr = d_result.dataPtr();
                amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) {
                    res_ptr[0] = IntSuperAccumulatorFab::finalize(arr, i, j, k, entry);
                });
                amrex::Gpu::streamSynchronize();
                return d_result.dataValue();
            } else {
                for (float v : values) {
                    IntSuperAccumulatorFab::accumulate(digits, i, j, k, entry, v);
                }
                return IntSuperAccumulatorFab::finalize(digits, i, j, k, entry);
            }
        };

        for (const auto& cfg : ranges) {
            std::mt19937 gen(cfg.seed);
            std::uniform_real_distribution<float> dist(cfg.min_value, cfg.max_value);

            amrex::Gpu::HostVector<float> values(cfg.count);
            for (int n = 0; n < cfg.count; ++n) {
                values[n] = dist(gen);
            }

            const float reference = compute_reference(values);
            const float result = accumulate_values(values);

            if (result != reference) {
                amrex::Print() << "Mismatch for range [" << cfg.min_value << ", "
                               << cfg.max_value << "] with count " << cfg.count
                               << ": reference=" << reference
                               << " accumulator=" << result << '\n';
                status = 1;
                break;
            }
        }

        if (status == 0) {
            constexpr float value = 0.125f;
            constexpr int repeats = 200000;

            amrex::Gpu::HostVector<float> values(repeats);
            for (int n = 0; n < repeats; ++n) {
                values[n] = value;
            }

            const float reference = compute_reference(values);
            const float result = accumulate_values(values);

            if (result != reference) {
                amrex::Print() << "Mismatch for constant accumulation: reference="
                               << reference << " accumulator=" << result << '\n';
                status = 1;
            }
        }
    }

    amrex::Finalize();
    return status;
}
