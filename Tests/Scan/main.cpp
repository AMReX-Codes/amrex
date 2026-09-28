#include <AMReX_Scan.H>

#include <AMReX.H>
#include <AMReX_GpuContainers.H>

#include <cstdint>
#include <limits>
#include <type_traits>

using namespace amrex;

namespace {

template <typename N, typename Type>
void test_scan (N n, Type type)
{
    Gpu::DeviceVector<int> result(n, -1);
    auto* p = result.data();
    int total = Scan::PrefixSum<int>(n,
        [=] AMREX_GPU_DEVICE (N i) -> int {
            AMREX_ALWAYS_ASSERT(Long(i) >= 0 && Long(i) < Long(n));
            return i % 3 + 1;
        },
        [=] AMREX_GPU_DEVICE (N i, int value) {
            AMREX_ALWAYS_ASSERT(Long(i) >= 0 && Long(i) < Long(n));
            p[i] = value;
        }, type);

    Vector<int> host(n);
    Gpu::copy(Gpu::deviceToHost, result.begin(), result.end(), host.begin());
    int expected = 0;
    for (int i = 0; i < int(n); ++i) {
        if constexpr (std::is_same_v<Type, Scan::Type::Inclusive>) {
            expected += i % 3 + 1;
        }
        AMREX_ALWAYS_ASSERT(host[i] == expected);
        if constexpr (std::is_same_v<Type, Scan::Type::Exclusive>) {
            expected += i % 3 + 1;
        }
    }
    AMREX_ALWAYS_ASSERT(total == expected);
}

template <typename N>
void test_index_type ()
{
    // Narrow indices reach the overflow boundary with small allocations.
    for (N n : {N(0), N(1), N(255), N(3072),
                N(std::numeric_limits<N>::max()-100), std::numeric_limits<N>::max()}) {
        test_scan(n, Scan::Type::inclusive);
        test_scan(n, Scan::Type::exclusive);
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    test_index_type<std::int16_t>();
    test_index_type<std::uint16_t>();
    amrex::Finalize();
}
