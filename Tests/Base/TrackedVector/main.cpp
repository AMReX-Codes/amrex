#include <AMReX.H>
#include <AMReX_TrackedVector.H>

#include <iostream>
#include <memory>
#include <vector>

using namespace amrex;

using TVec = Gpu::TrackedVector<int>;
using Status = TVec::Status;

namespace {

void verify_host_device_match (TVec const & v)
{
    // Calling host_const() and device_const() auto-syncs both sides
    const auto n = v.host_const().size();
    AMREX_ALWAYS_ASSERT(v.device_const().size() == n);
    if (n == 0) { return; }
    std::vector<int> tmp(n);
    Gpu::copy(Gpu::deviceToHost, v.device_const().begin(), v.device_const().end(), tmp.begin());
    for (std::size_t i = 0; i < n; ++i) {
        AMREX_ALWAYS_ASSERT(tmp[i] == v.host_const()[i]);
    }
}

void fill_device_linear (TVec& v, int base)
{
    const int n = static_cast<int>(v.device().size());
    if (n == 0) { return; }
    int* dp = v.device().data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
        dp[i] = base + i;
    });
    Gpu::streamSynchronize();
}

void test_dirty_semantics ()
{
    TVec v;
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);

    // potential write (host): must mark dirty
    v.host().resize(3);
    v.host().at(2) = 42;
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(v.status() == Status::host_dirty);
#else
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);
#endif
    AMREX_ALWAYS_ASSERT(v.host_const().size() == 3);

    // read-only device access syncs host->device, status becomes up_to_date
    [[maybe_unused]] auto dsize = v.device_const().size();
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);

    // potential write (device): must mark dirty
    [[maybe_unused]] auto* dp = v.device().data();
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(v.status() == Status::device_dirty);
#else
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);
#endif

    // read-only host access syncs device->host, status becomes up_to_date
    auto first = v.host_const().at(2);
    AMREX_ALWAYS_ASSERT(first == 42);
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);

    // read-only (device): must not mark dirty
    [[maybe_unused]] auto const* dcp = v.device_const().data();
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);
}

void test_release_gpu ()
{
    TVec v;
    v.host().assign({1, 2, 3});
    verify_host_device_match(v);

    v.release_gpu();
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(v.status() == Status::host_dirty);
#else
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);
#endif
    AMREX_ALWAYS_ASSERT(v.host_const().size() == 3U);

    // device_const() auto-syncs from host
    verify_host_device_match(v);
}

void test_d2h ()
{
    TVec v;
    v.host().assign({0, 0, 0});

    // device() auto-syncs host->device, then fill_device_linear writes on device
    fill_device_linear(v, 100);

    // host_const() auto-syncs device->host
    AMREX_ALWAYS_ASSERT(v.host_const()[0] == 100 && v.host_const()[1] == 101 && v.host_const()[2] == 102);
    verify_host_device_match(v);
}

void test_empty ()
{
    TVec v;
    AMREX_ALWAYS_ASSERT(v.device_const().empty());
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);

    // device() on empty gives empty device vec, marks device_dirty
    // then host_const() syncs back (empty)
    fill_device_linear(v, 100);

    // Overwrite host (discards any device data)
    v.host().clear();
    AMREX_ALWAYS_ASSERT(v.device_const().empty());
    AMREX_ALWAYS_ASSERT(v.status() == Status::up_to_date);
}

void test_aggregate ()
{
    struct S {
        double x;
        TVec a, b;
    };

    std::vector<int> i = {1, 2, 3};
    std::vector<int> j = {4, 5, 6};

    [[maybe_unused]] auto ptr1 = std::shared_ptr<S>(new S{42.0, i, j});  // NOLINT(modernize-make-shared)
}

void test_copy_constructor ()
{
    TVec a;
    a.host().assign({1, 2, 3});
    // Sync to device via read-only access
    verify_host_device_match(a);

    TVec b(a);  // copy construct

    // Both have same data
    AMREX_ALWAYS_ASSERT(b.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(b.host_const()[0] == 1 && b.host_const()[1] == 2 && b.host_const()[2] == 3);
    AMREX_ALWAYS_ASSERT(b.status() == a.status());
    verify_host_device_match(b);

    // Modifying b doesn't affect a
    b.host()[0] = 99;
    AMREX_ALWAYS_ASSERT(a.host_const()[0] == 1);
}

void test_move_constructor ()
{
    TVec a;
    a.host().assign({4, 5, 6});
    // Sync to device via read-only access
    [[maybe_unused]] auto ds = a.device_const().size();

    TVec b(std::move(a));  // move construct

    // b has the data
    AMREX_ALWAYS_ASSERT(b.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(b.host_const()[0] == 4 && b.host_const()[1] == 5 && b.host_const()[2] == 6);
    AMREX_ALWAYS_ASSERT(b.status() == Status::up_to_date);
    verify_host_device_match(b);

    // a is valid but empty (got fresh shared_ptrs via swap)
    // NOLINTBEGIN(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
    AMREX_ALWAYS_ASSERT(a.host_const().empty());
    AMREX_ALWAYS_ASSERT(a.status() == Status::up_to_date);
    AMREX_ALWAYS_ASSERT(a.device_const().empty());

    // Modifying b must not affect a's status or device
    auto a_status_before = a.status();
    auto a_device_size_before = a.device_const().size();
    b.host()[0] = 99;
    // device_const() on b auto-syncs
    [[maybe_unused]] auto bds = b.device_const().size();
    AMREX_ALWAYS_ASSERT(a.status() == a_status_before);
    AMREX_ALWAYS_ASSERT(a.device_const().size() == a_device_size_before);

    // a can be reused
    a.host().assign({10, 20});
    verify_host_device_match(a);
    AMREX_ALWAYS_ASSERT(a.host_const()[0] == 10);
    // NOLINTEND(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
}

void test_copy_assignment ()
{
    TVec a;
    a.host().assign({7, 8, 9});
    // Sync via read-only access
    [[maybe_unused]] auto ds = a.device_const().size();

    TVec b;
    b.host().assign({0});
    [[maybe_unused]] auto ds2 = b.device_const().size();

    b = a;  // copy assign

    AMREX_ALWAYS_ASSERT(b.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(b.host_const()[0] == 7 && b.host_const()[1] == 8 && b.host_const()[2] == 9);
    AMREX_ALWAYS_ASSERT(b.status() == a.status());
    verify_host_device_match(b);

    // Modifying b doesn't affect a
    b.host()[0] = 99;
    AMREX_ALWAYS_ASSERT(a.host_const()[0] == 7);
}

void test_move_assignment ()
{
    TVec a;
    a.host().assign({11, 12, 13});
    // Sync via read-only access
    [[maybe_unused]] auto ds = a.device_const().size();

    TVec b;
    b.host().assign({0});
    [[maybe_unused]] auto ds2 = b.device_const().size();

    b = std::move(a);  // move assign

    // b has a's original data
    AMREX_ALWAYS_ASSERT(b.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(b.host_const()[0] == 11 && b.host_const()[1] == 12 && b.host_const()[2] == 13);
    AMREX_ALWAYS_ASSERT(b.status() == Status::up_to_date);
    verify_host_device_match(b);

    // a is valid (has b's old data via swap)
    // NOLINTBEGIN(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
    AMREX_ALWAYS_ASSERT(a.host_const().size() == 1U);
    AMREX_ALWAYS_ASSERT(a.host_const()[0] == 0);
    AMREX_ALWAYS_ASSERT(a.status() == Status::up_to_date);
    verify_host_device_match(a);

    // Modifying b must not affect a's status or device
    auto a_status_before = a.status();
    auto a_device_size_before = a.device_const().size();
    b.host()[0] = 99;
    [[maybe_unused]] auto bds = b.device_const().size();
    AMREX_ALWAYS_ASSERT(a.status() == a_status_before);
    AMREX_ALWAYS_ASSERT(a.device_const().size() == a_device_size_before);
    AMREX_ALWAYS_ASSERT(a.host_const()[0] == 0);  // a's data unchanged

    // a can be reused
    a.host().assign({30, 40});
    verify_host_device_match(a);
    // NOLINTEND(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
}

void run_tests_before_finalize ()
{
    test_dirty_semantics();
    test_release_gpu();
    test_d2h();
    test_empty();
    test_aggregate();
    test_copy_constructor();
    test_move_constructor();
    test_copy_assignment();
    test_move_assignment();
}

} // namespace

int main (int argc, char* argv[])
{
    TVec cross_session;
    cross_session.host() = {7, 8, 9};

#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    amrex::Initialize(argc, argv);
    {
        run_tests_before_finalize();

        // device_const() auto-syncs host->device
        verify_host_device_match(cross_session);
        AMREX_ALWAYS_ASSERT(cross_session.host_const()[0] == 7);
    }
    amrex::Finalize();  // calls implicitly: cross_session.release_gpu();

#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(cross_session.status() == Status::host_dirty);
#endif

    amrex::Initialize(argc, argv);
    {
        AMREX_ALWAYS_ASSERT(cross_session.host_const().size() == 3U);
        AMREX_ALWAYS_ASSERT(cross_session.host_const()[0] == 7 && cross_session.host_const()[1] == 8 &&
                            cross_session.host_const()[2] == 9);
        // GPU builds:
        // Device was released on Finalize, device access will re-sync.
        // auto-syncs host->device in new AMReX session
        verify_host_device_match(cross_session);
    }
    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    std::cout << "TrackedVector tests passed.\n";
    return 0;
}
