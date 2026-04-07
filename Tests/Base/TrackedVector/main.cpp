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

void verify_host_values (TVec const & v, std::vector<int> const & expected)
{
    auto const & host = v.host_const();
    AMREX_ALWAYS_ASSERT(host.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        AMREX_ALWAYS_ASSERT(host[i] == expected[i]);
    }
}

void test_copy_without_amrex_session (TVec const & src, std::vector<int> const & expected)
{
    TVec copy_ctor(src); // NOLINT(performance-unnecessary-copy-initialization)

    TVec copy_assign;
    copy_assign.host().assign({-1});
    copy_assign = src;

    verify_host_values(copy_ctor, expected);
    verify_host_values(copy_assign, expected);
    AMREX_ALWAYS_ASSERT(copy_ctor.status() == src.status());
    AMREX_ALWAYS_ASSERT(copy_assign.status() == src.status());
}

void fill_device_linear (TVec& v, int base)
{
    const int n = static_cast<int>(v.device().size());
    if (n == 0) { return; }
    int* dp = v.device().data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept {
        dp[i] = base + i;
    });
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

void test_release_gpu_device_dirty ()
{
    // Device is modified, then release_gpu() must sync D->H before freeing.
    TVec v;
    v.host().assign({0, 0, 0});

    // Write on device: {100, 101, 102}
    fill_device_linear(v, 100);
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(v.status() == Status::device_dirty);
#endif

    // release_gpu() should sync device->host, then free device memory
    v.release_gpu();
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(v.status() == Status::host_dirty);
#endif

    // Host must have the device-written values
    AMREX_ALWAYS_ASSERT(v.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(v.host_const()[0] == 100);
    AMREX_ALWAYS_ASSERT(v.host_const()[1] == 101);
    AMREX_ALWAYS_ASSERT(v.host_const()[2] == 102);

    // Can round-trip back to device
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
    test_release_gpu_device_dirty();
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
    TVec pre_init_copy_source;
    pre_init_copy_source.host() = {21, 22, 23};
    test_copy_without_amrex_session(pre_init_copy_source, {21, 22, 23});

    TVec cross_session;
    cross_session.host() = {7, 8, 9};

    // This vector will be device-dirty when Finalize runs.
    // ExecOnFinalize must sync D->H so the data survives.
    TVec dirty_on_finalize;
    dirty_on_finalize.host().assign({0, 0, 0});

    // Cross-session copy/move tests (review point 2):
    // Copies must register their own finalize handler so the device mirror
    // is drained and released even when the copy's device() is never called.
    TVec copy_assign_cross;
    TVec move_assign_cross;
    TVec copy_assign_dirty_cross;
    std::unique_ptr<TVec> copy_ctor_cross;
    std::unique_ptr<TVec> move_ctor_cross;

#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    amrex::Initialize(argc, argv);
    {
        run_tests_before_finalize();

        // device_const() auto-syncs host->device
        verify_host_device_match(cross_session);
        AMREX_ALWAYS_ASSERT(cross_session.host_const()[0] == 7);

        // Write on device, leave device_dirty for Finalize to handle
        fill_device_linear(dirty_on_finalize, 200);
#ifdef AMREX_USE_GPU
        AMREX_ALWAYS_ASSERT(dirty_on_finalize.status() == Status::device_dirty);
#endif
    }
    amrex::Finalize();  // calls implicitly: release_gpu() with D->H sync

#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(cross_session.status() == Status::host_dirty);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.status() == Status::host_dirty);
#endif

    // ExecOnFinalize must have synced device data to host
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[0] == 200);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[1] == 201);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[2] == 202);

    test_copy_without_amrex_session(cross_session, {7, 8, 9});
    test_copy_without_amrex_session(dirty_on_finalize, {200, 201, 202});

    // --- Session 2 ---
    amrex::Initialize(argc, argv);
    {
        AMREX_ALWAYS_ASSERT(cross_session.host_const().size() == 3U);
        AMREX_ALWAYS_ASSERT(cross_session.host_const()[0] == 7 && cross_session.host_const()[1] == 8 &&
                            cross_session.host_const()[2] == 9);
        // GPU builds:
        // Device was released on Finalize, device access will re-sync.
        // auto-syncs host->device in new AMReX session
        verify_host_device_match(cross_session);

        // dirty_on_finalize also round-trips into the new session
        verify_host_device_match(dirty_on_finalize);
        AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[0] == 200);

        // Write new device values for the session-3 finalize test.
        // This re-arms the finalize hook via to_device().
        fill_device_linear(dirty_on_finalize, 300);
#ifdef AMREX_USE_GPU
        AMREX_ALWAYS_ASSERT(dirty_on_finalize.status() == Status::device_dirty);
#endif
    }
    amrex::Finalize();

#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.status() == Status::host_dirty);
#endif

    // Re-armed finalize hook must have synced {300,301,302} back to host
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[0] == 300);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[1] == 301);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[2] == 302);

    // --- Session 3 ---
    amrex::Initialize(argc, argv);
    {
        // Verify data round-trips into a 3rd session
        verify_host_device_match(cross_session);
        verify_host_device_match(dirty_on_finalize);
        AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[0] == 300);

        // Write yet again — exercises re-arming a second time
        fill_device_linear(dirty_on_finalize, 400);
    }
    amrex::Finalize();

    // Finalize hook fired a 3rd time
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[0] == 400);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[1] == 401);
    AMREX_ALWAYS_ASSERT(dirty_on_finalize.host_const()[2] == 402);

    // --- Session 4: read-only device access after 3 finalizations ---
    amrex::Initialize(argc, argv);
    {
        verify_host_device_match(cross_session);
        verify_host_device_match(dirty_on_finalize);
    }
    amrex::Finalize();

    // --- Session 5: cross-session copy and move ---
    // Exercises all four code paths: copy ctor, copy assignment,
    // move ctor, move assignment — each must arm finalize so the
    // device mirror is drained even when device() is never called
    // on the target before Finalize.
    amrex::Initialize(argc, argv);
    {
        // Source with up_to_date status (both sides synced)
        TVec src_synced;
        src_synced.host() = {60, 61, 62};
        verify_host_device_match(src_synced);

        // Copy assignment
        copy_assign_cross = src_synced;

        // Copy constructor
        copy_ctor_cross = std::make_unique<TVec>(src_synced);

        // Move assignment
        TVec move_assign_src;
        move_assign_src.host() = {70, 71, 72};
        verify_host_device_match(move_assign_src);
        move_assign_cross = std::move(move_assign_src);

        // Move constructor
        TVec move_ctor_src;
        move_ctor_src.host() = {80, 81, 82};
        verify_host_device_match(move_ctor_src);
        move_ctor_cross = std::make_unique<TVec>(std::move(move_ctor_src));

        // Copy assignment from a device-dirty source
        TVec src_dirty;
        src_dirty.host().assign({0, 0, 0});
        fill_device_linear(src_dirty, 90);
#ifdef AMREX_USE_GPU
        AMREX_ALWAYS_ASSERT(src_dirty.status() == Status::device_dirty);
#endif
        copy_assign_dirty_cross = src_dirty;
#ifdef AMREX_USE_GPU
        AMREX_ALWAYS_ASSERT(copy_assign_dirty_cross.status() == Status::device_dirty);
#endif

        // Intentionally do NOT call device() on any target before Finalize
    }
    amrex::Finalize();

    // All five must have host_dirty status after finalize (GPU builds)
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT(copy_assign_cross.status() == Status::host_dirty);
    AMREX_ALWAYS_ASSERT(copy_ctor_cross->status() == Status::host_dirty);
    AMREX_ALWAYS_ASSERT(move_assign_cross.status() == Status::host_dirty);
    AMREX_ALWAYS_ASSERT(move_ctor_cross->status() == Status::host_dirty);
    AMREX_ALWAYS_ASSERT(copy_assign_dirty_cross.status() == Status::host_dirty);
#endif

    // Copy assignment: host data intact
    AMREX_ALWAYS_ASSERT(copy_assign_cross.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(copy_assign_cross.host_const()[0] == 60);
    AMREX_ALWAYS_ASSERT(copy_assign_cross.host_const()[1] == 61);
    AMREX_ALWAYS_ASSERT(copy_assign_cross.host_const()[2] == 62);

    // Copy constructor: host data intact
    AMREX_ALWAYS_ASSERT(copy_ctor_cross->host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(copy_ctor_cross->host_const()[0] == 60);
    AMREX_ALWAYS_ASSERT(copy_ctor_cross->host_const()[1] == 61);
    AMREX_ALWAYS_ASSERT(copy_ctor_cross->host_const()[2] == 62);

    // Move assignment: host data intact
    AMREX_ALWAYS_ASSERT(move_assign_cross.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(move_assign_cross.host_const()[0] == 70);
    AMREX_ALWAYS_ASSERT(move_assign_cross.host_const()[1] == 71);
    AMREX_ALWAYS_ASSERT(move_assign_cross.host_const()[2] == 72);

    // Move constructor: host data intact
    AMREX_ALWAYS_ASSERT(move_ctor_cross->host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(move_ctor_cross->host_const()[0] == 80);
    AMREX_ALWAYS_ASSERT(move_ctor_cross->host_const()[1] == 81);
    AMREX_ALWAYS_ASSERT(move_ctor_cross->host_const()[2] == 82);

    // Copy assignment (device-dirty): finalize synced D->H
    AMREX_ALWAYS_ASSERT(copy_assign_dirty_cross.host_const().size() == 3U);
    AMREX_ALWAYS_ASSERT(copy_assign_dirty_cross.host_const()[0] == 90);
    AMREX_ALWAYS_ASSERT(copy_assign_dirty_cross.host_const()[1] == 91);
    AMREX_ALWAYS_ASSERT(copy_assign_dirty_cross.host_const()[2] == 92);

    // --- Session 6: verify all five round-trip into a new session ---
    amrex::Initialize(argc, argv);
    {
        verify_host_device_match(copy_assign_cross);
        verify_host_device_match(*copy_ctor_cross);
        verify_host_device_match(move_assign_cross);
        verify_host_device_match(*move_ctor_cross);
        verify_host_device_match(copy_assign_dirty_cross);
    }
    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    std::cout << "TrackedVector tests passed.\n";
    return 0;
}
