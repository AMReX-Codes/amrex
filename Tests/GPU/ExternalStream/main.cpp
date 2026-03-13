#include <AMReX.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

using namespace amrex;

namespace {

#ifdef AMREX_USE_GPU

gpuStream_t make_external_stream ()
{
    gpuStream_t s;
    AMREX_HIP_OR_CUDA_OR_SYCL(
        AMREX_HIP_SAFE_CALL(hipStreamCreateWithFlags(&s, hipStreamNonBlocking)),
        AMREX_CUDA_SAFE_CALL(cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking)),
        s.queue = new sycl::queue(Gpu::Device::streamQueue()));
    return s;
}

void destroy_external_stream (gpuStream_t s)
{
#if defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaStreamDestroy(s));
#elif defined(AMREX_USE_HIP)
    AMREX_HIP_SAFE_CALL(hipStreamDestroy(s));
#elif defined(AMREX_USE_SYCL)
    if (s.queue != nullptr) {
        s.queue->wait_and_throw();
        delete s.queue;
    }
#endif
}

void test_external_stream_region (gpuStream_t external)
{
    constexpr int n = 64;
    Gpu::DeviceVector<int> device_data(n, -1);
    Vector<int> host(n, 0);

    {
        Gpu::ExternalGpuStreamRegion guard(
            external, Gpu::ExternalStreamSync::No);
        AMREX_ALWAYS_ASSERT(Gpu::Device::usingExternalStream());
        AMREX_ALWAYS_ASSERT(Gpu::Device::gpuStream() == external);
        AMREX_ALWAYS_ASSERT(Gpu::Device::numGpuStreams() == 1);

        auto ptr = device_data.dataPtr();
        ParallelFor(n, [ptr] AMREX_GPU_DEVICE (int i) noexcept {
            ptr[i] = 2*i;
        });
    }

    Gpu::streamSynchronize(external);

    Gpu::dtoh_memcpy(host.data(), device_data.dataPtr(), n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        AMREX_ALWAYS_ASSERT(host[i] == 2*i);
    }

    AMREX_ALWAYS_ASSERT(!Gpu::Device::usingExternalStream());

    bool cleared_before = Gpu::clearFreeAsyncBuffer();
    AMREX_ALWAYS_ASSERT(!cleared_before);

    {
        Gpu::ExternalGpuStreamRegion guard(
            external, Gpu::ExternalStreamSync::No);
        Gpu::AsyncVector<int> async_device_data(n);
        Gpu::dtod_memcpy_async(async_device_data.data(), device_data.data(),
                               n * sizeof(int));
    }
    // When The_Async_Arena is involved, the dtor of
    // Gpu::ExternalGpuStreamRegion forces a stream sync. Thus the stream
    // should not be active.
    AMREX_ALWAYS_ASSERT(!Gpu::isStreamActive(external));

    bool cleared_after = Gpu::clearFreeAsyncBuffer();
    AMREX_ALWAYS_ASSERT(!cleared_after);
}

void test_explicit_reset (gpuStream_t external)
{
    auto* arena = The_Async_Arena();
    void* tmp = arena->alloc(512);

    Gpu::setExternalGpuStream(external);
    AMREX_ALWAYS_ASSERT(Gpu::Device::usingExternalStream());
    AMREX_ALWAYS_ASSERT(Gpu::Device::gpuStream() == external);

    Gpu::freeAsync(arena, tmp);
    Gpu::resetExternalGpuStream(Gpu::ExternalStreamSync::No);
    AMREX_ALWAYS_ASSERT(!Gpu::Device::usingExternalStream());

    Gpu::streamSynchronize(external);
}

void test_nested_external_stream_region (gpuStream_t outer, gpuStream_t inner)
{
    Gpu::ExternalGpuStreamRegion outer_guard(
        outer, Gpu::ExternalStreamSync::No);
    AMREX_ALWAYS_ASSERT(Gpu::Device::gpuStream() == outer);

    {
        Gpu::ExternalGpuStreamRegion inner_guard(
            inner, Gpu::ExternalStreamSync::No);
        AMREX_ALWAYS_ASSERT(Gpu::Device::gpuStream() == inner);
    }

    AMREX_ALWAYS_ASSERT(Gpu::Device::usingExternalStream());
    AMREX_ALWAYS_ASSERT(Gpu::Device::gpuStream() == outer);
}

void test_mfiter_and_reducer (gpuStream_t external)
{
    constexpr int ncell = 16;
    Box domain(IntVect(AMREX_D_DECL(0,0,0)),
               IntVect(AMREX_D_DECL(ncell-1, ncell-1, ncell-1)));
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    iMultiFab mf(ba, dm, 1, 0);

    Long expected = 0;
    LoopOnCpu(domain, [&expected] (int i, int j, int k) noexcept {
        expected += Long(i + j + k);
    });

    Gpu::ExternalGpuStreamRegion guard(
        external, Gpu::ExternalStreamSync::Yes);
    MFItInfo info;
    info.DisableDeviceSync();
    for (MFIter mfi(mf, info); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k) = i + j + k;
        });
    }

    Long device_sum = mf.sum(0);
    AMREX_ALWAYS_ASSERT(expected == device_sum);
}

void run_external_stream_tests ()
{
    gpuStream_t external = make_external_stream();
    gpuStream_t nested_external = make_external_stream();
    test_external_stream_region(external);
    test_explicit_reset(external);
    test_nested_external_stream_region(external, nested_external);
    test_mfiter_and_reducer(external);
    destroy_external_stream(nested_external);
    destroy_external_stream(external);
    amrex::Print() << "External GPU stream override test completed.\n";
}

#endif // AMREX_USE_GPU

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
#ifdef AMREX_USE_GPU
        run_external_stream_tests();
#else
        amrex::Print() << "External GPU stream test requires AMReX built with GPU support.\n";
#endif
    }
    amrex::Finalize();
    return 0;
}
