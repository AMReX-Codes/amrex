// SPDX-FileCopyrightText: Copyright (c) 2025-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
// Device algorithm adapted from NVIDIA NCCL 74738463791127c33abef2e8f7da3a585e025d48:
// docs/examples/06_device_api/02_alltoall_gin/main.cu.
// Independent data generation and result verification live in nccl_gin.cpp (C++23). See NCCL-GIN-LICENSE.txt.
#include <cuda_runtime.h>
#include <nccl_device.h>
#include "nccl_gin.H"
#include <cstdio>
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

extern "C" void amrex_test_gin_initialize (int* argc, char*** argv)
{
    amrex::Initialize(*argc, *argv);
}

extern "C" void amrex_test_gin_finalize ()
{
    amrex::Finalize();
}

extern "C" MPI_Comm amrex_test_gin_communicator ()
{
    return amrex::ParallelDescriptor::Communicator();
}

__global__ void gin_alltoall (ncclWindow_t send, ncclWindow_t recv,
                             size_t count, ncclDevComm dev)
{
    ncclGin gin{dev, 0};
    uint64_t signal = gin.readSignal(0);
    ncclGinBarrierSession<ncclCoopCta> barrier{
        ncclCoopCta(), gin, ncclTeamTagWorld(), 0};
    barrier.sync(ncclCoopCta(), cuda::memory_order_acquire,
                 ncclGinFenceLevel::Relaxed);
    size_t bytes = count * sizeof(float);
    for (int peer = threadIdx.x; peer < dev.nRanks; peer += blockDim.x) {
        gin.put(ncclTeamWorld(dev), peer, recv, dev.rank * bytes,
                send, peer * bytes, bytes, ncclGin_SignalInc{0});
    }
    gin.waitSignal(ncclCoopCta(), 0, signal + dev.nRanks);
    gin.flush(ncclCoopCta());
    barrier.sync(ncclCoopCta(), cuda::memory_order_release,
                 ncclGinFenceLevel::Relaxed);
}

struct Fixture {
    ncclComm_t comm;
    size_t count;
    size_t bytes;
    void* send{};
    void* recv{};
    ncclWindow_t sendwin{};
    ncclWindow_t recvwin{};
    ncclDevComm dev{};
};

#define NC(call) do { auto rc = (call); if (rc != ncclSuccess) { \
  std::fprintf(stderr, "%s: %s\n", #call, ncclGetErrorString(rc)); return int(rc); } } while (0)
#define CU(call) do { auto rc = (call); if (rc != cudaSuccess) { \
  std::fprintf(stderr, "%s: %s\n", #call, cudaGetErrorString(rc)); return 10000 + int(rc); } } while (0)

extern "C" int amrex_test_gin_open (ncclComm_t comm, int count, void** output)
{
    auto* fixture = new Fixture{};
    fixture->comm = comm;
    fixture->count = count;
    fixture->bytes = 2 * count * sizeof(float);
    NC(ncclMemAlloc(&fixture->send, fixture->bytes));
    NC(ncclMemAlloc(&fixture->recv, fixture->bytes));
    NC(ncclCommWindowRegister(comm, fixture->send, fixture->bytes,
        &fixture->sendwin, NCCL_WIN_COLL_SYMMETRIC));
    NC(ncclCommWindowRegister(comm, fixture->recv, fixture->bytes,
        &fixture->recvwin, NCCL_WIN_COLL_SYMMETRIC));
    ncclDevCommRequirements requirements = NCCL_DEV_COMM_REQUIREMENTS_INITIALIZER;
    requirements.worldGinBarrierCount = 1;
    requirements.ginSignalCount = 1;
    requirements.ginConnectionType = NCCL_GIN_CONNECTION_FULL;
    NC(ncclDevCommCreate(comm, &requirements, &fixture->dev));
    std::puts("GIN_DEVICE_COMM_CREATED");
    std::fflush(stdout);
    *output = fixture;
    return 0;
}

extern "C" int amrex_test_gin_run (void* opaque, float const* input, float* output)
{
    auto* fixture = static_cast<Fixture*>(opaque);
    CU(cudaMemcpy(fixture->send, input, fixture->bytes, cudaMemcpyHostToDevice));
    gin_alltoall<<<1, 128>>>(fixture->sendwin, fixture->recvwin,
                            fixture->count, fixture->dev);
    CU(cudaGetLastError());
    CU(cudaDeviceSynchronize());
    CU(cudaMemcpy(output, fixture->recv, fixture->bytes, cudaMemcpyDeviceToHost));
    return 0;
}

extern "C" int amrex_test_gin_close (void* opaque)
{
    auto* fixture = static_cast<Fixture*>(opaque);
    NC(ncclDevCommDestroy(fixture->comm, &fixture->dev));
    NC(ncclCommWindowDeregister(fixture->comm, fixture->sendwin));
    NC(ncclCommWindowDeregister(fixture->comm, fixture->recvwin));
    NC(ncclMemFree(fixture->send));
    NC(ncclMemFree(fixture->recv));
    delete fixture;
    return 0;
}

extern "C" int amrex_test_gin_properties (ncclComm_t comm, int* support, int* teams, int* gin)
{
    ncclCommProperties_t props = NCCL_COMM_PROPERTIES_INITIALIZER;
    ncclResult_t result = ncclCommQueryProperties(comm, &props);
    if (result == ncclSuccess) {
        *support = props.deviceApiSupport;
        *teams = props.nLsaTeams;
        *gin = static_cast<int>(props.ginType);
    }
    return static_cast<int>(result);
}
