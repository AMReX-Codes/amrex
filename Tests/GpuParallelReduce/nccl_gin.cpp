#include "nccl_gin.H"
#include <cuda_runtime.h>
#include <mpi.h>
#include <nccl.h>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static_assert(__cplusplus >= 202100L, "Compile the numerical reference as C++23");
static_assert(NCCL_VERSION_CODE >= 23004, "The GIN test requires NCCL 2.30.4 or newer");

namespace {
void fail (char const* operation)
{
    std::fprintf(stderr, "FAILED: %s\n", operation);
    MPI_Abort(amrex_test_gin_communicator(), 1);
    std::abort();
}
void check_mpi (int error, char const* operation)
{
    if (error != MPI_SUCCESS) { fail(operation); }
}
void check_cuda (cudaError_t error, char const* operation)
{
    if (error != cudaSuccess) {
        std::fprintf(stderr, "%s\n", cudaGetErrorString(error));
        fail(operation);
    }
}
void check_nccl (ncclResult_t error, char const* operation)
{
    if (error != ncclSuccess) {
        std::fprintf(stderr, "%s\n", ncclGetErrorString(error));
        fail(operation);
    }
}
}

int main (int argc, char** argv)
{
    amrex_test_gin_initialize(&argc, &argv);
    int rank = 0, ranks = 0;
    check_mpi(MPI_Comm_rank(amrex_test_gin_communicator(), &rank), "rank");
    check_mpi(MPI_Comm_size(amrex_test_gin_communicator(), &ranks), "size");
    if (ranks != 2) { fail("exactly two MPI ranks required"); }
    int selected = 0;
    check_cuda(cudaGetDevice(&selected), "AMReX-selected GPU");
    cudaDeviceProp device{};
    check_cuda(cudaGetDeviceProperties(&device, selected), "device properties");
    std::array<cudaUUID_t, 2> uuids{};
    check_mpi(MPI_Allgather(&device.uuid, static_cast<int>(sizeof(device.uuid)), MPI_BYTE,
                     uuids.data(), static_cast<int>(sizeof(device.uuid)), MPI_BYTE,
                     amrex_test_gin_communicator()), "gather GPU UUIDs");
    if (std::memcmp(&uuids[0], &uuids[1], sizeof(cudaUUID_t)) == 0) {
        fail("two distinct physical GPU UUIDs required");
    }
    std::array<char, MPI_MAX_PROCESSOR_NAME> hostname{};
    int length = 0;
    check_mpi(MPI_Get_processor_name(hostname.data(), &length), "hostname");
    std::array<std::array<char, MPI_MAX_PROCESSOR_NAME>, 2> hosts{};
    check_mpi(MPI_Allgather(hostname.data(), static_cast<int>(hostname.size()), MPI_CHAR,
                     hosts.data(), static_cast<int>(hostname.size()), MPI_CHAR,
                     amrex_test_gin_communicator()), "gather hostnames");
    if (std::strcmp(hosts[0].data(), hosts[1].data()) == 0) {
        fail("two distinct execution nodes required");
    }
    std::printf("rank=%d host=%s GPU=%s uuid=", rank, hostname.data(), device.name);
    for (unsigned char byte : device.uuid.bytes) { std::printf("%02x", byte); }
    std::printf(" CXX=%ld\n", static_cast<long>(__cplusplus));
    std::fflush(stdout);

    int version = 0;
    check_nccl(ncclGetVersion(&version), "NCCL version");
    if (version < 23004) { fail("NCCL 2.30.4 or newer is required"); }
    ncclUniqueId id{};
    if (rank == 0) { check_nccl(ncclGetUniqueId(&id), "unique ID"); }
    check_mpi(MPI_Bcast(&id, static_cast<int>(sizeof(id)), MPI_BYTE, 0, amrex_test_gin_communicator()), "broadcast ID");
    ncclComm_t comm{};
    check_nccl(ncclCommInitRank(&comm, ranks, id, rank), "initialize NCCL");
    int support = 0, teams = 0, gin = 0;
    check_nccl(static_cast<ncclResult_t>(amrex_test_gin_properties(comm, &support, &teams, &gin)),
         "query Device API capabilities");
    std::printf("rank=%d NCCL=%d deviceApiSupport=%d nLsaTeams=%d ginType=%d\n",
                rank, version, support, teams, gin);
    std::fflush(stdout);

    if (!support || gin == 0) { fail("Device API and GIN must both be supported"); }
    constexpr int count = 32;
    std::array<float, 2 * count> input{}, output{};
    void* fixture = nullptr;
    if (amrex_test_gin_open(comm, count, &fixture) != 0) { fail("create GIN fixture"); }
    for (int iteration = 0; iteration < 2; ++iteration) {
        for (int destination = 0; destination < ranks; ++destination) {
            for (int i = 0; i < count; ++i) {
                input[destination * count + i] = static_cast<float>(
                    rank * 1000 + destination * 100 + i + iteration * 10000);
            }
        }
        if (amrex_test_gin_run(fixture, input.data(), output.data()) != 0) { fail("GIN device kernel"); }
        for (int source = 0; source < ranks; ++source) {
            for (int i = 0; i < count; ++i) {
                float expected = static_cast<float>(
                    source * 1000 + rank * 100 + i + iteration * 10000);
                if (output[source * count + i] != expected) { fail("exact C++23 GIN result check"); }
            }
        }
        std::printf("GIN_VALIDATED rank=%d iteration=%d src0_first=%.0f src1_first=%.0f src0_last=%.0f src1_last=%.0f\n",
                    rank, iteration, output[0], output[count], output[count-1], output[2*count-1]);
        std::fflush(stdout);
    }
    if (amrex_test_gin_close(fixture) != 0) { fail("close GIN fixture"); }
    check_nccl(ncclCommFinalize(comm), "finalize NCCL");
    check_nccl(ncclCommDestroy(comm), "destroy NCCL");
    check_mpi(MPI_Barrier(amrex_test_gin_communicator()), "final barrier");
    if (rank == 0) {
        std::puts("NCCL_GIN_DEVICE_KERNEL_TWO_DISTINCT_GPUS_TWO_NODES_CXX23_TWO_ITERATIONS_EXACT_PASS");
    }
    amrex_test_gin_finalize();
}
