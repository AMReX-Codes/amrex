#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

#include <cuda_runtime.h>
#include <nccl.h>
#include <nccl_device.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr int reused_iterations = 2;

void check_cuda (cudaError_t error, char const* what)
{
    if (error != cudaSuccess) {
        std::fprintf(stderr, "%s: %s\n", what, cudaGetErrorString(error));
        amrex::Abort(what);
    }
}

void check_nccl (ncclResult_t result, char const* what)
{
    if (result != ncclSuccess) {
        std::fprintf(stderr, "%s: %s\n", what, ncclGetErrorString(result));
        amrex::Abort(what);
    }
}

void check_mpi (int error, char const* what)
{
    if (error != MPI_SUCCESS) {
        amrex::Abort(what);
    }
}

__global__ void seed_kernel (float* data, int count, int rank)
{
    for (int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < count; i += blockDim.x * gridDim.x) {
        data[i] = static_cast<float>(rank + i);
    }
}

__global__ void allreduce_kernel (ncclWindow_t send_window,
                                  ncclWindow_t recv_window,
                                  int count,
                                  ncclDevComm dev_comm)
{
    ncclLsaBarrierSession<ncclCoopCta> barrier(
        ncclCoopCta(), dev_comm, ncclTeamTagLsa(), blockIdx.x);
    barrier.sync(ncclCoopCta(), cuda::memory_order_acquire);

    int const rank = dev_comm.rank;
    int const nranks = dev_comm.nRanks;
    int const global_thread = threadIdx.x + blockDim.x *
        (rank + blockIdx.x * nranks);
    int const global_threads = blockDim.x * gridDim.x * nranks;

    for (int i = global_thread; i < count; i += global_threads) {
        float value = 0.0f;
        for (int peer = 0; peer < nranks; ++peer) {
            auto* peer_data = static_cast<float*>(
                ncclGetLsaPointer(send_window, 0, peer));
            value += peer_data[i];
        }
        for (int peer = 0; peer < nranks; ++peer) {
            auto* peer_data = static_cast<float*>(
                ncclGetLsaPointer(recv_window, 0, peer));
            peer_data[i] = value;
        }
    }

    barrier.sync(ncclCoopCta(), cuda::memory_order_release);
}

class NCCLCommCache {
public:
    ~NCCLCommCache ()
    {
        for (auto entry = m_entries.rbegin(); entry != m_entries.rend(); ++entry) {
            check_nccl(ncclCommFinalize(entry->nccl), "ncclCommFinalize");
            check_nccl(ncclCommDestroy(entry->nccl), "ncclCommDestroy");
            check_mpi(MPI_Comm_free(&entry->mpi), "MPI_Comm_free");
        }
    }

    ncclComm_t get (MPI_Comm comm)
    {
        for (auto const& entry : m_entries) {
            int relation = MPI_UNEQUAL;
            check_mpi(MPI_Comm_compare(comm, entry.mpi, &relation),
                      "MPI_Comm_compare");
            if (relation == MPI_IDENT || relation == MPI_CONGRUENT) {
                return entry.nccl;
            }
        }

        Entry entry;
        check_mpi(MPI_Comm_dup(comm, &entry.mpi), "MPI_Comm_dup");

        int rank = 0;
        int nranks = 0;
        check_mpi(MPI_Comm_rank(entry.mpi, &rank), "MPI_Comm_rank");
        check_mpi(MPI_Comm_size(entry.mpi, &nranks), "MPI_Comm_size");

        ncclUniqueId id;
        if (rank == 0) {
            check_nccl(ncclGetUniqueId(&id), "ncclGetUniqueId");
        }
        check_mpi(MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, entry.mpi),
                  "MPI_Bcast(ncclUniqueId)");
        check_nccl(ncclCommInitRank(&entry.nccl, nranks, id, rank),
                   "ncclCommInitRank");

        m_entries.push_back(entry);
        return m_entries.back().nccl;
    }

private:
    struct Entry {
        MPI_Comm mpi = MPI_COMM_NULL;
        ncclComm_t nccl = nullptr;
    };

    std::vector<Entry> m_entries;
};

class NCCLDeviceResources {
public:
    NCCLDeviceResources (ncclComm_t comm, int count, int blocks)
        : m_comm(comm), m_count(count), m_blocks(blocks)
    {
        check_nccl(ncclMemAlloc(&m_send_data, m_count * sizeof(float)),
                   "ncclMemAlloc(send)");
        check_nccl(ncclMemAlloc(&m_recv_data, m_count * sizeof(float)),
                   "ncclMemAlloc(recv)");
        check_nccl(ncclCommWindowRegister(
            m_comm, m_send_data, m_count * sizeof(float), &m_send_window,
            NCCL_WIN_COLL_SYMMETRIC), "ncclCommWindowRegister(send)");
        check_nccl(ncclCommWindowRegister(
            m_comm, m_recv_data, m_count * sizeof(float), &m_recv_window,
            NCCL_WIN_COLL_SYMMETRIC), "ncclCommWindowRegister(recv)");

        ncclDevCommRequirements requirements =
            NCCL_DEV_COMM_REQUIREMENTS_INITIALIZER;
        requirements.lsaBarrierCount = m_blocks;
        check_nccl(ncclDevCommCreate(m_comm, &requirements, &m_dev_comm),
                   "ncclDevCommCreate");
    }

    ~NCCLDeviceResources ()
    {
        check_nccl(ncclDevCommDestroy(m_comm, &m_dev_comm),
                   "ncclDevCommDestroy");
        check_nccl(ncclCommWindowDeregister(m_comm, m_send_window),
                   "ncclCommWindowDeregister(send)");
        check_nccl(ncclCommWindowDeregister(m_comm, m_recv_window),
                   "ncclCommWindowDeregister(recv)");
        check_nccl(ncclMemFree(m_send_data), "ncclMemFree(send)");
        check_nccl(ncclMemFree(m_recv_data), "ncclMemFree(recv)");
    }

    void allreduce (amrex::Gpu::DeviceVector<float>& values,
                    std::vector<float>& host, int rank, int iteration)
    {
        seed_kernel<<<m_blocks, threads>>>(values.data(), m_count,
                                            rank + iteration);
        check_cuda(cudaGetLastError(), "seed_kernel launch");
        check_cuda(cudaDeviceSynchronize(), "seed_kernel synchronize");
        check_cuda(cudaMemcpy(m_send_data, values.data(),
                              m_count * sizeof(float),
                              cudaMemcpyDeviceToDevice),
                   "cudaMemcpy(input scratch)");

        allreduce_kernel<<<m_blocks, threads>>>(m_send_window, m_recv_window,
                                                m_count, m_dev_comm);
        check_cuda(cudaGetLastError(), "allreduce_kernel launch");
        check_cuda(cudaDeviceSynchronize(), "allreduce_kernel synchronize");
        check_cuda(cudaMemcpy(values.data(), m_recv_data,
                              m_count * sizeof(float),
                              cudaMemcpyDeviceToDevice),
                   "cudaMemcpy(output scratch)");
        check_cuda(cudaMemcpy(host.data(), values.data(),
                              m_count * sizeof(float),
                              cudaMemcpyDeviceToHost),
                   "cudaMemcpy(result)");
    }

private:
    static constexpr int threads = 128;
    ncclComm_t m_comm;
    int m_count;
    int m_blocks;
    void* m_send_data = nullptr;
    void* m_recv_data = nullptr;
    ncclWindow_t m_send_window = nullptr;
    ncclWindow_t m_recv_window = nullptr;
    ncclDevComm m_dev_comm{};
};

int local_rank (MPI_Comm comm)
{
    MPI_Comm local_comm = MPI_COMM_NULL;
    int rank = 0;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &local_comm);
    MPI_Comm_rank(local_comm, &rank);
    MPI_Comm_free(&local_comm);
    return rank;
}

} // namespace

int main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);

    int result = 0;
    {
        MPI_Comm comm = amrex::ParallelDescriptor::Communicator();
        int rank = amrex::ParallelDescriptor::MyProc(comm);
        int nranks = amrex::ParallelDescriptor::NProcs(comm);
        int device = local_rank(comm);

        int device_count = 0;
        check_cuda(cudaGetDeviceCount(&device_count), "cudaGetDeviceCount");
        if (device >= device_count) {
            std::fprintf(stderr, "rank %d selected GPU %d, only %d visible\n",
                         rank, device, device_count);
            result = 1;
        } else {
            check_cuda(cudaSetDevice(device), "cudaSetDevice");

            NCCLCommCache comm_cache;
            ncclComm_t nccl_comm = comm_cache.get(comm);
            if (comm_cache.get(comm) != nccl_comm) {
                std::fprintf(stderr,
                    "rank %d: equivalent MPI communicator did not hit NCCL cache\n",
                    rank);
                result = 4;
            }

            ncclCommProperties_t props = NCCL_COMM_PROPERTIES_INITIALIZER;
            check_nccl(ncclCommQueryProperties(nccl_comm, &props),
                       "ncclCommQueryProperties");
            if (!props.deviceApiSupport || props.nLsaTeams != 1) {
                if (rank == 0) {
                    std::fprintf(stderr,
                        "NCCL Device API unavailable: support=%d lsa_teams=%d\n",
                        props.deviceApiSupport, props.nLsaTeams);
                }
                result = 2;
            } else {
                constexpr int count = 32;
                constexpr int blocks = 1;
                std::vector<float> host(count);
                amrex::Gpu::DeviceVector<float> values(count);
                NCCLDeviceResources resources(nccl_comm, count, blocks);

                float const expected_rank_sum =
                    static_cast<float>(nranks * (nranks - 1) / 2);
                for (int iteration = 0; iteration < reused_iterations; ++iteration) {
                    resources.allreduce(values, host, rank, iteration);

                    for (int i = 0; i < count; ++i) {
                        float const expected = expected_rank_sum +
                            nranks * (i + iteration);
                        if (std::fabs(host[i] - expected) > 1.0e-6f) {
                            std::fprintf(stderr,
                                "rank %d iteration %d: result[%d]=%g expected=%g\n",
                                rank, iteration, i, host[i], expected);
                            result = 3;
                            break;
                        }
                    }
                }
            }

        }

        MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MAX, comm);
        if (rank == 0) {
            if (result == 0) {
                std::printf("NCCL Device API AllReduce: passed on %d rank(s), "
                            "with %d reused iterations.\n", nranks,
                            reused_iterations);
            } else {
                std::printf("NCCL Device API AllReduce: failed (%d).\n", result);
            }
        }
    }

    amrex::Finalize();
    return result;
}
