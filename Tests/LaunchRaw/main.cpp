#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_Reduce.H>

using namespace amrex;

void test3d () {

    const IntVectND<3> num_blocks {31, 23, 11};
#ifdef AMREX_USE_GPU
    static constexpr IntVectND<3> blockdim {2, 8, 16};
#else
    static constexpr IntVectND<3> blockdim {1, 1, 1};
#endif
    static constexpr int num_threads = blockdim[0] * blockdim[1] * blockdim[2];

    Gpu::DeviceVector<int> vect(num_threads * num_blocks[0] * num_blocks[1] * num_blocks[2], -999);

    auto data = vect.dataPtr();

    LaunchRaw<num_threads>(num_blocks,
        [=](auto lh){
            data[lh.globalIdx1D()] = lh.blockIdx1D();
        });

    LaunchRaw<num_threads>(num_blocks,
        [=](auto lh){
            data[lh.blockIdx1D() * lh.blockDim1D() + lh.threadIdx1D()] += lh.threadIdx1D();
        });

    LaunchRaw<num_threads>(num_blocks,
        [=](auto lh){
            auto block = lh.blockIdxND();
            auto thread = lh.template threadIdxND<blockdim[0], blockdim[1], blockdim[2]>();
            auto tmp = data[
                (block[0] + block[1] * num_blocks[0] +
                block[2] * num_blocks[0] * num_blocks[1]) * num_threads +
                thread[2] + thread[1] * blockdim[2] +
                thread[0] * blockdim[2] * blockdim[1]
            ];
            lh.syncthreads();
            data[lh.blockIdx1D() * lh.blockDim1D() + lh.threadIdx1D()] = tmp;
        });

    LaunchRaw<num_threads, int>(num_blocks, num_threads,
        [=](auto lh){
            auto smem = lh.shared_memory();
            auto thread1 = lh.template threadIdxND<blockdim[2], blockdim[1], blockdim[0]>();
            auto locid1 = thread1[2] + thread1[1] * blockdim[0] +
                          thread1[0] * blockdim[0] * blockdim[1];
            auto thread2 = lh.template threadIdxND<blockdim[2], blockdim[0], blockdim[1]>();
            auto locid2 = thread2[0] + thread2[2] * blockdim[2] +
                          thread2[1] * blockdim[2] * blockdim[1];
            smem[lh.threadIdx1D()] = data[lh.blockIdx1D() * lh.blockDim1D() + locid1];
            lh.syncthreads();
            data[lh.blockIdx1D() * lh.blockDim1D() + locid2] = smem[locid2];
        });

    LaunchRaw<num_threads>(num_blocks,
        [=](auto lh){
            data[lh.globalIdx1D()] =
                data[lh.globalIdx1D()] == static_cast<int>(lh.blockIdx1D() + lh.threadIdx1D());
        });

    AMREX_ALWAYS_ASSERT(Reduce::Sum(vect.size(), data, 0) == static_cast<int>(vect.size()));
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        test3d();

        amrex::Print() << "Passed! \n";
    }
    amrex::Finalize();
}
