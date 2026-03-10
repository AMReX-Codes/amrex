#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_Reduce.H>

using namespace amrex;

void test3d () {

    constexpr int num_threads = 256;
    IntVectND<3> num_blocks {31, 23, 11};

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
            auto thread = lh.template threadIdxND<2, 8, 16>();
            auto tmp = data[
                (block[0] + block[1] * num_blocks[0] +
                block[2] * num_blocks[0] * num_blocks[1]) * num_threads +
                thread[2] + thread[1] * 16 +
                thread[0] * 16 * 8
            ];
            lh.synctheads();
            data[lh.blockIdx1D() * lh.blockDim1D() + lh.threadIdx1D()] = tmp;
        });

    LaunchRaw<num_threads, int>(num_blocks, num_threads,
        [=](auto lh){
            auto smem = lh.shared_memory();
            auto thread1 = lh.template threadIdxND<16, 8, 2>();
            auto locid1 = thread1[2] + thread1[1] * 2 + thread1[0] * 2 * 8;
            auto thread2 = lh.template threadIdxND<16, 2, 8>();
            auto locid2 = thread2[0] + thread2[2] * 16 + thread2[1] * 16 * 8;

            smem[lh.threadIdx1D()] = data[lh.blockIdx1D() * lh.blockDim1D() + locid1];
            lh.syncthreads();
            data[lh.blockIdx1D() * lh.blockDim1D() + locid2] = smem[locid2];
        });

    LaunchRaw<num_threads>(num_blocks,
        [=](auto lh){
            data[lh.globalIdx1D()] = data[lh.globalIdx1D()] == (lh.blockIdx1D() + lh.threadIdx1D());
        });

    AMREX_ALWAYS_ASSERT(Reduce::Sum(vect.size(), data, 0) == vect.size());
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
