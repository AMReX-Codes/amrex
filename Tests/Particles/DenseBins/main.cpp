#include <AMReX.H>
#include <AMReX_DenseBins.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

using namespace amrex;

void checkAnswer (const amrex::DenseBins<int>& bins)
{
    BL_PROFILE("checkAnswer");
    const auto* const perm = bins.permutationPtr();
    const auto* const bins_ptr = bins.binsPtr();
    const auto* const offsets = bins.offsetsPtr();

#ifdef AMREX_USE_GPU
    amrex::ParallelFor(bins.numItems(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        AMREX_ALWAYS_ASSERT(bins_ptr[perm[i]] <= bins_ptr[perm[i+1]]);
    });
#else
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < bins.numItems()-1; ++i)
    {
        AMREX_ALWAYS_ASSERT(bins_ptr[perm[i]] <= bins_ptr[perm[i+1]]);
    }
#endif

#ifdef AMREX_USE_GPU
    amrex::ParallelFor(bins.numItems(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        auto start = offsets[i  ];
        auto stop  = offsets[i+1];
        if (start < stop) {
            for (auto j = start+1; j < stop; ++j)
            {
                AMREX_ALWAYS_ASSERT(bins_ptr[perm[start]] == bins_ptr[perm[j]]);
            }
        }
    });
#else
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < bins.numBins(); ++i) {
        auto start = offsets[i  ];
        auto stop  = offsets[i+1];
        if (start == stop) continue;
        for (auto j = start+1; j < stop; ++j)
        {
            AMREX_ALWAYS_ASSERT(bins_ptr[perm[start]] == bins_ptr[perm[j]]);
        }
    }
#endif
}

void testGPU (int nbins, const amrex::Vector<int>& items)
{
    // copy to device
    Gpu::DeviceVector<int> items_d(items.size());
    Gpu::copyAsync(Gpu::hostToDevice, items.begin(), items.end(), items_d.begin());
    Gpu::Device::streamSynchronize();

    amrex::DenseBins<int> bins;
    bins.build(BinPolicy::GPU, items_d.size(), items_d.data(), nbins, [=] AMREX_GPU_DEVICE (int j) noexcept -> unsigned int { return j ; });

    checkAnswer(bins);
}

void testOpenMP (int nbins, const amrex::Vector<int>& items)
{
    amrex::DenseBins<int> bins;
    bins.build(BinPolicy::OpenMP, items.size(), items.data(), nbins, [=] (int j) noexcept -> unsigned int { return j ; });

    checkAnswer(bins);
}

void testSerial (int nbins, const amrex::Vector<int>& items)
{
    amrex::DenseBins<int> bins;
    bins.build(BinPolicy::Serial, items.size(), items.data(), nbins, [=] (int j) noexcept -> unsigned int { return j ; });

    checkAnswer(bins);
}

void initData (int nbins, amrex::Vector<int>& items)
{
    BL_PROFILE("init");

    const auto nitems = int(items.size());

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nitems; ++i) { items[i] = int(amrex::Random_int(nbins)); }
}

void testDenseBins ()
{
    int nitems;
    int nbins;

    ParmParse pp;
    pp.get("nitems", nitems);
    pp.get("nbins" , nbins);

    amrex::Vector<int> items(nitems);
    initData(nbins, items);

#ifndef AMREX_USE_SYCL
    testSerial(nbins, items);
#endif
#ifdef AMREX_USE_OMP
    testOpenMP(nbins, items);
#endif
#ifdef AMREX_USE_GPU
    testGPU(nbins, items);
#endif
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    testDenseBins();

    amrex::Finalize();
}
