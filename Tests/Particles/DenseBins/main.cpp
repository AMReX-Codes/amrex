#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DenseBins.H>

using namespace amrex;

void testDenseBins ()
{
    int nitems;
    int nbins;

    ParmParse pp;
    pp.get("nitems", nitems);
    pp.get("nbins" , nbins);

    amrex::Vector<int> items(nitems);
    {
        BL_PROFILE("init");
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < nitems; ++i) { items[i] = amrex::Random_int(nbins); }
    }

    amrex::DenseBins<int> bins;
    bins.buildOpenMP(nitems, items.data(), nbins, [=] (int j) noexcept -> unsigned int { return j ; });

    {
        BL_PROFILE("checkAnswer");
        const auto perm = bins.permutationPtr();
        const auto bins_ptr = bins.binsPtr();
        const auto offsets = bins.offsetsPtr();
        const auto counts = bins.countsPtr();
        for (int i = 0; i < nitems-1; ++i)
        {
            AMREX_ALWAYS_ASSERT(bins_ptr[perm[i]] <= bins_ptr[perm[i+1]]);
        }

        for (int i = 0; i < bins.numBins(); ++i) {
            auto start = offsets[i  ];
            auto stop  = offsets[i+1];
            AMREX_ALWAYS_ASSERT(counts[i] == stop - start);
            if (start == stop) continue;
            for (int j = start+1; j < stop; ++j)
            {
                AMREX_ALWAYS_ASSERT(bins_ptr[perm[start]] == bins_ptr[perm[j]]);
            }
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    testDenseBins();

    amrex::Finalize();
}
