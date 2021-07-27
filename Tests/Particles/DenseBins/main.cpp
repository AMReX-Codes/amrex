#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DenseBins.H>

using namespace amrex;

void checkAnswer (const amrex::DenseBins<int>& bins,
                  const amrex::Vector<int>& items)
{
    BL_PROFILE("checkAnswer");
    const auto perm = bins.permutationPtr();
    const auto bins_ptr = bins.binsPtr();
    const auto offsets = bins.offsetsPtr();
    for (int i = 0; i < bins.numItems()-1; ++i)
    {
        AMREX_ALWAYS_ASSERT(bins_ptr[perm[i]] <= bins_ptr[perm[i+1]]);
    }

    for (int i = 0; i < bins.numBins(); ++i) {
        auto start = offsets[i  ];
        auto stop  = offsets[i+1];
        if (start == stop) continue;
        for (int j = start+1; j < stop; ++j)
        {
            AMREX_ALWAYS_ASSERT(bins_ptr[perm[start]] == bins_ptr[perm[j]]);
        }
    }
}

void testOpenMP (int nbins, const amrex::Vector<int>& items)
{
    amrex::DenseBins<int> bins;
    bins.buildOpenMP(items.size(), items.data(), nbins, [=] (int j) noexcept -> unsigned int { return j ; });

    checkAnswer(bins, items);
}

void testSerial (int nbins, const amrex::Vector<int>& items)
{
    amrex::DenseBins<int> bins;
    bins.build(items.size(), items.data(), nbins, [=] (int j) noexcept -> unsigned int { return j ; });

    checkAnswer(bins, items);
}

void initData (int nbins, amrex::Vector<int>& items)
{
    BL_PROFILE("init");

    const int nitems = items.size();

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nitems; ++i) { items[i] = amrex::Random_int(nbins); }
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

    testSerial(nbins, items);
    testOpenMP(nbins, items);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    testDenseBins();

    amrex::Finalize();
}
