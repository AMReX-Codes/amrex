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
    {
        BL_PROFILE("build");
        bins.buildOpenMP(nitems, items.data(), nbins, [=] (int j) noexcept -> int { return j; });
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    testDenseBins();

    amrex::Finalize();
}
