// Benchmark for DenseBins::build.
//
// Times DenseBins::build for several distributions of items over bins and
// verifies each result. The distributions range from fully sorted by bin
// (adjacent items share a bin, the common case for particles, which are
// periodically sorted by cell) to uniformly random.
//
// Run with the default inputs, or override on the command line, e.g.
//   ./main3d.hip.HIP.ex inputs tests="sorted random" nitems=33554432 nbins=65536
//   ./main3d.hip.HIP.ex inputs tests=cells ppc=125

#include <AMReX.H>
#include <AMReX_DenseBins.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <string>

using namespace amrex;

namespace
{

// Every bin index is assigned exactly once, in increasing order of bin.
Vector<int> makeSorted (int nitems, int nbins)
{
    Vector<int> items(nitems);
    for (int i = 0; i < nitems; ++i) {
        items[i] = static_cast<int>((Long(i) * nbins) / nitems);
    }
    return items;
}

// Sorted, then a fraction of the items are moved to a random bin.
Vector<int> makeNearlySorted (int nitems, int nbins, Real frac_moved)
{
    Vector<int> items = makeSorted(nitems, nbins);
    for (auto& item : items) {
        if (Random() < frac_moved) { item = static_cast<int>(Random_int(nbins)); }
    }
    return items;
}

Vector<int> makeRandom (int nitems, int nbins)
{
    Vector<int> items(nitems);
    for (auto& item : items) { item = static_cast<int>(Random_int(nbins)); }
    return items;
}

// ppc items per cell of an ncell^3 grid, in cell order (x fastest), binned
// into tiles of tile[0] x tile[1] x tile[2] cells, as done by shared-memory
// particle deposition.
Vector<int> makeCells (int ncell, int ppc, const Vector<int>& tile, int& nbins)
{
    const int ntx = (ncell + tile[0] - 1) / tile[0];
    const int nty = (ncell + tile[1] - 1) / tile[1];
    const int ntz = (ncell + tile[2] - 1) / tile[2];
    nbins = ntx * nty * ntz;

    Vector<int> items;
    items.reserve(std::size_t(ncell) * ncell * ncell * ppc);
    for (int k = 0; k < ncell; ++k) {
    for (int j = 0; j < ncell; ++j) {
    for (int i = 0; i < ncell; ++i) {
        const int bin = (k / tile[2] * nty + j / tile[1]) * ntx + i / tile[0];
        for (int p = 0; p < ppc; ++p) { items.push_back(bin); }
    }}}
    return items;
}

// Check that the permutation is a valid bin sort of the items.
bool checkBins (const DenseBins<int>& bins, const Vector<int>& items, int nbins)
{
    const auto nitems = static_cast<int>(items.size());

    Vector<int> perm(nitems);
    Vector<int> offsets(nbins+1);
    Gpu::copyAsync(Gpu::deviceToHost, bins.permutationPtr(), bins.permutationPtr()+nitems,
                   perm.begin());
    Gpu::copyAsync(Gpu::deviceToHost, bins.offsetsPtr(), bins.offsetsPtr()+nbins+1,
                   offsets.begin());
    Gpu::streamSynchronize();

    if (offsets[0] != 0 || offsets[nbins] != nitems) { return false; }

    Vector<char> seen(nitems, 0);
    for (int b = 0; b < nbins; ++b) {
        if (offsets[b] > offsets[b+1]) { return false; }
        for (int j = offsets[b]; j < offsets[b+1]; ++j) {
            const int i = perm[j];
            if (i < 0 || i >= nitems || seen[i] || items[i] != b) { return false; }
            seen[i] = 1;
        }
    }
    return true;
}

struct Result
{
    double min_ms = std::numeric_limits<double>::max();
    double avg_ms = 0.;
    bool correct = false;
};

Result benchmark (const Vector<int>& items, int nbins, int nwarmup, int nrepeat)
{
    const auto nitems = static_cast<int>(items.size());
    Gpu::DeviceVector<int> items_d(nitems);
    Gpu::copyAsync(Gpu::hostToDevice, items.begin(), items.end(), items_d.begin());
    Gpu::streamSynchronize();
    const int* pitems = items_d.data();

    DenseBins<int> bins;
    auto build = [&] () {
        bins.build(BinPolicy::Default, nitems, pitems, nbins,
                   [=] AMREX_GPU_HOST_DEVICE (int j) noexcept -> unsigned int { return j; });
    };

    Result r;
    for (int n = 0; n < nwarmup; ++n) { build(); }
    for (int n = 0; n < nrepeat; ++n) {
        Gpu::streamSynchronize();
        const double t0 = amrex::second();
        build();
        Gpu::streamSynchronize();
        const double dt = (amrex::second() - t0) * 1.e3;
        r.min_ms = std::min(r.min_ms, dt);
        r.avg_ms += dt / nrepeat;
    }
    r.correct = checkBins(bins, items, nbins);
    return r;
}

void report (const std::string& name, int nitems, int nbins, const Result& r)
{
    amrex::Print() << std::left << std::setw(16) << name << std::right
                   << std::setw(12) << nitems
                   << std::setw(10) << nbins
                   << std::setw(12) << nitems / nbins
                   << std::fixed << std::setprecision(3)
                   << std::setw(11) << r.min_ms
                   << std::setw(11) << r.avg_ms
                   << std::setw(11) << nitems / r.min_ms * 1.e-6
                   << std::setw(9) << (r.correct ? "ok" : "FAIL") << "\n";
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int nitems = 1 << 24;
        int nbins = 4096;
        Real frac_moved = 0.05;
        int ncell = 64;
        int ppc = 8;
        Vector<int> tile = {6, 6, 8};
        int nwarmup = 2;
        int nrepeat = 20;
        ParmParse pp;
        pp.query("nitems", nitems);
        pp.query("nbins", nbins);
        pp.query("frac_moved", frac_moved);
        pp.query("ncell", ncell);
        pp.query("ppc", ppc);
        pp.queryarr("tile", tile);
        pp.query("nwarmup", nwarmup);
        pp.query("nrepeat", nrepeat);
        Vector<std::string> test_list;
        pp.queryarr("tests", test_list);
        if (test_list.empty()) {
            test_list = {"sorted", "nearly_sorted", "random", "cells"};
        }

        amrex::Print() << "DenseBins::build benchmark, " << nrepeat << " repetitions\n"
                       << std::left << std::setw(16) << "distribution" << std::right
                       << std::setw(12) << "nitems" << std::setw(10) << "nbins"
                       << std::setw(12) << "items/bin" << std::setw(11) << "min [ms]"
                       << std::setw(11) << "avg [ms]" << std::setw(11) << "Gitems/s"
                       << std::setw(9) << "check" << "\n";

        bool all_correct = true;
        for (const auto& test : test_list) {
            Vector<int> items;
            int nb = nbins;
            std::string name = test;
            if (test == "sorted") {
                items = makeSorted(nitems, nbins);
            } else if (test == "nearly_sorted") {
                items = makeNearlySorted(nitems, nbins, frac_moved);
            } else if (test == "random") {
                items = makeRandom(nitems, nbins);
            } else if (test == "cells") {
                items = makeCells(ncell, ppc, tile, nb);
                name = "cells_ppc" + std::to_string(ppc);
            } else {
                amrex::Abort("Unknown test " + test);
            }
            const Result r = benchmark(items, nb, nwarmup, nrepeat);
            report(name, static_cast<int>(items.size()), nb, r);
            all_correct = all_correct && r.correct;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(all_correct, "DenseBins produced wrong bins");
    }
    amrex::Finalize();
}
