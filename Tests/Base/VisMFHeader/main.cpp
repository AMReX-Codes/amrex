#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <bit>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <type_traits>

//
// Writes MultiFabs whose min and max in the VisMF header are inf, tiny or the
// largest numbers, reads them back, and checks the data and the header values.
//
// Fast math lets the compiler assume there is no inf or NaN.  So the special
// values are bit patterns that depend on a run-time value, and they become
// floating-point numbers in memory only, as in Tests/Base/NaN.
//

using namespace amrex;

namespace {

using U = std::conditional_t<sizeof(Real) == 8, std::uint64_t, std::uint32_t>;

struct Bits
{
    U inf, minf, nan, tiny, mtiny, big, mbig;
};

// The values depend only on the cell, so a re-read MultiFab can be checked on
// any layout. There is one box for every 8 cells in x.
void fill (FArrayBox& fab, int pattern, Bits const& b)
{
    auto const& a = fab.array();
    amrex::LoopOnCpu(fab.box(), fab.nComp(), [&] (int i, int j, int k, int n)
    {
        U v = std::bit_cast<U>(Real(1.) + Real(0.5)*Real(i) - Real(2*n));
        if (pattern == 1) {
            // Every min and max in the header is tiny.
            v = (n == 0) ? b.tiny : b.mtiny;
        } else {
            const int ibox = i/8;
            if (ibox == 1) {
                v = b.inf;
            } else if (ibox == 2) {
                v = b.minf;
            } else if (ibox == 3 && i == 28) {
                v = b.inf;
            } else if (ibox == 4) {
                v = b.nan;
            } else if (ibox == 5) {
                v = (n == 0) ? b.tiny : b.mtiny;
            } else if (ibox == 6) {
                v = (n == 0) ? b.big : b.mbig;
            }
        }
        std::memcpy(&a(i,j,k,n), &v, sizeof(Real));
    });
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        Box domain(IntVect(0), IntVect(AMREX_D_DECL(63,7,7)));
        BoxArray ba(domain);
        ba.maxSize(8);
        DistributionMapping dm(ba);
        const int ncomp = 2;

        const std::string dir = "vismf_header_test";
        amrex::UtilCreateCleanDirectory(dir, true);

        constexpr int nman = std::numeric_limits<Real>::digits - 1;
        constexpr U sign = U(1) << (sizeof(Real)*8-1);
        constexpr U expo = ((U(1) << (sizeof(Real)*8-1-nman)) - 1) << nman;
        volatile U vzero = 0;
        const U zero = vzero;

        // Use a subnormal number, or the smallest normal number if the build
        // flushes subnormal numbers to zero (e.g., with fast math).
        volatile Real vtiny = std::numeric_limits<Real>::denorm_min();
        volatile Real vone = Real(1.);
        volatile Real vprod = vtiny*vone;
        const Real prod = vprod;
        const U tiny = (std::bit_cast<U>(prod) == U(0)) ? (U(1) << nman) : U(1);

        const Bits b{.inf = expo + zero, .minf = (sign | expo) + zero,
                     .nan = (expo | (U(1) << (nman-1))) + zero,
                     .tiny = tiny + zero, .mtiny = (sign | tiny) + zero,
                     .big = (expo - U(1)) + zero, .mbig = (sign | (expo - U(1))) + zero};

        Long nfail = 0;
        auto check = [&] (Real const& a, U bits)
        {
            if (std::bit_cast<U>(a) != bits) { ++nfail; }
        };

        // Write and read back with every header version.
        for (int version = 1; version <= 4; ++version) {
            VisMF::SetHeaderVersion(static_cast<VisMF::Header::Version>(version));
            for (int pattern = 0; pattern < 2; ++pattern) {
                MultiFab mf(ba, dm, ncomp, 0);
                for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    FArrayBox host(mfi.fabbox(), ncomp, The_Pinned_Arena());
                    fill(host, pattern, b);
                    Gpu::htod_memcpy(mf[mfi].dataPtr(), host.dataPtr(), host.nBytes());
                }

                const std::string name = dir + "/mf_v" + std::to_string(version)
                    + "_p" + std::to_string(pattern);
                VisMF::Write(mf, name);

                MultiFab mf2;
                VisMF::Read(mf2, name);
                AMREX_ALWAYS_ASSERT(mf2.boxArray() == ba && mf2.nComp() == ncomp
                                    && mf2.nGrow() == 0);

                for (MFIter mfi(mf2); mfi.isValid(); ++mfi) {
                    FArrayBox expected(mfi.fabbox(), ncomp, The_Pinned_Arena());
                    fill(expected, pattern, b);
                    FArrayBox host(mfi.fabbox(), ncomp, The_Pinned_Arena());
                    Gpu::dtoh_memcpy(host.dataPtr(), mf2[mfi].dataPtr(), host.nBytes());
                    if (std::memcmp(host.dataPtr(), expected.dataPtr(), host.nBytes()) != 0) {
                        ++nfail;
                    }
                }

                // The min and max read from the header
                VisMF vmf(name);
                if (version == VisMF::Header::NoFabHeaderFAMinMax_v1) {
                    if (pattern == 0) {
                        check(vmf.min(0), b.minf);
                        check(vmf.max(0), b.inf);
                    } else {
                        check(vmf.min(0), b.tiny);
                        check(vmf.max(0), b.tiny);
                        check(vmf.min(1), b.mtiny);
                        check(vmf.max(1), b.mtiny);
                    }
                } else if (version != VisMF::Header::NoFabHeader_v1) {
                    for (int ifab = 0; ifab < static_cast<int>(ba.size()); ++ifab) {
                        const int ibox = ba[ifab].smallEnd(0)/8;
                        if (pattern == 1) {
                            check(vmf.min(ifab,0), b.tiny);
                            check(vmf.max(ifab,0), b.tiny);
                            check(vmf.min(ifab,1), b.mtiny);
                            check(vmf.max(ifab,1), b.mtiny);
                        } else if (ibox == 1) {
                            check(vmf.max(ifab,0), b.inf);
                        } else if (ibox == 2) {
                            check(vmf.min(ifab,0), b.minf);
                        } else if (ibox == 3) {
                            check(vmf.min(ifab,0), std::bit_cast<U>(Real(13.)));
                            check(vmf.max(ifab,0), b.inf);
                        } else if (ibox == 5) {
                            check(vmf.min(ifab,0), b.tiny);
                            check(vmf.max(ifab,0), b.tiny);
                            check(vmf.min(ifab,1), b.mtiny);
                        } else if (ibox == 6) {
                            check(vmf.max(ifab,0), b.big);
                            check(vmf.min(ifab,1), b.mbig);
                        }
                    }
                }
            }
        }

        // A header edited by hand. In the min of the first and the last box,
        // it has nan, a blank before the comma, which takes the reader's
        // second path, and a number beyond the range of float.
        {
            if (ParallelDescriptor::IOProcessor()) {
                std::ifstream ifs(dir + "/mf_v1_p0_H");
                std::string hdr((std::istreambuf_iterator<char>(ifs)),
                                std::istreambuf_iterator<char>());
                const std::string from[3] = {"1.00000000000000000e+00,",
                                             "-1.00000000000000000e+00,",
                                             "2.90000000000000000e+01,"};
                const std::string to[3] = {"nan,",
                                           "-1.00000000000000000e+00 ,",
                                           "1.00000000000000000e+300,"};
                for (int i = 0; i < 3; ++i) {
                    auto pos = hdr.find(from[i]);
                    AMREX_ALWAYS_ASSERT(pos != std::string::npos);
                    hdr.replace(pos, from[i].size(), to[i]);
                }
                std::ofstream ofs(dir + "/edited_H");
                ofs << hdr;
            }
            ParallelDescriptor::Barrier();

            VisMF vmf(dir + "/edited");
            for (int ifab = 0; ifab < static_cast<int>(ba.size()); ++ifab) {
                const int ibox = ba[ifab].smallEnd(0)/8;
                if (ibox == 0) {
                    if (!amrex::isnan(vmf.min(ifab,0))) { ++nfail; }
                    check(vmf.min(ifab,1), std::bit_cast<U>(Real(-1.)));
                } else if (ibox == 7) {
#ifdef BL_USE_FLOAT
                    check(vmf.min(ifab,0), b.inf);
#else
                    check(vmf.min(ifab,0), std::bit_cast<U>(Real(1.e300)));
#endif
                }
            }
        }

        ParallelDescriptor::ReduceLongSum(nfail);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfail == 0,
            "VisMF::Read did not return the data written by VisMF::Write");
        amrex::Print() << "VisMF header test passed\n";
    }
    amrex::Finalize();
}
