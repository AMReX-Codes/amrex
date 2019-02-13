#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PPMUtil.H>
#include <limits>
#include <iterator>
#include <fstream>
#include <cstdlib>

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string home(std::getenv("HOME"));
    if (home.empty()) {
        amrex::Abort("Failed to get environment variable HOME is");
    }
    Array<std::string,3> pfname = {"", home+"/.amrvis.Palette", home+"/.amrvis.palette"};

    std::string pltfile;
    std::string compname = "density";
    int max_level = -1;
    bool ldef_mx = false;
    bool ldef_mn = false;
    bool do_log = false;
    Real def_mx = std::numeric_limits<Real>::lowest();
    Real def_mn = std::numeric_limits<Real>::max();

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-p" or name == "--palette") {
            pfname[0] = amrex::get_command_argument(++farg);
        } else if (name == "-v" or name == "--variable") {
            compname = amrex::get_command_argument(++farg);
        } else if (name == "-M" or name == "--max") {
            def_mx = std::stod(amrex::get_command_argument(++farg));
            ldef_mx = true;
        } else if (name == "-m" or name == "--min") {
            def_mn = std::stod(amrex::get_command_argument(++farg));
            ldef_mn = true;
        } else if (name == "--max_level") {
            max_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-l" or name == "--log") {
            do_log = true;
        } else {
            break;
        }
        ++farg;
    }

    if (pltfile.empty() and farg <= narg) {
        pltfile = amrex::get_command_argument(farg);
    }

    if (pltfile.empty() or compname.empty()) {
        amrex::Print()
            << "\n"
            << " produce an image of a 2-d dataset\n"
            << "\n"
            << " usage: fsnapshot2d [args] plotfile\n"
            << " args [-v|--variable] name  : variable to plot (default: density)\n"
            << "      [-p|--palette] pfname : use the file pfname as the Palette (default ~/amrvis.Palette)\n"
            << "      -m val                : set the minimum value of the data to val\n"
            << "      -M val                : set the maximum value of the data to val\n"
            << "      [-l|--log]            : toggle log plot\n";
        return;
    }

    // make sure we have valid options set
    if (do_log) {
        if (ldef_mx and def_mx < 0.) amrex::Abort("ERROR: log plot specified with negative maximum");
        if (ldef_mn and def_mn < 0.) amrex::Abort("ERROR: log plot specified with negative minimum");
    }

    // get the palette
    Array<unsigned char,256> r, g, b, a;
    int numElements = 0;
    for (int i = 0; i < 3; ++i) {
        if (!pfname[i].empty()) {
            bool exist = false;
            {
                std::ifstream ifs(pfname[i]);
                exist = ifs.good();
            }
            if (exist) {
                numElements = loadPalette(pfname[i], r, g, b, a);
                break;
            }
        }
    }

    PlotFileData pf(pltfile);
    int dim = pf.spaceDim();

    if (dim != 2) {
        amrex::Abort("ERROR: not a 2D pltfile");
    }

    if (max_level < 0) {
        max_level = pf.finestLevel();
    } else {
        if (max_level < 0 or max_level > pf.finestLevel()) {
            amrex::Abort("ERROR: specified level not allowed");
        }
    }

    const auto& var_names = pf.varNames();
    if (std::find(var_names.begin(), var_names.end(), compname) == var_names.end()) {
        amrex::Abort("ERROR: " + compname + " not found in pltfile " + pltfile);
    }

    const Box& finebox = pf.probDomain(max_level);
    BoxArray ba(finebox);
    MultiFab datamf(ba, DistributionMapping{ba}, 1, 0);

    Vector<int> rr(max_level+1,1);
    for (int ilev = max_level-1; ilev >= 0; --ilev) {
        rr[ilev] = rr[ilev+1] * pf.refRatio(ilev);
    }

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        const auto& data = datamf.array(0); // there is only one box
        const MultiFab& pltmf = pf.get(ilev, compname);
        if (ilev < max_level) {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pltmf, pf.boxArray(ilev+1), ratio);
            for (MFIter mfi(pltmf); mfi.isValid(); ++mfi) {
                const auto& m = mask.array(mfi);
                const auto& plt = pltmf.array(mfi);
                const Box& bx = mfi.validbox();
                IntVect rrlev {rr[ilev]};
                for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                    rrlev[idim] = 1;
                }
                amrex::For(bx, [=] (int i, int j, int k)
                {
                    if (m(i,j,k) == 0) { // not covered by fine
                        const Real d = plt(i,j,k);
                        for         (int koff = 0; koff < rrlev[2]; ++koff) {
                            int kk = k*rrlev[2] + koff;
                            for     (int joff = 0; joff < rrlev[1]; ++joff) {
                                int jj = j*rrlev[1] + joff;
                                for (int ioff = 0; ioff < rrlev[0]; ++ioff) {
                                    int ii = i*rrlev[0] + ioff;
                                    data(ii,jj,kk) = d;
                                }
                            }
                        }
                    }
                });
            }
        } else {
            for (MFIter mfi(pltmf); mfi.isValid(); ++mfi) {
                const auto& plt = pltmf.array(mfi);
                const Box& bx = mfi.validbox();
                amrex::ParallelFor(bx, [=] (int i, int j, int k)
                {
                    data(i,j,k) = plt(i,j,k);
                });
            }
        }
    }

    Real gmx = datamf.max(0);
    Real gmn = datamf.min(0);

    amrex::Print() << " plotfile variable maximum = " << gmx << "\n"
                   << " plotfile variable minimum = " << gmn << "\n";

    if (ldef_mx) {
        amrex::Print() << " resetting variable maximum to " << def_mx << "\n";
        gmx = def_mx;
    }

    if (ldef_mn) {
        amrex::Print() << " resetting variable minimum to " << def_mn << "\n";
        gmn = def_mn;
    }

    if (do_log) {
        gmn = std::log10(gmn);
        gmx = std::log10(gmx);
    }

    BaseFab<unsigned char> intdat(finebox,1);
    const int width = finebox.length(0);
    const int height = finebox.length(1);
    const auto& intarr = intdat.array();
    const auto& realarr = datamf.array(0);
    Real fac = 253.999 / (gmx-gmn);
    amrex::ParallelFor(finebox, [=] (int i, int j, int k)
    {
        int jj = height - 1 - j;
        Real rd = realarr(i,jj,k);
        if (do_log) rd = std::log10(rd);
        int id = std::max(0,std::min(255,static_cast<int>((rd-gmn)*fac)));
        unsigned char c = static_cast<unsigned char>(id);
        constexpr unsigned char cmn = static_cast<unsigned char>(1);  // avoid zero
        constexpr unsigned char cmx = static_cast<unsigned char>(255);
        intarr(i,j,k) = std::max(cmn,std::min(cmx,c));
    });

    std::string ppmfile = pltfile;
    if (ppmfile.back() == '/') {
        ppmfile.pop_back();
    }
    ppmfile += "." + compname + ".ppm";

    storePPM(ppmfile, intdat.dataPtr(), width, height, r, g, b);
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
#ifdef BL_USE_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int color = (rank == 0) ? 0 : 1;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    if (rank == 0) {
        // Only process 0 is doing the work
        amrex::Initialize(argc, argv, false, comm);
        main_main();
        amrex::Finalize();
    }
    MPI_Comm_free(&comm);
    MPI_Finalize();
#else
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
#endif
}
