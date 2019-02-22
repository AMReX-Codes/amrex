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
    int ndir_pass = 3;
    bool origin = false;
    Array<Real,3> location = {std::numeric_limits<Real>::lowest(),
                              std::numeric_limits<Real>::lowest(),
                              std::numeric_limits<Real>::lowest()};
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
        if (name == "-n" or name == "--normaldir") {
            ndir_pass = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-p" or name == "--palette") {
            pfname[0] = amrex::get_command_argument(++farg);
        } else if (name == "-v" or name == "--variable") {
            compname = amrex::get_command_argument(++farg);
        } else if (name == "-M" or name == "--max") {
            def_mx = std::stod(amrex::get_command_argument(++farg));
            ldef_mx = true;
        } else if (name == "-m" or name == "--min") {
            def_mn = std::stod(amrex::get_command_argument(++farg));
            ldef_mn = true;
        } else if (name == "-L" or name == "--max_level") {
            max_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-l" or name == "--log") {
            do_log = true;
        } else if (name == "-g" or name == "--origin") {
            origin = true;
        } else if (name == "-c" or name == "--coordinates") {
            location[0] = std::stod(amrex::get_command_argument(++farg));
            location[1] = std::stod(amrex::get_command_argument(++farg));
            location[2] = std::stod(amrex::get_command_argument(++farg));
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
            << " produce an image of 2-d plotfile or slice of 3-d plotfile\n"
            << "\n"
            << " usage: fsnapshot [args] plotfile\n"
            << " args [-v|--variable] name       : variable to plot (default: density)\n"
            << "      [-p|--palette] pfname      : use the file pfname as the Palette (default ~/amrvis.Palette)\n"
            << "      -m val                     : set the minimum value of the data to val\n"
            << "      -M val                     : set the maximum value of the data to val\n"
            << "      [-L|--max_level] n         : max fine level to get data from (default: finest)\n"
            << "      [-l|--log]                 : toggle log plot\n"
            << "      [-n|--normaldir] {0,1,2,3} : direction normal to slice. (default: 3, i.e., all directions)\n"
            << "                                   This option is for 3d plotfile only.\n"
            << "      [-g|--origin]              : slice through origin (i.e., lower-left corner) or center (default: center)\n"
            << "                                   This option is for 3d plotfile only.\n"
            << "      [-c|--coordinates] x y z   : coordinates on slice.  If specified, this will override the -g option\n"
            << "                                   For slices of all three directions, this is the intersection of the three\n"
            << "                                   planes.  Otherwise, the transverse coordinates are not used.\n"
            << "                                   If specified, all three coordinates must be provided.\n"
            << "                                   This option is for 3d plotfile only.\n"
            << "\n";
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

    if (dim == 1) {
        amrex::Abort("ERROR: This is a 1D pltfile");
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

    const Box& finedomainbox = pf.probDomain(max_level);
    const auto flo = amrex::lbound(finedomainbox);
    const auto fhi = amrex::ubound(finedomainbox);
    int iloc[3];
    if (origin) {
        iloc[0] = flo.x;
        iloc[1] = flo.y;
        iloc[2] = flo.z;
    } else {
        iloc[0] = (fhi.x-flo.x+1)/2 + flo.x;
        iloc[1] = (fhi.y-flo.y+1)/2 + flo.y;
        iloc[2] = (fhi.z-flo.z+1)/2 + flo.z;
    }

    if (location[0] > -1.e36 or location[1] > -1.e36 or location[2] > -1.e36) {
        Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
        Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(max_level);
        for (int idim = 0; idim < dim; ++idim) {
            iloc[idim] = (location[idim]-problo[idim]) / dx[idim];
        }
    }

    iloc[0] = std::max(flo.x,std::min(fhi.x,iloc[0]));
    iloc[1] = std::max(flo.y,std::min(fhi.y,iloc[1]));
    iloc[2] = std::max(flo.z,std::min(fhi.z,iloc[2]));

    int ndir_begin, ndir_end, ndirs;
    if (dim == 2) {
        ndir_begin = 2;
        ndir_end = 3;
        ndirs = 1;
    } else if (ndir_pass == 3) {
        ndir_begin = 0;
        ndir_end = 3;
        ndirs = 3;
    } else {
        ndir_begin = ndir_pass;
        ndir_end = ndir_begin+1;
        ndirs = 1;
    }

    Vector<Box> finebox(3, finedomainbox);
    Vector<MultiFab> datamf(3);
    for (int idir = ndir_begin; idir < ndir_end; ++idir) {
        finebox[idir].setSmall(idir, iloc[idir]);
        finebox[idir].setBig(idir, iloc[idir]);
        BoxArray ba(finebox[idir]);
        datamf[idir].define(ba, DistributionMapping{ba}, 1, 0);
    }

    Vector<int> rr(max_level+1,1);
    for (int ilev = max_level-1; ilev >= 0; --ilev) {
        rr[ilev] = rr[ilev+1] * pf.refRatio(ilev);
    }

    Real gmx = std::numeric_limits<Real>::lowest();
    Real gmn = std::numeric_limits<Real>::max();

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        const MultiFab& pltmf = pf.get(ilev, compname);
        gmx = std::max(gmx, pltmf.max(0));
        gmn = std::min(gmn, pltmf.min(0));
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
                for (int idir = ndir_begin; idir < ndir_end; ++idir) {
                    const Box& crsebox = amrex::coarsen(finebox[idir], rrlev);
                    const Box& ibox = bx & crsebox;
                    if (ibox.ok()) {
                        const auto& data = datamf[idir].array(0); // there is only one box
                        IntVect rrslice = rrlev;
                        rrslice[idir] = 1;
                        amrex::For(ibox, [=] (int i, int j, int k)
                        {
                            if (m(i,j,k) == 0) { // not covered by fine
                                const Real d = plt(i,j,k);
                                for         (int koff = 0; koff < rrslice[2]; ++koff) {
                                    int kk = k*rrlev[2] + koff;
                                    for     (int joff = 0; joff < rrslice[1]; ++joff) {
                                        int jj = j*rrlev[1] + joff;
                                        for (int ioff = 0; ioff < rrslice[0]; ++ioff) {
                                            int ii = i*rrlev[0] + ioff;
                                            data(ii,jj,kk) = d;
                                        }
                                    }
                                }
                            }
                        });
                    }
                }
            }
        } else {
            for (MFIter mfi(pltmf); mfi.isValid(); ++mfi) {
                const auto& plt = pltmf.array(mfi);
                const Box& bx = mfi.validbox();
                for (int idir = ndir_begin; idir < ndir_end; ++idir) {
                    const Box& ibox = bx & finebox[idir];
                    if (ibox.ok()) {
                        const auto& data = datamf[idir].array(0); // there is only one box
                        amrex::ParallelFor(ibox, [=] (int i, int j, int k)
                        {
                            data(i,j,k) = plt(i,j,k);
                        });
                    }
                }
            }
        }
    }

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

    BaseFab<unsigned char> intdat;
    for (int idir = ndir_begin; idir < ndir_end; ++idir) {
        intdat.resize(finebox[idir],1);
        const int width = (idir == 0) ? finebox[idir].length(1) : finebox[idir].length(0);
        const int height = (idir == 2) ? finebox[idir].length(1) : finebox[idir].length(2);
        const auto& intarr = intdat.array();
        const auto& realarr = datamf[idir].array(0);
        Real fac = 253.999 / (gmx-gmn);
        amrex::ParallelFor(finebox[idir], [=] (int i, int j, int k)
        {
            int jj = (idir == 2) ? height - 1 - j : j;  // flip the data in second image direction
            int kk = (idir == 2) ? k : height - 1 - k;
            Real rd = realarr(i,jj,kk);
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
        if (dim == 2) {
            ppmfile += "." + compname + ".ppm";
        } else {
            if (idir == 0) {
                ppmfile += "." + compname + ".YZ" + ".ppm";
            } else if (idir == 1) {
                ppmfile += "." + compname + ".XZ" + ".ppm";
            } else {
                ppmfile += "." + compname + ".XY" + ".ppm";
            }
        }

        storePPM(ppmfile, intdat.dataPtr(), width, height, r, g, b);
    }
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
