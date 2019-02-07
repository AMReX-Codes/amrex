#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_AnyDimUtil.H>
#include <AMReX_PlotfileData.H>

using namespace amrex;

void main_main()
{
    bool b_full = false;
    bool b_gridfile = false;
    bool b_levels = false;

    const int narg = amrex::command_argument_count();

    int farg = 1;
    while (farg <= narg) {
        const auto fname = get_command_argument(farg);
        if (fname == "-f" or fname == "--full") {
            b_full = true;
        } else if (fname == "-g" or fname == "--gridfile") {
            b_gridfile = true;
        } else if (fname == "-l" or fname == "--levels") {
            b_levels = true;
        } else {
            break;
        }
        ++farg;
    }

    if (b_gridfile and b_full) {
        amrex::Abort("ERROR: cannot specify both full and gridfile modes");
    }

    if (farg > narg) {
        amrex::Print() << "\n"
                       << " Dump out information about the AMR levels and boxes\n"
                       << " Works with 1-, 2-, or 3-d datasets.\n"
                       << "\n"
                       << " usage:\n"
                       << "    fboxinfo [-f|--full] plotfile\n"
                       << "\n"
                       << " args:\n"
                       << "    [-f|--full]     output detailed information about the boxes\n"
                       << "    [-g|--gridfile] output a gridfile for use with test_average\n"
                       << "    [-l|--levels]   just output the number of levels\n"
                       << std::endl;
        return;
    }

    for (int f = farg; f <= narg; ++f) {
        const auto& fname = amrex::get_command_argument(f);
        PlotfileData plotfile(fname);

        if (!b_gridfile and !b_levels) {
            amrex::Print() << " plotfile: " << fname << "\n";
        }

        const int dim = plotfile.spaceDim();
        const int nlevels = plotfile.finestLevel()+1;

        if (b_levels) {
            amrex::Print() << " " << nlevels << std::endl;
            continue;
        }

        if (!b_gridfile)
        {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                const long nboxes = plotfile.boxArray(ilev).size();
                const long ncells = plotfile.boxArray(ilev).numPts();
                const Box prob_domain = plotfile.probDomain(ilev);
                const Real ncells_domain = prob_domain.d_numPts();
                amrex::Print() << " level " << std::setw(3) << ilev
                               << ": number of boxes = " << std::setw(6) << nboxes
                               << ", volume = "
                               << std::fixed << std::setw(6) << std::setprecision(2)
                               << 100.*(ncells/ncells_domain) << "%\n";
                if (dim == 1) {
                    amrex::Print() << "          maximum zones =   "
                                   << std::setw(7) << prob_domain.length(0) << "\n";
                } else if (dim == 2) {
                    amrex::Print() << "          maximum zones =   "
                                   << std::setw(7) << prob_domain.length(0)
                                   << " x "
                                   << std::setw(7) << prob_domain.length(1) << "\n";
                } else {
                    amrex::Print() << "          maximum zones =   "
                                   << std::setw(7) << prob_domain.length(0)
                                   << " x "
                                   << std::setw(7) << prob_domain.length(1)
                                   << " x "
                                   << std::setw(7) << prob_domain.length(2) << "\n";
                }
                amrex::Print() << std::endl;
            }
        }

        if (b_full) {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                amrex::Print() << "\n  level " << ilev << "\n";
                const BoxArray& ba = plotfile.boxArray(ilev);
                const long nboxes = ba.size();
                for (long ibox = 0; ibox < nboxes; ++ibox) {
                    const Box& b = ba[ibox];
                    if (dim == 1) {
                        amrex::Print() << "   box " << std::setw(5) << ibox
                                       << ":  (" << std::setw(5) << b.smallEnd(0)
                                       << ")"
                                       << "   (" << std::setw(5) << b.bigEnd(0)
                                       << ")"
                                       << "\n";
                    } else if (dim == 2) {
                        amrex::Print() << "   box " << std::setw(5) << ibox
                                       << ":  (" << std::setw(5) << b.smallEnd(0)
                                       <<    "," << std::setw(5) << b.smallEnd(1)
                                       << ")"
                                       << "   (" << std::setw(5) << b.bigEnd(0)
                                       <<    "," << std::setw(5) << b.bigEnd(1)
                                       << ")"
                                       << "\n";
                    } else {
                        amrex::Print() << "   box " << std::setw(5) << ibox
                                       << ":  (" << std::setw(5) << b.smallEnd(0)
                                       <<    "," << std::setw(5) << b.smallEnd(1)
                                       <<    "," << std::setw(5) << b.smallEnd(2)
                                       << ")"
                                       << "   (" << std::setw(5) << b.bigEnd(0)
                                       <<    "," << std::setw(5) << b.bigEnd(1)
                                       <<    "," << std::setw(5) << b.bigEnd(2)
                                       << ")"
                                       << "\n";
                    }
                }
            }
        }

        if (b_gridfile) {
            amrex::Print() << " " << std::setw(2) << nlevels << "\n";
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                const BoxArray& ba = plotfile.boxArray(ilev);
                const long nboxes = ba.size();
                const Box prob_domain = plotfile.probDomain(ilev);
                amrex::Print() << "   " << BoxND(prob_domain,dim)
                               << "  " << nboxes << "\n";
                for (int ibox = 0; ibox < nboxes; ++ibox) {
                    amrex::Print() << "      "  << BoxND(ba[ibox],dim) << "\n";
                }
            }
        }

        amrex::Print() << std::endl;
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
