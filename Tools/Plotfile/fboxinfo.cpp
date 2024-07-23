#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

namespace {

}

void main_main()
{
    bool b_full = false;
    bool b_gridfile = false;
    bool b_levels = false;

    const int narg = amrex::command_argument_count();

    int farg = 1;
    while (farg <= narg) {
        const auto fname = get_command_argument(farg);
        if (fname == "-f" || fname == "--full") {
            b_full = true;
        } else if (fname == "-g" || fname == "--gridfile") {
            b_gridfile = true;
        } else if (fname == "-l" || fname == "--levels") {
            b_levels = true;
        } else {
            break;
        }
        ++farg;
    }

    if (b_gridfile && b_full) {
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
                       << '\n';
        return;
    }

    for (int f = farg; f <= narg; ++f) {
        const auto& fname = amrex::get_command_argument(f);
        PlotFileData plotfile(fname);

        if (!b_gridfile && !b_levels) {
            amrex::Print() << " plotfile: " << fname << "\n";
        }

        const int dim = plotfile.spaceDim();
        const int nlevels = plotfile.finestLevel()+1;

        if (b_levels) {
            amrex::Print() << " " << nlevels << '\n';
            continue;
        }

        if (!b_gridfile)
        {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                const Long nboxes = plotfile.boxArray(ilev).size();
                const Long ncells = plotfile.boxArray(ilev).numPts();
                const Box prob_domain = plotfile.probDomain(ilev);
                const auto ncells_domain = prob_domain.d_numPts();
                amrex::Print() << " level " << std::setw(3) << ilev
                               << ": number of boxes = " << std::setw(6) << nboxes
                               << ", volume = "
                               << std::fixed << std::setw(6) << std::setprecision(2)
                               << 100.*static_cast<double>(ncells)/ncells_domain << "%\n";
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
                amrex::Print() << '\n';
            }
        }

        if (b_full) {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                amrex::Print() << "\n  level " << ilev << "\n";
                const BoxArray& ba = plotfile.boxArray(ilev);
                const auto nboxes = static_cast<int>(ba.size());
                for (int ibox = 0; ibox < nboxes; ++ibox) {
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
                const Long nboxes = ba.size();
                const Box prob_domain = plotfile.probDomain(ilev);
                amrex::Print() << "   ";
                amrex::IntVect prob_domain_type = prob_domain.type();
                amrex::detail::box_write(amrex::Print(), prob_domain.smallEnd().begin(),
                                                         prob_domain.bigEnd().begin(),
                                                         prob_domain_type.begin(), dim);
                amrex::Print() << "  " << nboxes << "\n";
                for (int ibox = 0; ibox < nboxes; ++ibox) {
                    amrex::Print() << "      ";
                    amrex::IntVect ba_type = ba[ibox].type();
                    amrex::detail::box_write(amrex::Print(), ba[ibox].smallEnd().begin(),
                                                             ba[ibox].bigEnd().begin(),
                                                             ba_type.begin(), dim);
                    amrex::Print() << "\n";
                }
            }
        }

        amrex::Print() << '\n';
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
