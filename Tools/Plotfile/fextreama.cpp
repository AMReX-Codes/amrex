#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <iterator>

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string varnames_arg;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-v" or name == "--variable") {
            varnames_arg = amrex::get_command_argument(++farg);
        } else {
            break;
        }
        ++farg;
    }

    if (farg > narg) {
        amrex::Print() << "\n"
                       << " Report the extrema (min/max) for each variable in a plotfile\n"
                       << " usage: \n"
                       << "    fextrema {[-v|--variable] name} plotfiles\n"
                       << "\n"
                       << "   -v names    : output information only for specified variables, given\n"
                       << "                 as a space-spearated string\n"
                       << std::endl;
        return;
    }

    int ntime = narg - farg + 1;
    Vector<std::string> var_names;

    for (int f = 0; f < ntime; ++f) {
        const std::string& filename = amrex::get_command_argument(f+farg);
        PlotFileData pf(filename);
        Vector<std::string> const& var_names_pf = pf.varNames();

        if (f == 0) {
            if (varnames_arg.empty()) { // all variables
                var_names = var_names_pf;
            } else {
                std::istringstream is(varnames_arg);
                var_names.assign(std::istream_iterator<std::string>{is},
                                 std::istream_iterator<std::string>{  });
            }
        }

        // get the extrema
        Vector<Real> vvmin(var_names.size());
        Vector<Real> vvmax(var_names.size());

        const int dim = pf.spaceDim();

        for (int ilev = pf.finestLevel(); ilev >= 0; --ilev) {
            if (ilev == pf.finestLevel()) {
                for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                    const MultiFab& mf = pf.get(ilev, var_names[ivar]);
                    vvmin[ivar] = mf.min(0,0,false);
                    vvmax[ivar] = mf.max(0,0,false);
                }
            } else {
                IntVect ratio{pf.refRatio(ilev)};
                for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                    ratio[idim] = 1;
                }
                iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                              pf.boxArray(ilev+1), ratio);
                for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                    const MultiFab& mf = pf.get(ilev, var_names[ivar]);
                    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        const auto lo = amrex::lbound(bx);
                        const auto hi = amrex::ubound(bx);
                        const auto& ifab = mask.array(mfi);
                        const auto& fab = mf.array(mfi);
                        for         (int k = lo.z; k <= hi.z; ++k) {
                            for     (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {
                                    if (ifab(i,j,k) == 0) {
                                        vvmin[ivar] = std::min(fab(i,j,k),vvmin[ivar]);
                                        vvmax[ivar] = std::max(fab(i,j,k),vvmax[ivar]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        ParallelDescriptor::ReduceRealMin(vvmin.data(), vvmin.size());
        ParallelDescriptor::ReduceRealMax(vvmax.data(), vvmax.size());

        if (ntime == 1) {
            amrex::Print() << " plotfile = " << filename << "\n"
                           << " time = " << std::setprecision(17) << pf.time() << "\n"
                           << " " << std::setw(22) << "variables"
                           << " " << std::setw(22) << "minimum value"
                           << " " << std::setw(22) << "maximum value"
                           << "\n";
            for (int i = 0; i < var_names.size(); ++i) {
                amrex::Print().SetPrecision(11)
                    << " " << std::setw(22) << std::left << var_names[i]
                    << " " << std::setw(22) << std::right << vvmin[i]
                    << " " << std::setw(22) << std::right << vvmax[i]
                    << "\n";
            }
            amrex::Print() << std::endl;
        } else {
            if (f == 0) {
                amrex::Print() << "# " << std::setw(23) << std::right << "time ";
                for (int i = 0; i < var_names.size(); ++i) {
                    amrex::Print() << "|" << std::setw(44) << std::left << var_names[i];
                }
                amrex::Print() << "|\n";
                amrex::Print() << "# " << std::setw(23) << " ";
                for (int i = 0; i < var_names.size(); ++i) {
                    amrex::Print() << "|" << std::setw(12) << std::right << "min" << std::setw(10) << " "
                                   << std::setw(12) << "max" << std::setw(10) << " ";
                }
                amrex::Print() << "|\n";
            }
            amrex::Print() << std::setw(23) << std::setprecision(11) << pf.time() << "   ";
            for  (int i = 0; i < var_names.size(); ++i) {
                amrex::Print() << std::setw(21) << std::setprecision(10) << vvmin[i]
                               << std::setw(21) << std::setprecision(10) << vvmax[i]
                               << "   ";
            }
            amrex::Print() << "\n";
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
