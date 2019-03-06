#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <algorithm>

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    if (narg == 0) {
        amrex::Print()
            << "\n"
            << " Usage:\n"
            << "      fnan plotfile\n"
            << "\n"
            << " Description:\n"
            << "      This program takes a single plotfile and reports for each\n"
            << "      variable whether there is a NaN."
            << std::endl;
        return;
    }

    const auto& fname = amrex::get_command_argument(1);
    PlotFileData plotfile(fname);
    const auto& names = plotfile.varNames();
    const int ncomp = plotfile.nComp();
    const int nlevels = plotfile.finestLevel() + 1;
    int nwidth = 0;
    for (auto const& name : names) {
        nwidth = std::max(nwidth, static_cast<int>(name.size()));
    }
    for (int n = 0; n < ncomp; ++n) {
        const std::string& varname = names[n];
        Vector<int> has_nan(nlevels);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            const MultiFab& mf = plotfile.get(ilev,varname);
            has_nan[ilev] = mf.contains_nan(0,1,0);
        }

        int num_nans = 0;
        for (int b : has_nan) {
            if (b) ++num_nans;
        }
        if (num_nans == 0) {
            amrex::Print() << " " << std::setw(nwidth+1) << std::left << varname << ": clean" << "\n";
        } else {
            amrex::Print() << " " << varname << ": has NaNs on level(s)";
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                if (has_nan[ilev]) {
                    amrex::Print() << "  " << ilev;
                }
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
