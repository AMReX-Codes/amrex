#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    if (narg == 0) {
        amrex::Print()
            << "\n"
            << " Usage:\n"
            << "      fvarnames plotfile\n"
            << "\n"
            << " Description:\n"
            << "      This program takes a single plotfile and dumps out the list of variables\n"
            << std::endl;
        return;
    }

    const auto& fname = amrex::get_command_argument(1);
    PlotFileData plotfile(fname);
    const auto& names = plotfile.varNames();
    int n = 0;
    for (auto const& name : names) {
        amrex::Print() << std::setw(4) << n++ << "   " << name << "\n";
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
