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
            << "      ftime plotfile\n"
            << "\n"
            << " Description:\n"
            << "      This program takes a whitespace-separated list of plotfiles and\n"
            << "      returns the time for each plotfile.\n"
            << std::endl;
        return;
    }

    for (int f = 1; f <= narg; ++f) {
        const auto& fname = amrex::get_command_argument(f);
        PlotFileData plotfile(fname);
        amrex::Print().SetPrecision(17) << fname << "    " << plotfile.time() << std::endl;
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
