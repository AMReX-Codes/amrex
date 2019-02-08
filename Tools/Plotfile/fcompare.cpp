#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotfileData.H>
#include <algorithm>
#include <limits>

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    // defaults
    int norm = 0;
    std::string plotfile_a;
    std::string plotfile_b;
    std::string diffvar;
    int do_ghost = false;
    int zone_info = false;
    int allow_diff_grids = false;
    std::string zone_info_var_name;
    Real global_error = 0.0;
    Vector<std::string> plot_names(1);

    int farg = 1;
    while (farg <= narg) {
        const std::string fname = amrex::get_command_argument(farg);
        if (fname == "--infile1") {
            plotfile_a = amrex::get_command_argument(++farg);
        } else if (fname == "--infile2") {
            plotfile_b = amrex::get_command_argument(++farg);
        } else if (fname == "-n" or fname == "--norm") {
            norm = std::stoi(amrex::get_command_argument(++farg));
        } else if (fname == "-g" or fname == "--ghost") {
            do_ghost = true;
        } else if (fname == "-z" or fname == "--zone_info") {
            zone_info_var_name = amrex::get_command_argument(++farg);
            zone_info = true;
        } else if (fname == "-d" or fname == "--diffvar") {
            diffvar = amrex::get_command_argument(++farg);
            plot_names[0] = diffvar;
        } else if (fname == "-a" or fname == "--allow_diff_ba") {
            allow_diff_grids = true;
        } else {
            break;
        }
        ++farg;
    };

    if (plotfile_a.empty()) {
        plotfile_a = amrex::get_command_argument(farg++);
    }
    if (plotfile_b.empty()) {
        plotfile_b = amrex::get_command_argument(farg++);
    }

    if (plotfile_a.empty() and plotfile_b.empty()) {
        amrex::Print()
            << "\n"
            << " Compare two plotfiles, zone by zone, to machine precision\n"
            << " and report the maximum absolute and relative errors for each\n"
            << " variable.\n"
            << "\n"
            << " usage:\n"
            << "    fcompare [-g|--ghost] [-n|--norm num] [-d|--diffvar var] [-z|--zone_info var] [-a|--allow_diff_grids] file1 file2\n"
            << "\n"
            << " optional arguments:\n"
            << "    -g|--ghost         : compare the ghost cells too (if stored)\n"
            << "    -n|--norm num      : what norm to use (default is 0 for inf norm)\n"
            << "    -d|--diffvar var   : output a plotfile showing the differences for\n"
            << "                         variable var\n"
            << "    -z|--zone_info var : output the information for a zone corresponding\n"
            << "                         to the maximum error for the given variable\n"
            << "    -a|--allow_diff_ba : allow different BoxArrays covering the same domain\n"
            << std::endl;
        return;
    }

    PlotfileData pf_a(plotfile_a);
    PlotfileData pf_b(plotfile_b);

    const int dm = pf_a.spaceDim();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pf_a.spaceDim() == pf_b.spaceDim(),
                                     "ERROR: plotfiles have different numbers of spatial dimensions");

    const int finest_level = pf_a.finestLevel();
    const int nlevels = finest_level+1;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pf_a.finestLevel() == pf_b.finestLevel(),
                                     "ERROR: number of levels do not match");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pf_a.boxArray(finest_level) == pf_b.boxArray(finest_level),
                                     "ERROR: grids do not match");

    const int ncomp_a = pf_a.nComp();
    const int ncomp_b = pf_b.nComp();
    if (ncomp_a != ncomp_b) {
        amrex::Print() << "\n WARNING: number of variables do not match\n";
    }

    bool all_variables_found = true;

    int save_var_a = -1;
    int zone_info_var_a = -1;

    const Vector<std::string>& names_a = pf_a.varNames();
    const Vector<std::string>& names_b = pf_b.varNames();

    Vector<int> ivar_b(ncomp_a); // in case the variables are not in the same order
    for (int n_a = 0; n_a < ncomp_a; ++n_a) {
        auto r = std::find(std::begin(names_b), std::end(names_b), names_a[n_a]);
        if (r == std::end(names_b)) {
            amrex::Print() << " WARNING: variable " << names_a[n_a] << " not found in plotfile 2\n";
            all_variables_found = false;
        } else {
            ivar_b[n_a] = std::distance(std::begin(names_b), r);
        }

        if (names_a[n_a] == diffvar) {
            save_var_a = n_a;
        }

        if (names_a[n_a] == zone_info_var_name) {
            zone_info_var_a = n_a;
        }
    }

    // also print out, as a diagnostic, those variables in plotfile 1 that
    // are not in plotfile 2
    for (int n_b = 0; n_b < ncomp_b; ++n_b) {
        auto r = std::find(std::begin(names_a),std::end(names_a),names_b[n_b]);
        if (r == std::end(names_a)) {
            amrex::Print() << " WARNING: variable " << names_b[n_b] << " not found in plotfile 1\n";
            all_variables_found = false;            
        }
    }

    // create a multifab to store the difference for output, if desired
    Vector<MultiFab> mf_array(nlevels);
    if (save_var_a > 0) {
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            mf_array[ilev].define(pf_a.boxArray(ilev),
                                  pf_a.DistributionMap(ilev),
                                  1, 0);
        }
    }

    // go level-by-level and patch-by-patch and compare the data
    bool any_nans = false;
    bool gc_warn = false; // move this above xxxxx
    Real max_abs_err = std::numeric_limits<Real>::lowest();
    Vector<Real> aerror(ncomp_a);
    Vector<Real> rerror(ncomp_a);
    Vector<Real> rerror_denom(ncomp_a);
    Vector<int> has_nan_a(ncomp_a,false);
    Vector<int> has_nan_b(ncomp_a,false);
    for (int icomp = 0; icomp < ncomp_a; ++icomp) {
#if 0
        const std::string 
        for (int ilev = 0; ilev < nlevels; ++ilev) {
        }
        plotfile_a.flush(names_a[icomp]);
        plotfile_b.flush();
#endif
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
