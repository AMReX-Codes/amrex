#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdlib>

using namespace amrex;

struct ErrZone {
    Real max_abs_err = std::numeric_limits<Real>::lowest();
    int level;
    int grid_index;
    IntVect cell;
};

int main_main()
{
    const int narg = amrex::command_argument_count();

    Real global_error = 0.0;
    bool any_nans = false;
    ErrZone err_zone;
    bool all_variables_found = true;

    // defaults
    int norm = 0;
    std::string plotfile_a;
    std::string plotfile_b;
    std::string diffvar;
    int zone_info = false;
    int allow_diff_grids = false;
    std::string zone_info_var_name;
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
        } else if (fname == "-z" or fname == "--zone_info") {
            zone_info_var_name = amrex::get_command_argument(++farg);
            zone_info = true;
        } else if (fname == "-d" or fname == "--diffvar") {
            diffvar = amrex::get_command_argument(++farg);
            plot_names[0] = diffvar;
        } else if (fname == "-a" or fname == "--allow_diff_grids") {
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
            << "    -g|--ghost            : compare the ghost cells too (if stored)\n"
            << "    -n|--norm num         : what norm to use (default is 0 for inf norm)\n"
            << "    -d|--diffvar var      : output a plotfile showing the differences for\n"
            << "                            variable var\n"
            << "    -z|--zone_info var    : output the information for a zone corresponding\n"
            << "                            to the maximum error for the given variable\n"
            << "    -a|--allow_diff_grids : allow different BoxArrays covering the same domain\n"
            << std::endl;
        return 0;
    }

    PlotFileData pf_a(plotfile_a);
    PlotFileData pf_b(plotfile_b);
    pf_b.syncDistributionMap(pf_a);

    const int dm = pf_a.spaceDim();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pf_a.spaceDim() == pf_b.spaceDim(),
                                     "ERROR: plotfiles have different numbers of spatial dimensions");

    const int finest_level = pf_a.finestLevel();
    const int nlevels = finest_level+1;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pf_a.finestLevel() == pf_b.finestLevel(),
                                     "ERROR: number of levels do not match");

    const int ncomp_a = pf_a.nComp();
    const int ncomp_b = pf_b.nComp();
    if (ncomp_a != ncomp_b) {
        amrex::Print() << "\n WARNING: number of variables do not match\n";
    }

    int save_var_a = -1;
    int zone_info_var_a = -1;

    const Vector<std::string>& names_a = pf_a.varNames();
    const Vector<std::string>& names_b = pf_b.varNames();

    Vector<int> ivar_b(ncomp_a,-1); // in case the variables are not in the same order
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

    for (int ilev = 0; ilev < nlevels; ++ilev) {
        const auto& dx_a = pf_a.cellSize(ilev);
        const auto& dx_b = pf_b.cellSize(ilev);
        bool not_match = AMREX_D_TERM(   dx_a[0] != dx_b[0],
                                      || dx_a[1] != dx_b[1],
                                      || dx_a[2] != dx_b[2] );
        if (not_match) {
            amrex::Abort("ERROR: grid dx does not match");
        }
    }

    // create a multifab to store the difference for output, if desired
    Vector<MultiFab> mf_array(nlevels);
    if (save_var_a >= 0) {
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            mf_array[ilev].define(pf_a.boxArray(ilev),
                                  pf_a.DistributionMap(ilev),
                                  1, 0);
        }
    }

    amrex::Print() << "\n"
                   << " " << std::setw(24) << std::right << "variable name"
                   << "  " << std::setw(24) << "absolute error"
                   << "  " << std::setw(24) << "relative error" << "\n"
                   << " " << std::setw(24) << " "
                   << "  " << std::setw(24) << "(||A - B||)"
                   << "  " << std::setw(24) << "(||A - B||/||A||)" << "\n"
                   << " " << std::string(76,'-') << "\n";

    // go level-by-level and patch-by-patch and compare the data
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        if (pf_a.boxArray(ilev).empty() && pf_b.boxArray(ilev).empty()) {
            continue;
        }
        bool grids_match = pf_a.boxArray(ilev) == pf_b.boxArray(ilev);
        if (!grids_match && !allow_diff_grids) {
            amrex::Abort("ERROR: grids do not match");
        } else if (!grids_match) {
            // do they cover the same domain?
            if (!pf_a.boxArray(ilev).contains(pf_b.boxArray(ilev)) ||
                !pf_a.boxArray(ilev).contains(pf_b.boxArray(ilev))) {
                amrex::Abort("ERROR: grids do not cover same domain");
            }
        }

        Vector<Real> aerror(ncomp_a, 0.0);
        Vector<Real> rerror(ncomp_a, 0.0);
        Vector<Real> rerror_denom(ncomp_a, 0.0);
        Vector<int> has_nan_a(ncomp_a, false);
        Vector<int> has_nan_b(ncomp_a, false);
        for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
            if (ivar_b[icomp_a] >= 0) {
                const MultiFab& mf_a = pf_a.get(ilev, names_a[icomp_a]);
                MultiFab mf_b;
                if (grids_match) {
                    mf_b = pf_b.get(ilev, names_b[ivar_b[icomp_a]]);
                } else {
                    mf_b.define(mf_a.boxArray(), mf_a.DistributionMap(), 1, 0);
                    MultiFab tmp = pf_b.get(ilev, names_b[ivar_b[icomp_a]]);
                    mf_b.ParallelCopy(tmp);
                }
                has_nan_a[icomp_a] = mf_a.contains_nan();
                has_nan_b[icomp_a] = mf_b.contains_nan();
                MultiFab::Subtract(mf_b,mf_a,0,0,1,0); // b = b - a
                Real max_err = mf_b.norm0();
                if (norm == 1) {
                    aerror[icomp_a] = mf_b.norm1();
                    rerror[icomp_a] = aerror[icomp_a];
                    rerror_denom[icomp_a] = mf_a.norm1();
                } else if (norm == 2) {
                    aerror[icomp_a] = mf_b.norm2();
                    rerror[icomp_a] = aerror[icomp_a];
                    rerror_denom[icomp_a] = mf_a.norm2();
                } else {
                    aerror[icomp_a] = max_err;
                    rerror[icomp_a] = aerror[icomp_a];
                    rerror_denom[icomp_a] = mf_a.norm0();
                }

                if (norm == 0) {
                    rerror[icomp_a] /= rerror_denom[icomp_a];
                } else {
                    const auto& dx = pf_a.cellSize(ilev);
                    Real dv = 1.0;
                    for (int idim = 0; idim < dm; ++idim) {
                        dv *= dx[idim];
                    }
                    aerror[icomp_a] *= std::pow(dv,1./static_cast<Real>(norm));
                    rerror[icomp_a] = rerror[icomp_a]/rerror_denom[icomp_a];
                }

                if (icomp_a == save_var_a or icomp_a == zone_info_var_a) {
                    mf_b.abs(0,1);
                }

                if (icomp_a == save_var_a) {
                    MultiFab::Copy(mf_array[ilev], mf_b, 0, 0, 1, 0);
                }

                if (icomp_a == zone_info_var_a) {
                    if (max_err > err_zone.max_abs_err) {
                        err_zone.max_abs_err = max_err;
                        err_zone.level = ilev;
                        err_zone.cell = mf_b.maxIndex(0);
                        auto isects = pf_a.boxArray(ilev).intersections
                            (Box(err_zone.cell,err_zone.cell), true, 0);
                        err_zone.grid_index = isects[0].first;
                    }
                }
            }
        }

        amrex::Print() << " level = " << ilev << "\n";
        for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
            if (ivar_b[icomp_a] < 0) {
                amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                               << "  " << std::setw(50)
                               << "< variable not present in both files > \n";
            } else if (has_nan_a[icomp_a] or has_nan_b[icomp_a]) {
                amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                               << "  " << std::setw(50)
                               << "< NaN present > \n";
            } else {
                Real aerr = 0., rerr = 0.;
                if (aerror[icomp_a] > 0.) {
                    aerr = std::min(std::max(aerror[icomp_a], 1.e-99), 1.e98);
                }
                if (rerror[icomp_a] > 0.) {
                    rerr = std::min(std::max(rerror[icomp_a], 1.e-99), 1.e98);
                }
                amrex::Print() << " " << std::setw(24) << std::left << names_a[icomp_a]
                               << std::right
                               << "  " << std::setw(24) << std::setprecision(10) << aerr
                               << "  " << std::setw(24) << std::setprecision(10) << rerr
                               << "\n";
            }
        }

        global_error = std::max(global_error,
                                *(std::max_element(aerror.begin(),
                                                   aerror.end())));
        for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
            any_nans = any_nans or has_nan_a[icomp_a] or has_nan_b[icomp_a];
        }
    }

    if (save_var_a >= 0) {
        Vector<Geometry> geom;
        Vector<int> levsteps;
        Vector<IntVect> rr;
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            RealBox rb(pf_a.probLo(), pf_a.probHi());
            Array<int,3> isper = {0,0,0};
            geom.emplace_back(pf_a.probDomain(ilev), &rb,
                              pf_a.coordSys(), isper.data());
            levsteps.push_back(pf_a.levelStep(ilev));
            rr.emplace_back(pf_a.refRatio(ilev));
        }
        WriteMultiLevelPlotfile ("diffs", nlevels, GetVecOfConstPtrs(mf_array),
                                 plot_names, geom, pf_a.time(),
                                 levsteps, rr);
    }

    if (zone_info) {
        if (err_zone.max_abs_err > 0.) {
            ParallelDescriptor::Barrier();
            const DistributionMapping& dmap = pf_a.DistributionMap(err_zone.level);
            bool owner_proc = ParallelDescriptor::MyProc() == dmap[err_zone.grid_index];

            if (owner_proc) {
                amrex::AllPrint() << std::endl
                                  << " maximum error in " << zone_info_var_name << "\n"
                                  << "   level = " << err_zone.level << " (i,j,k) = " << err_zone.cell << "\n";
            }

            for (int icomp_a = 0; icomp_a < ncomp_a; ++icomp_a) {
                const MultiFab& mf = pf_a.get(err_zone.level,names_a[icomp_a]);
                if (owner_proc) {
                    Real v = mf[err_zone.grid_index](err_zone.cell);
                    amrex::AllPrint() << " " << std::setw(24)
                                      << names_a[icomp_a] << "  "
                                      << std::setw(24) << std::right
                                      << v << "\n";
                }
            }
        }
    }

    if (global_error == 0.0 and !any_nans) {
        if (! all_variables_found) {
            amrex::Print() << " WARNING: not all variables present in both files\n";
        }
        amrex::Print() << " PLOTFILE AGREE" << std::endl;
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    int r = main_main();
    ignore_unused(r);
    amrex::Finalize();
#ifndef BL_USE_MPI
    return r;
#endif
}
