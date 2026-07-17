#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

#include <memory>

using namespace amrex;

namespace
{
    void PrintUsage ()
    {
        amrex::Print()
            << "\n"
            << " Convert one or more native AMReX plotfiles to HDF5 plotfiles.\n"
            << "\n"
            << " usage:\n"
            << "    fconvert [-c|--compression descriptor] [-o|--output output] plotfile [plotfile ...]\n"
            << "\n"
            << " optional arguments:\n"
            << "    -c|--compression descriptor : HDF5 compression descriptor.\n"
            << "                                  Default: ZFP_ACCURACY@0.001\n"
            << "                                  Examples: None@0, ZLIB@5, ZFP_RATE@8,\n"
            << "                                            ZFP_PRECISION@12,\n"
            << "                                            ZFP_ACCURACY@0.001,\n"
            << "                                            ZFP_REVERSIBLE@1\n"
            << "    -o|--output output         : output name for one plotfile, or output prefix\n"
            << "                                  for multiple plotfiles\n"
            << '\n';
    }

    [[nodiscard]] std::string BaseName (std::string const& path)
    {
        auto pos = path.find_last_of("/\\");
        if (pos == std::string::npos) {
            return path;
        } else {
            return path.substr(pos+1);
        }
    }

#ifdef AMREX_USE_HDF5
    [[nodiscard]] std::string OutputName (std::string const& input,
                                          std::string const& output,
                                          bool multiple_inputs)
    {
        if (output.empty()) {
            return input;
        }

        if (!multiple_inputs) {
            return output;
        }

        if (output.back() == '/' || output.back() == '\\') {
            return output + BaseName(input);
        } else {
            return output + "_" + BaseName(input);
        }
    }
#endif
}

void main_main ()
{
    const int narg = amrex::command_argument_count();

    if (narg == 0) {
        PrintUsage();
        return;
    }

    std::string compression = "ZFP_ACCURACY@0.001";
    std::string output_name;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-h" || name == "--help") {
            PrintUsage();
            return;
        } else if (name == "-c" || name == "--compression") {
            if (farg == narg) {
                amrex::Abort("fconvert requires an argument after --compression");
            }
            compression = amrex::get_command_argument(++farg);
        } else if (name == "-o" || name == "--output") {
            if (farg == narg) {
                amrex::Abort("fconvert requires an argument after --output");
            }
            output_name = amrex::get_command_argument(++farg);
        } else {
            break;
        }
        ++farg;
    }

    if (farg > narg) {
        PrintUsage();
        return;
    }

#ifndef AMREX_USE_HDF5
    amrex::Abort("fconvert requires AMReX to be built with AMReX_HDF5=ON. "
                 "Build with AMReX_HDF5_ZFP=ON as well to use ZFP_* "
                 "compression descriptors.");
#else
    Vector<std::string> infiles;
    infiles.reserve(narg-farg+1);
    for (int i = farg; i <= narg; ++i) {
        infiles.push_back(amrex::get_command_argument(i));
    }

    bool const multiple_inputs = infiles.size() > 1;

    for (auto const& infile : infiles) {
        auto const start_time = amrex::second();

        PlotFileData plotfile(infile);
        int const finest_level = plotfile.finestLevel();
        int const nlevels = finest_level + 1;

        Vector<std::unique_ptr<MultiFab> > mf(nlevels);
        Vector<const MultiFab*> mf_ptrs(nlevels);
        Vector<Geometry> geom;
        Vector<int> level_steps(nlevels);
        Vector<IntVect> ref_ratio(finest_level);

        geom.reserve(nlevels);
        RealBox rb(plotfile.probLo(), plotfile.probHi());
        Array<int,AMREX_SPACEDIM> is_per{AMREX_D_DECL(0,0,0)};

        for (int ilev = 0; ilev < nlevels; ++ilev) {
            mf[ilev] = std::make_unique<MultiFab>(plotfile.get(ilev));
            mf_ptrs[ilev] = mf[ilev].get();
            geom.push_back(Geometry(plotfile.probDomain(ilev), rb,
                                    plotfile.coordSys(), is_per));
            level_steps[ilev] = plotfile.levelStep(ilev);
            if (ilev < finest_level) {
                ref_ratio[ilev] = plotfile.refRatioVect(ilev);
            }
        }

        auto const outfile = OutputName(infile, output_name, multiple_inputs);

        amrex::Print() << "fconvert: converting " << infile
                       << " -> " << outfile << ".h5"
                       << " with compression " << compression << '\n';

        WriteMultiLevelPlotfileHDF5(outfile, nlevels, mf_ptrs,
                                    plotfile.varNames(), geom, plotfile.time(),
                                    level_steps, ref_ratio, compression);

        amrex::Print() << "fconvert: finished " << outfile << ".h5 in "
                       << amrex::second() - start_time << " seconds\n";
    }
#endif
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
