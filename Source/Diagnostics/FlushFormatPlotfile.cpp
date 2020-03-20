#include "FlushFormatPlotfile.H"
#include "WarpX.H"

#include <AMReX_buildInfo.H>

using namespace amrex;

void
FlushFormatPlotfile::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<const amrex::MultiFab*> mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    MultiParticleContainer& mpc, int nlev,
    const std::string prefix) const
{
    auto & warpx = WarpX::GetInstance();
    const auto step = iteration[0];
    const std::string& filename = amrex::Concatenate(prefix, step);
    amrex::Print() << "  Writing plotfile " << filename << "\n";

    Vector<std::string> rfs;
    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::Version_v1);
    amrex::WriteMultiLevelPlotfile(filename, nlev,
                                   mf,
                                   varnames, geom,
                                   time, iteration, warpx.refRatio(),
                                   "HyperCLaw-V1.1",
                                   "Level_",
                                   "Cell",
                                   rfs
                                   );

    mpc.WritePlotFile(filename);

    WriteJobInfo(filename);

    WriteWarpXHeader(filename);

    VisMF::SetHeaderVersion(current_version);
}

void
FlushFormatPlotfile::WriteJobInfo(const std::string& dir) const
{

    auto & warpx = WarpX::GetInstance();

    if (ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = dir;

        std::string PrettyLine = std::string(78, '=') + "\n";
//        std::string OtherLine = std::string(78, '-') + "\n";
//        std::string SkipSpace = std::string(8, ' ') + "\n";

        FullPathJobInfoFile += "/warpx_job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " WarpX Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

        jobInfoFile << "\n\n";

        // build information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
        jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
        jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
        jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
        jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
        jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
        jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
        jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

        jobInfoFile << "\n";

        const char* githash1 = buildInfoGetGitHash(1);
        const char* githash2 = buildInfoGetGitHash(2);
        const char* githash3 = buildInfoGetGitHash(3);
        if (strlen(githash1) > 0) {
            jobInfoFile << "WarpX  git describe: " << githash1 << "\n";
        }
        if (strlen(githash2) > 0) {
            jobInfoFile << "AMReX  git describe: " << githash2 << "\n";
        }
        if (strlen(githash3) > 0) {
            jobInfoFile << "PICSAR git describe: " << githash3 << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= warpx.finestLevel(); i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << warpx.boxArray(i).size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                jobInfoFile << warpx.Geom(i).Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << " Boundary conditions\n";

        jobInfoFile << "   -x: " << "interior" << "\n";
        jobInfoFile << "   +x: " << "interior" << "\n";
        if (AMREX_SPACEDIM >= 2) {
            jobInfoFile << "   -y: " << "interior" << "\n";
            jobInfoFile << "   +y: " << "interior" << "\n";
        }
        if (AMREX_SPACEDIM == 3) {
            jobInfoFile << "   -z: " << "interior" << "\n";
            jobInfoFile << "   +z: " << "interior" << "\n";
        }

        jobInfoFile << "\n\n";


        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
    }
}

void
FlushFormatPlotfile::WriteWarpXHeader(const std::string& name) const
{
    auto & warpx = WarpX::GetInstance();
    if (ParallelDescriptor::IOProcessor())
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(name + "/WarpXHeader");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << "Checkpoint version: 1\n";

        const int nlevels = warpx.finestLevel()+1;
        HeaderFile << nlevels << "\n";

        for (int i = 0; i < warpx.getistep().size(); ++i) {
            HeaderFile << warpx.getistep(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.getnsubsteps().size(); ++i) {
            HeaderFile << warpx.getnsubsteps(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.gett_new().size(); ++i) {
            HeaderFile << warpx.gett_new(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.gett_old().size(); ++i) {
            HeaderFile << warpx.gett_old(i) << " ";
        }
        HeaderFile << "\n";

        for (int i = 0; i < warpx.getdt().size(); ++i) {
            HeaderFile << warpx.getdt(i) << " ";
        }
        HeaderFile << "\n";

        HeaderFile << warpx.getmoving_window_x() << "\n";

        HeaderFile << warpx.getis_synchronized() << "\n";

        // Geometry
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << warpx.Geom(0).ProbLo(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << warpx.Geom(0).ProbHi(i) << ' ';
        }
        HeaderFile << '\n';

        // BoxArray
        for (int lev = 0; lev < nlevels; ++lev) {
            warpx.boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }

        warpx.GetPartContainer().WriteHeader(HeaderFile);
    }
}
