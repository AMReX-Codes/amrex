
#ifdef _OPENMP
#include <omp.h>
#endif

#include <MultiFabUtil.H>
#include <PlotFileUtil.H>

#include <WarpX.H>

#include "buildInfo.H"

void
WarpX::WritePlotFile (int istep, Real t) const
{
    BL_PROFILE("WarpX::WritePlotFile()");

    const int finest_lev = 0;

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "    Writing plotfile " << istep << std::endl;
    }
    
    int ncomp = 3*3;
    int ngrow = 0;
    int lev = 0;
    MultiFab mf(ba_arr[lev], ncomp, ngrow, dmap_arr[lev]);

    std::vector<std::string> varnames {"jx", "jy", "jz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"};

    std::vector<MultiFab*> srcmf(3);
    for (int i = 0; i < 3; ++i) {
	srcmf[i] = current[i].get();
    }
    int dcomp = 0;
    BoxLib::average_edge_to_cellcenter(mf, dcomp, srcmf);

    for (int i = 0; i < 3; ++i) {
	srcmf[i] = Efield[i].get();
    }
    dcomp += 3;
    BoxLib::average_edge_to_cellcenter(mf, dcomp, srcmf);

    for (int i = 0; i < 3; ++i) {
	srcmf[i] = Bfield[i].get();
    }
    dcomp += 3;
    BoxLib::average_face_to_cellcenter(mf, dcomp, srcmf);
    
    const std::string& plotfilename = BoxLib::Concatenate("plt",istep);
    BoxLib::WriteSingleLevelPlotfile(plotfilename, mf, varnames, geom_arr[lev], t);

    mypc->Checkpoint(plotfilename, "particle");

    if (false && ParallelDescriptor::IOProcessor()) {
	// job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = plotfilename;
	std::string PrettyLine = "===============================================================================\n";

	FullPathJobInfoFile += "/job_info";
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
	jobInfoFile << "BoxLib dir:    " << buildInfoGetBoxlibDir() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
	jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
	jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
	jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	const char* githash3 = buildInfoGetGitHash(3);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "WarpX  git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "BoxLib git hash: " << githash2 << "\n";
	}
	if (strlen(githash3) > 0) {
	  jobInfoFile << "PICSAR git hash: " << githash3 << "\n";
	}

	jobInfoFile << "\n\n";

	// grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= finest_lev; i++)
	{
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << ba_arr[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
	    {
                jobInfoFile << geom_arr[i].Domain().length(n) << " ";
	    }
            jobInfoFile << "\n\n";
	}

        jobInfoFile << " Boundary conditions\n";

        jobInfoFile << "   -x: " << "interior" << "\n";
        jobInfoFile << "   +x: " << "interior" << "\n";
        if (BL_SPACEDIM >= 2) {
	    jobInfoFile << "   -y: " << "interior" << "\n";
	    jobInfoFile << "   +y: " << "interior" << "\n";
        }
        if (BL_SPACEDIM == 3) {
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
