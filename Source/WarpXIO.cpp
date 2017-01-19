
#include <MultiFabUtil.H>
#include <PlotFileUtil.H>

#include <WarpX.H>

#include "buildInfo.H"

namespace
{
    const std::string level_prefix {"Level_"};
    constexpr std::streamsize bl_ignore_max { 100000 };
}

void
WarpX::WriteCheckPointFile() const
{
    BL_PROFILE("WarpX::WriteCheckPointFile()");

    const std::string& checkpointname = BoxLib::Concatenate(check_file,istep[0]);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Writing checkpoint " << checkpointname << std::endl;
    }

    const int checkpoint_nfiles = 64;  // could make this parameter
    VisMF::SetNOutFiles(checkpoint_nfiles);
    
    const int nlevels = finestLevel()+1;
    BoxLib::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    if (ParallelDescriptor::IOProcessor())
    {
	std::string HeaderFileName(checkpointname + "/WarpXHeader");
	std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
				                         std::ofstream::trunc |
				                         std::ofstream::binary);
	if( ! HeaderFile.good()) {
	    BoxLib::FileOpenFailed(HeaderFileName);
	}

	HeaderFile.precision(17);

	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
	HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

	HeaderFile << "Checkpoint version: 1\n";

	HeaderFile << nlevels << "\n";

	for (int i = 0; i < istep.size(); ++i) {
	    HeaderFile << istep[i] << " ";
	}
	HeaderFile << "\n";

	for (int i = 0; i < nsubsteps.size(); ++i) {
	    HeaderFile << nsubsteps[i] << " ";
	}
	HeaderFile << "\n";
	
	for (int i = 0; i < t_new.size(); ++i) {
	    HeaderFile << t_new[i] << " ";
	}
	HeaderFile << "\n";

	for (int i = 0; i < t_old.size(); ++i) {
	    HeaderFile << t_old[i] << " ";
	}
	HeaderFile << "\n";

	for (int i = 0; i < dt.size(); ++i) {
	    HeaderFile << dt[i] << " ";
	}
	HeaderFile << "\n";
	
	HeaderFile << moving_window_x << "\n";

	// Geometry
	for (int i = 0; i < BL_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbLo(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbHi(i) << ' ';
	}
        HeaderFile << '\n';

	// BoxArray
	for (int lev = 0; lev < nlevels; ++lev) {
	    boxArray(lev).writeOn(HeaderFile);
	    HeaderFile << '\n';
	}
    }
    
    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlevels; ++lev)
    {
	VisMF::Write(*Efield[lev][0],
		     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ex"));
	VisMF::Write(*Efield[lev][1],
		     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey"));
	VisMF::Write(*Efield[lev][2],
		     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez"));
	VisMF::Write(*Bfield[lev][0],
		     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx"));
	VisMF::Write(*Bfield[lev][1],
		     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By"));
	VisMF::Write(*Bfield[lev][2],
		     BoxLib::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz"));
    }

    mypc->Checkpoint(checkpointname, "particle");
}


void
WarpX::InitFromCheckpoint ()
{
    BL_PROFILE("WarpX::InitFromCheckpoint()");

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Restart from checkpoint " << restart_chkfile << std::endl;
    }

    const int checkpoint_nfiles = 64;  // could make this parameter
    VisMF::SetNOutFiles(checkpoint_nfiles);
    
    // Header
    {
	std::string File(restart_chkfile + "/WarpXHeader");

	VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

	Array<char> fileCharPtr;
	ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
	std::string fileCharPtrString(fileCharPtr.dataPtr());
	std::istringstream is(fileCharPtrString, std::istringstream::in);

	std::string line, word;

	std::getline(is, line);

	int nlevs;
	is >> nlevs;
	is.ignore(bl_ignore_max, '\n');
	finest_level = nlevs-1;

	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		istep[i++] = std::stoi(word);
	    }
	}

	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		nsubsteps[i++] = std::stoi(word);
	    }
	}

	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		t_new[i++] = std::stod(word);
	    }
	}

	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		t_old[i++] = std::stod(word);
	    }
	}

	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		dt[i++] = std::stod(word);
	    }
	}

	is >> moving_window_x;
	is.ignore(bl_ignore_max, '\n');

	Real prob_lo[BL_SPACEDIM];
	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		prob_lo[i++] = std::stod(word);
	    }
	}
	
	Real prob_hi[BL_SPACEDIM];
	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		prob_hi[i++] = std::stod(word);
	    }
	}

	Geometry::ProbDomain(RealBox(prob_lo,prob_hi));

	for (int lev = 0; lev < nlevs; ++lev) {
	    BoxArray ba;
	    ba.readFrom(is);
	    is.ignore(bl_ignore_max, '\n');
	    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
	    MakeNewLevel(lev, ba, dm);
	}
    }

    // Initilize the field data
    for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
    {
	for (int i = 0; i < 3; ++i) {
	    Efield[lev][i]->setVal(0.0);
	    Bfield[lev][i]->setVal(0.0);
	    current[lev][i]->setVal(0.0);
	}

	// xxxxx This will be done differently in amrex!
	{
	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex"));
	    Efield[lev][0]->copy(mf, 0, 0, 1, mf.nGrow(), Efield[lev][0]->nGrow());
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey"));
	    Efield[lev][1]->copy(mf, 0, 0, 1, mf.nGrow(), Efield[lev][1]->nGrow());
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez"));
	    Efield[lev][2]->copy(mf, 0, 0, 1, mf.nGrow(), Efield[lev][2]->nGrow());
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx"));
	    Bfield[lev][0]->copy(mf, 0, 0, 1, mf.nGrow(), Bfield[lev][0]->nGrow());
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By"));
	    Bfield[lev][1]->copy(mf, 0, 0, 1, mf.nGrow(), Bfield[lev][1]->nGrow());
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, BoxLib::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz"));
	    Bfield[lev][2]->copy(mf, 0, 0, 1, mf.nGrow(), Bfield[lev][2]->nGrow());
	}
    }

    // Initilize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile, "particle");
    mypc->Redistribute(false);
}


void
WarpX::WritePlotFile () const
{
    BL_PROFILE("WarpX::WritePlotFile()");

    const std::string& plotfilename = BoxLib::Concatenate(plot_file,istep[0]);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Writing plotfile " << plotfilename << std::endl;
    }
    
    {
	Array<std::string> varnames {"jx", "jy", "jz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"};

	Array<std::unique_ptr<MultiFab> > mf(finest_level+1);
    
	for (int lev = 0; lev <= finest_level; ++lev)
	{
	    const int ncomp = 3*3;
	    const int ngrow = 0;
	    mf[lev].reset(new MultiFab(grids[lev], ncomp, ngrow, dmap[lev]));

	    std::vector<MultiFab*> srcmf(BL_SPACEDIM);
	    PackPlotDataPtrs(srcmf, current[lev]);
	    int dcomp = 0;
	    BoxLib::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
	    WarpX::Copy(*mf[lev], dcomp+1, 1, *current[lev][1], 0);
#endif

	    PackPlotDataPtrs(srcmf, Efield[lev]);
	    dcomp += 3;
	    BoxLib::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
	    WarpX::Copy(*mf[lev], dcomp+1, 1, *Efield[lev][1], 0);
#endif

	    PackPlotDataPtrs(srcmf, Bfield[lev]);
	    dcomp += 3;
	    BoxLib::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
	    WarpX::Copy(*mf[lev], dcomp+1, 1, *Bfield[lev][1], 0);
#endif
	}
    
	Array<const MultiFab*> mf2(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    mf2[lev] = mf[lev].get();
	}

	BoxLib::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf2, varnames,
					Geom(), t_new[0], istep, refRatio());
    }

    mypc->Checkpoint(plotfilename, "particle");

    WriteJobInfo(plotfilename);
}

void
WarpX::WriteJobInfo (const std::string& dir) const
{
    if (ParallelDescriptor::IOProcessor())
    {
	// job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	std::string PrettyLine = "===============================================================================\n";

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

        for (int i = 0; i <= finest_level; i++)
	{
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
	    {
                jobInfoFile << geom[i].Domain().length(n) << " ";
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

void
WarpX::PackPlotDataPtrs(std::vector<MultiFab*>& pmf,
			const Array<std::unique_ptr<MultiFab> >& data)
{
    BL_ASSERT(pmf.size() == BL_SPACEDIM);
#if (BL_SPACEDIM == 3)
    pmf[0] = data[0].get();
    pmf[1] = data[1].get();
    pmf[2] = data[2].get();
#elif (BL_SPACEDIM == 2)
    pmf[0] = data[0].get();
    pmf[1] = data[2].get();
#endif
}
