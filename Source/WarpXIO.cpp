
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <WarpX.H>

#include "AMReX_buildInfo.H"

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

void
WarpX::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void
WarpX::WriteWarpXHeader(const std::string& name) const
{
   if (ParallelDescriptor::IOProcessor())
    {
	std::string HeaderFileName(name + "/WarpXHeader");
	std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
				                         std::ofstream::trunc |
				                         std::ofstream::binary);
	if( ! HeaderFile.good()) {
	    amrex::FileOpenFailed(HeaderFileName);
	}

	HeaderFile.precision(17);

	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
	HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

	HeaderFile << "Checkpoint version: 1\n";

	const int nlevels = finestLevel()+1;
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

	mypc->WriteHeader(HeaderFile);
    }
}

void
WarpX::WriteCheckPointFile() const
{
    BL_PROFILE("WarpX::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0]);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Writing checkpoint " << checkpointname << std::endl;
    }

    const int checkpoint_nfiles = 64;  // could make this parameter
    VisMF::SetNOutFiles(checkpoint_nfiles);
    
    const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteWarpXHeader(checkpointname);
    
    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlevels; ++lev)
    {
	VisMF::Write(*Efield[lev][0],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ex"));
	VisMF::Write(*Efield[lev][1],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey"));
	VisMF::Write(*Efield[lev][2],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez"));
	VisMF::Write(*Bfield[lev][0],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx"));
	VisMF::Write(*Bfield[lev][1],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By"));
	VisMF::Write(*Bfield[lev][2],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz"));
    }

    mypc->Checkpoint(checkpointname, "particle", true);
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
	GotoNextLine(is);
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
	GotoNextLine(is);

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
	    GotoNextLine(is);
	    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
	    MakeNewLevel(lev, ba, dm);
	}

	mypc->ReadHeader(is);
    }

    // Initialize the field data
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
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex"));
	    Efield[lev][0]->copy(mf, 0, 0, 1, 0, 0);
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey"));
	    Efield[lev][1]->copy(mf, 0, 0, 1, 0, 0);
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez"));
	    Efield[lev][2]->copy(mf, 0, 0, 1, 0, 0);
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx"));
	    Bfield[lev][0]->copy(mf, 0, 0, 1, 0, 0);
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By"));
	    Bfield[lev][1]->copy(mf, 0, 0, 1, 0, 0);
	}
	{
	    MultiFab mf;
	    VisMF::Read(mf, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz"));
	    Bfield[lev][2]->copy(mf, 0, 0, 1, 0, 0);
	}
    }

    // Initilize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile, "particle");
}


void
WarpX::WritePlotFile () const
{
    BL_PROFILE("WarpX::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file,istep[0]);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "  Writing plotfile " << plotfilename << std::endl;
    }
    
    {
	Array<std::string> varnames;
	Array<std::unique_ptr<MultiFab> > mf(finest_level+1);
    
	for (int lev = 0; lev <= finest_level; ++lev)
	{
            int ncomp;
            if (ParallelDescriptor::NProcs() > 1) {
	       ncomp = 3*3 + 4; // The +4 is for the number of particles per cell/grid/process + Processor ID
            } else {
	       ncomp = 3*3 + 3; // The +4 is for the number of particles per cell/grid/process + Processor ID
            }

	    const int ngrow = 0;
	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

	    Array<const MultiFab*> srcmf(BL_SPACEDIM);
	    PackPlotDataPtrs(srcmf, current[lev]);
	    int dcomp = 0;
	    amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
	    WarpX::Copy(*mf[lev], dcomp+1, 1, *current[lev][1], 0);
#endif
            if (lev == 0)
            {
                varnames.push_back("jx");
                varnames.push_back("jy");
                varnames.push_back("jz");
            }

	    PackPlotDataPtrs(srcmf, Efield[lev]);
	    dcomp += 3;
	    amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
	    WarpX::Copy(*mf[lev], dcomp+1, 1, *Efield[lev][1], 0);
#endif
            if (lev == 0)
            {
                varnames.push_back("Ex");
                varnames.push_back("Ey");
                varnames.push_back("Ez");
            }

	    PackPlotDataPtrs(srcmf, Bfield[lev]);
	    dcomp += 3;
	    amrex::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
	    WarpX::Copy(*mf[lev], dcomp+1, 1, *Bfield[lev][1], 0);
#endif
            if (lev == 0)
            {
                varnames.push_back("Bx");
                varnames.push_back("By");
                varnames.push_back("Bz");
            }

            MultiFab temp_dat(grids[lev],mf[lev]->DistributionMap(),1,0);
            temp_dat.setVal(0);

            // MultiFab containing number of particles in each cell
	    dcomp += 3;
            mypc->Increment(temp_dat, lev);
            MultiFab::Copy(*mf[lev], temp_dat, 0, dcomp, 1, 0);
            if (lev == 0)
            {
                varnames.push_back("part_per_cell");
            }

            // MultiFab containing number of particles per grid (stored as constant for all cells in each grid)
            const Array<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);

	    dcomp += 1;
            for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
		(*mf[lev])[mfi].setVal(static_cast<Real>(npart_in_grid[mfi.index()]), dcomp);
            }
            if (lev == 0)
            {
                varnames.push_back("part_per_grid");
            }

            // MultiFab containing number of particles per process (stored as constant for all cells in each grid)
	    dcomp += 1;
            long n_per_proc = 0;
            for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
                n_per_proc += npart_in_grid[mfi.index()];
            }
            for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
		(*mf[lev])[mfi].setVal(static_cast<Real>(n_per_proc), dcomp);
            }
            if (lev == 0)
            {
                varnames.push_back("part_per_proc");
            }

            if (ParallelDescriptor::NProcs() > 1) {
                // MultiFab containing the Processor ID 
                dcomp += 1;
                for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
                    (*mf[lev])[mfi].setVal(static_cast<Real>(ParallelDescriptor::MyProc()), dcomp);
                }
                if (lev == 0)
                {
                    varnames.push_back("proc_number");
                }
            }
	}
    
	Array<const MultiFab*> mf2(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    mf2[lev] = mf[lev].get();
	}

	amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf2, varnames,
					Geom(), t_new[0], istep, refRatio());
    }

    mypc->Checkpoint(plotfilename, "particle", false);

    WriteJobInfo(plotfilename);

    WriteWarpXHeader(plotfilename);
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
	jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

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
	  jobInfoFile << "AMReX git hash: " << githash2 << "\n";
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
WarpX::PackPlotDataPtrs(Array<const MultiFab*>& pmf,
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
