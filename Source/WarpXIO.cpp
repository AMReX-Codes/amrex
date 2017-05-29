
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

        HeaderFile << is_synchronized << "\n";

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
// xxxxx
#if 0
    BL_PROFILE("WarpX::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0]);

    amrex::Print() << "  Writing checkpoint " << checkpointname << "\n";

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
        if (is_synchronized) {
            // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
            VisMF::Write(*current[lev][0],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jx"));
            VisMF::Write(*current[lev][1],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jy"));
            VisMF::Write(*current[lev][2],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jz"));
        }

//xxxxx
#if 0
        if (do_pml) {
            const int lev = 0;
            VisMF::Write(*pml_E[0],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml_Ex"));
            VisMF::Write(*pml_E[1],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml_Ey"));
            VisMF::Write(*pml_E[2],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml_Ez"));
            VisMF::Write(*pml_B[0],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml_Bx"));
            VisMF::Write(*pml_B[1],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml_By"));
            VisMF::Write(*pml_B[2],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml_Bz"));
        }
#endif
    }

    mypc->Checkpoint(checkpointname, "particle", true);
#endif
}


void
WarpX::InitFromCheckpoint ()
{
// xxxxx
#if 0

    BL_PROFILE("WarpX::InitFromCheckpoint()");

    amrex::Print() << "  Restart from checkpoint " << restart_chkfile << "\n";

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

        is >> is_synchronized;
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
            SetBoxArray(lev, ba);
            SetDistributionMap(lev, dm);
	    AllocLevelData(lev, ba, dm);
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

        VisMF::Read(*Efield[lev][0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex"));
        VisMF::Read(*Efield[lev][1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey"));
        VisMF::Read(*Efield[lev][2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez"));

        VisMF::Read(*Bfield[lev][0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx"));
        VisMF::Read(*Bfield[lev][1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By"));
        VisMF::Read(*Bfield[lev][2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz"));

        if (is_synchronized) {
            VisMF::Read(*current[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jx"));
            VisMF::Read(*current[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jy"));
            VisMF::Read(*current[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jz"));    
        }
    }

// xxxxx
#if 0
    if (do_pml)
    {
        InitPML();

        const int lev = 0;
        VisMF::Read(*pml_E[0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_Ex"));
        VisMF::Read(*pml_E[1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_Ey"));
        VisMF::Read(*pml_E[2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_Ez"));

        VisMF::Read(*pml_B[0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_Bx"));
        VisMF::Read(*pml_B[1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_By"));
        VisMF::Read(*pml_B[2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_Bz"));
    }
#endif

    // Initilize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile, "particle");
#endif
}


void
WarpX::WritePlotFile () const
{
    BL_PROFILE("WarpX::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file,istep[0]);

    amrex::Print() << "  Writing plotfile " << plotfilename << "\n";
    
    {
	Array<std::string> varnames;
	Array<std::unique_ptr<MultiFab> > mf(finest_level+1);

        const int ncomp = 3*3 
            + static_cast<int>(plot_part_per_cell)
            + static_cast<int>(plot_part_per_grid)
            + static_cast<int>(plot_part_per_proc)
            + static_cast<int>(plot_proc_number)
            + static_cast<int>(plot_divb)
            + static_cast<int>(plot_finepatch)*6;

	for (int lev = 0; lev <= finest_level; ++lev)
	{
	    const int ngrow = 0;
	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

	    Array<const MultiFab*> srcmf(BL_SPACEDIM);
	    PackPlotDataPtrs(srcmf, current_fp[lev]);
	    int dcomp = 0;
	    amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter(*mf[lev], dcomp+1, *current[lev][1], 0, 1);
#endif
            if (lev == 0)
            {
                varnames.push_back("jx");
                varnames.push_back("jy");
                varnames.push_back("jz");
            }
	    dcomp += 3;

	    PackPlotDataPtrs(srcmf, Efield_aux[lev]);
	    amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter(*mf[lev], dcomp+1, *Efield_aux[lev][1], 0, 1);
#endif
            if (lev == 0)
            {
                varnames.push_back("Ex");
                varnames.push_back("Ey");
                varnames.push_back("Ez");
            }
	    dcomp += 3;

	    PackPlotDataPtrs(srcmf, Bfield_aux[lev]);
	    amrex::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
            MultiFab::Copy(*mf[lev], *Bfield_aux[lev][1], 0, dcomp+1, 1, ngrow);
#endif
            if (lev == 0)
            {
                varnames.push_back("Bx");
                varnames.push_back("By");
                varnames.push_back("Bz");
            }
	    dcomp += 3;

            if (plot_part_per_cell)
            {
                MultiFab temp_dat(grids[lev],mf[lev]->DistributionMap(),1,0);
                temp_dat.setVal(0);
                
                // MultiFab containing number of particles in each cell
                mypc->Increment(temp_dat, lev);
                MultiFab::Copy(*mf[lev], temp_dat, 0, dcomp, 1, 0);
                if (lev == 0)
                {
                    varnames.push_back("part_per_cell");
                }
                dcomp += 1;
            }

            if (plot_part_per_grid || plot_part_per_proc)
            {
                const Array<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);

                if (plot_part_per_grid)
                {
                    // MultiFab containing number of particles per grid (stored as constant for all cells in each grid)
                    for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
                        (*mf[lev])[mfi].setVal(static_cast<Real>(npart_in_grid[mfi.index()]), dcomp);
                    }
                    if (lev == 0)
                    {
                        varnames.push_back("part_per_grid");
                    }
                    dcomp += 1;
                }

                if (plot_part_per_proc)
                {
                    // MultiFab containing number of particles per process (stored as constant for all cells in each grid)
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
                    dcomp += 1;
                }
            }

            if (plot_proc_number)
            {
                // MultiFab containing the Processor ID 
                for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
                    (*mf[lev])[mfi].setVal(static_cast<Real>(ParallelDescriptor::MyProc()), dcomp);
                }
                if (lev == 0)
                {
                    varnames.push_back("proc_number");
                }
                dcomp += 1;
            }

            if (plot_divb)
            {
                PackPlotDataPtrs(srcmf, Bfield_aux[lev]);
                ComputeDivB(*mf[lev], dcomp, srcmf, Geom(lev).CellSize());
                if (lev == 0)
                {
                    varnames.push_back("divB");
                }                
                dcomp += 1;
            }

            if (plot_finepatch)
            {
                PackPlotDataPtrs(srcmf, Efield_fp[lev]);
                amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
                MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
                amrex::average_node_to_cellcenter(*mf[lev], dcomp+1, *Efield_fp[lev][1], 0, 1);
#endif
                if (lev == 0)
                {
                    varnames.push_back("Ex_fp");
                    varnames.push_back("Ey_fp");
                    varnames.push_back("Ez_fp");
                }
                dcomp += 3;
                
                PackPlotDataPtrs(srcmf, Bfield_fp[lev]);
                amrex::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (BL_SPACEDIM == 2)
                MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
                MultiFab::Copy(*mf[lev], *Bfield_fp[lev][1], 0, dcomp+1, 1, ngrow);
#endif
                if (lev == 0)
                {
                    varnames.push_back("Bx_fp");
                    varnames.push_back("By_fp");
                    varnames.push_back("Bz_fp");
                }
                dcomp += 3;
            }

            BL_ASSERT(dcomp == ncomp);
	}

#if 0
        for (int lev = finest_level; lev > 0; --lev)
        {
            amrex::average_down(*mf[lev], *mf[lev-1], 0, ncomp, refRatio(lev-1));
        }
#endif
    
	amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, 
                                       amrex::GetArrOfConstPtrs(mf),
                                       varnames, Geom(), t_new[0], istep, refRatio());
    }

    if (plot_raw_fields || plot_crsepatch)
    {
        const int raw_plot_nfiles = 64;  // could make this parameter
        VisMF::SetNOutFiles(raw_plot_nfiles);

        const int nlevels = finestLevel()+1;
        const std::string raw_plotfilename = plotfilename + "/raw_fields";
        amrex::PreBuildDirectorHierarchy(raw_plotfilename, level_prefix, nlevels, true);

        for (int lev = 0; lev < nlevels; ++lev)
        {
            const DistributionMapping& dm = DistributionMap(lev);

            MultiFab Ex(Efield_aux[lev][0]->boxArray(), dm, 1, 0);
            MultiFab Ey(Efield_aux[lev][1]->boxArray(), dm, 1, 0);
            MultiFab Ez(Efield_aux[lev][2]->boxArray(), dm, 1, 0);
            MultiFab Bx(Bfield_aux[lev][0]->boxArray(), dm, 1, 0);
            MultiFab By(Bfield_aux[lev][1]->boxArray(), dm, 1, 0);
            MultiFab Bz(Bfield_aux[lev][2]->boxArray(), dm, 1, 0);
            MultiFab jx(current_fp[lev][0]->boxArray(), dm, 1, 0);
            MultiFab jy(current_fp[lev][1]->boxArray(), dm, 1, 0);
            MultiFab jz(current_fp[lev][2]->boxArray(), dm, 1, 0);

            MultiFab::Copy(Ex, *Efield_aux[lev][0], 0, 0, 1, 0);
            MultiFab::Copy(Ey, *Efield_aux[lev][1], 0, 0, 1, 0);
            MultiFab::Copy(Ez, *Efield_aux[lev][2], 0, 0, 1, 0);
            MultiFab::Copy(Bx, *Bfield_aux[lev][0], 0, 0, 1, 0);
            MultiFab::Copy(By, *Bfield_aux[lev][1], 0, 0, 1, 0);
            MultiFab::Copy(Bz, *Bfield_aux[lev][2], 0, 0, 1, 0);
            MultiFab::Copy(jx, *current_fp[lev][0], 0, 0, 1, 0);
            MultiFab::Copy(jy, *current_fp[lev][1], 0, 0, 1, 0);
            MultiFab::Copy(jz, *current_fp[lev][2], 0, 0, 1, 0);

            VisMF::Write(Ex, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex"));
            VisMF::Write(Ey, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey"));
            VisMF::Write(Ez, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez"));
            VisMF::Write(Bx, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx"));
            VisMF::Write(By, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By"));
            VisMF::Write(Bz, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz"));
            VisMF::Write(jx, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jx"));
            VisMF::Write(jy, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jy"));
            VisMF::Write(jz, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jz"));

            if (plot_crsepatch && Bfield_cp[lev][0])
            {
                VisMF::Write(*Efield_cp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_cp"));
                VisMF::Write(*Efield_cp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_cp"));
                VisMF::Write(*Efield_cp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_cp"));
                VisMF::Write(*Bfield_cp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_cp"));
                VisMF::Write(*Bfield_cp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_cp"));
                VisMF::Write(*Bfield_cp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_cp"));
            }
        }
    }

    Array<std::string> particle_varnames;
    particle_varnames.push_back("weight");

    particle_varnames.push_back("velocity_x");
    particle_varnames.push_back("velocity_y");
    particle_varnames.push_back("velocity_z");

    particle_varnames.push_back("Ex");
    particle_varnames.push_back("Ey");
    particle_varnames.push_back("Ez");

    particle_varnames.push_back("Bx");
    particle_varnames.push_back("By");
    particle_varnames.push_back("Bz");

    mypc->Checkpoint(plotfilename, "particle", true, particle_varnames);

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
WarpX::PackPlotDataPtrs (Array<const MultiFab*>& pmf,
                         const std::array<std::unique_ptr<MultiFab>,3>& data)
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
