
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil_F.H>

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
	for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbLo(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
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

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(checkpoint_headerversion);

    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0]);

    amrex::Print() << "  Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteWarpXHeader(checkpointname);

    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlevels; ++lev)
    {
	VisMF::Write(*Efield_fp[lev][0],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ex_fp"));
	VisMF::Write(*Efield_fp[lev][1],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey_fp"));
	VisMF::Write(*Efield_fp[lev][2],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez_fp"));
	VisMF::Write(*Bfield_fp[lev][0],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx_fp"));
	VisMF::Write(*Bfield_fp[lev][1],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By_fp"));
	VisMF::Write(*Bfield_fp[lev][2],
		     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz_fp"));
        if (is_synchronized) {
            // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
            VisMF::Write(*current_fp[lev][0],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jx_fp"));
            VisMF::Write(*current_fp[lev][1],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jy_fp"));
            VisMF::Write(*current_fp[lev][2],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jz_fp"));
        }

        if (lev > 0)
        {
            VisMF::Write(*Efield_cp[lev][0],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ex_cp"));
            VisMF::Write(*Efield_cp[lev][1],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey_cp"));
            VisMF::Write(*Efield_cp[lev][2],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez_cp"));
            VisMF::Write(*Bfield_cp[lev][0],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx_cp"));
            VisMF::Write(*Bfield_cp[lev][1],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By_cp"));
            VisMF::Write(*Bfield_cp[lev][2],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz_cp"));
            if (is_synchronized) {
                // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
                VisMF::Write(*current_cp[lev][0],
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jx_cp"));
                VisMF::Write(*current_cp[lev][1],
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jy_cp"));
                VisMF::Write(*current_cp[lev][2],
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jz_cp"));
            }
        }

        if (do_pml && pml[lev]) {
            pml[lev]->CheckPoint(amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml"));
        }

        if (costs[lev]) {
            VisMF::Write(*costs[lev],
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "costs"));
        }
    }

    mypc->Checkpoint(checkpointname, "particle", true);

    VisMF::SetHeaderVersion(current_version);
}


void
WarpX::InitFromCheckpoint ()
{
    BL_PROFILE("WarpX::InitFromCheckpoint()");

    amrex::Print() << "  Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    {
	std::string File(restart_chkfile + "/WarpXHeader");

	VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

	Vector<char> fileCharPtr;
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

	Real prob_lo[AMREX_SPACEDIM];
	std::getline(is, line);
	{
	    std::istringstream lis(line);
	    int i = 0;
	    while (lis >> word) {
		prob_lo[i++] = std::stod(word);
	    }
	}

	Real prob_hi[AMREX_SPACEDIM];
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

    const int nlevs = finestLevel()+1;

    // Initialize the field data
    for (int lev = 0; lev < nlevs; ++lev)
    {
        for (int i = 0; i < 3; ++i) {
            current_fp[lev][i]->setVal(0.0);
            Efield_fp[lev][i]->setVal(0.0);
            Bfield_fp[lev][i]->setVal(0.0);
        }

        if (lev > 0) {
            for (int i = 0; i < 3; ++i) {
                Efield_aux[lev][i]->setVal(0.0);
                Bfield_aux[lev][i]->setVal(0.0);

                current_cp[lev][i]->setVal(0.0);
                Efield_cp[lev][i]->setVal(0.0);
                Bfield_cp[lev][i]->setVal(0.0);
            }
        }

        VisMF::Read(*Efield_fp[lev][0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex_fp"));
        VisMF::Read(*Efield_fp[lev][1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey_fp"));
        VisMF::Read(*Efield_fp[lev][2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez_fp"));

        VisMF::Read(*Bfield_fp[lev][0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx_fp"));
        VisMF::Read(*Bfield_fp[lev][1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By_fp"));
        VisMF::Read(*Bfield_fp[lev][2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz_fp"));

        if (is_synchronized) {
            VisMF::Read(*current_fp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jx_fp"));
            VisMF::Read(*current_fp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jy_fp"));
            VisMF::Read(*current_fp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jz_fp"));
        }

        if (lev > 0)
        {
            VisMF::Read(*Efield_cp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex_cp"));
            VisMF::Read(*Efield_cp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey_cp"));
            VisMF::Read(*Efield_cp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez_cp"));

            VisMF::Read(*Bfield_cp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx_cp"));
            VisMF::Read(*Bfield_cp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By_cp"));
            VisMF::Read(*Bfield_cp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz_cp"));

            if (is_synchronized) {
                VisMF::Read(*current_cp[lev][0],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jx_cp"));
                VisMF::Read(*current_cp[lev][1],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jy_cp"));
                VisMF::Read(*current_cp[lev][2],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jz_cp"));
            }
        }

        if (costs[lev]) {
            const auto& cost_mf_name = 
                amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "costs");
            if (VisMF::Exist(cost_mf_name)) {
                VisMF::Read(*costs[lev], cost_mf_name);
            } else {
                costs[lev]->setVal(0.0);
            }
        }
    }

    if (do_pml)
    {
        InitPML();
        for (int lev = 0; lev < nlevs; ++lev) {
            pml[lev]->Restart(amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml"));
        }
    }

    // Initilize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile, "particle");

#ifdef WARPX_DO_ELECTROSTATIC
    if (do_electrostatic) {
        getLevelMasks(masks);

        // the plus one is to convert from num_cells to num_nodes
        getLevelMasks(gather_masks, 4 + 1);
    }
#endif // WARPX_DO_ELECTROSTATIC
}


std::unique_ptr<MultiFab>
WarpX::GetCellCenteredData() {

    BL_PROFILE("WarpX::GetCellCenteredData");

    const int ng =  1;
    const int nc = 10;
    const int lev = 0;
    auto cc = std::unique_ptr<MultiFab>( new MultiFab(boxArray(lev),
                                                      DistributionMap(lev),
                                                      nc, ng) );

    Vector<const MultiFab*> srcmf(AMREX_SPACEDIM);
    int dcomp = 0;

    // first the electric field
    PackPlotDataPtrs(srcmf, Efield_aux[lev]);
    amrex::average_edge_to_cellcenter(*cc, dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
    MultiFab::Copy(*cc, *cc, dcomp+1, dcomp+2, 1, ng);
    amrex::average_node_to_cellcenter(*cc, dcomp+1, *Efield_aux[lev][1], 0, 1);
#endif
    dcomp += 3;
    
    // then the magnetic field
    PackPlotDataPtrs(srcmf, Bfield_aux[lev]);
    amrex::average_face_to_cellcenter(*cc, dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
    MultiFab::Copy(*cc, *cc, dcomp+1, dcomp+2, 1, ng);
    MultiFab::Copy(*cc, *Bfield_aux[lev][1], 0, dcomp+1, 1, ng);
#endif
    dcomp += 3;

    // then the current density
    PackPlotDataPtrs(srcmf, current_fp[lev]);
    amrex::average_edge_to_cellcenter(*cc, dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
    MultiFab::Copy(*cc, *cc, dcomp+1, dcomp+2, 1, ng);
    amrex::average_node_to_cellcenter(*cc, dcomp+1, *current_fp[lev][1], 0, 1);
#endif
    dcomp += 3;

    const std::unique_ptr<MultiFab>& charge_density = mypc->GetChargeDensity(lev);
    amrex::average_node_to_cellcenter(*cc, dcomp, *charge_density, 0, 1);

    cc->FillBoundary(geom[lev].periodicity());
    
    return cc;
}

void
WarpX::WritePlotFile () const
{
    BL_PROFILE("WarpX::WritePlotFile()");

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(plotfile_headerversion);

    const std::string& plotfilename = amrex::Concatenate(plot_file,istep[0]);

    amrex::Print() << "  Writing plotfile " << plotfilename << "\n";

    {
	Vector<std::string> varnames;
	Vector<std::unique_ptr<MultiFab> > mf(finest_level+1);

        const int ncomp = 3*3
            + static_cast<int>(plot_part_per_cell)
            + static_cast<int>(plot_part_per_grid)
            + static_cast<int>(plot_part_per_proc)
            + static_cast<int>(plot_proc_number)
            + static_cast<int>(plot_divb)
            + static_cast<int>(plot_dive)
            + static_cast<int>(plot_finepatch)*6
            + static_cast<int>(plot_crsepatch)*6
            + static_cast<int>(costs[0] != nullptr);

	for (int lev = 0; lev <= finest_level; ++lev)
	{
	    const int ngrow = 0;
	    mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

	    Vector<const MultiFab*> srcmf(AMREX_SPACEDIM);
	    PackPlotDataPtrs(srcmf, current_fp[lev]);
	    int dcomp = 0;
	    amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
	    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter(*mf[lev], dcomp+1, *current_fp[lev][1], 0, 1);
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
#if (AMREX_SPACEDIM == 2)
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
#if (AMREX_SPACEDIM == 2)
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
                const Vector<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);

                if (plot_part_per_grid)
                {
                    // MultiFab containing number of particles per grid (stored as constant for all cells in each grid)
#ifdef _OPENMP
#pragma omp parallel
#endif
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
#ifdef _OPENMP
#pragma omp parallel reduction(+:n_per_proc)
#endif
                    for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
                        n_per_proc += npart_in_grid[mfi.index()];
                    }
                    mf[lev]->setVal(static_cast<Real>(n_per_proc), dcomp,1);
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
#ifdef _OPENMP
#pragma omp parallel
#endif
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
                ComputeDivB(*mf[lev], dcomp,
                            {Bfield_aux[lev][0].get(),Bfield_aux[lev][1].get(),Bfield_aux[lev][2].get()},
                            WarpX::CellSize(lev));
                if (lev == 0)
                {
                    varnames.push_back("divB");
                }
                dcomp += 1;
            }

            if (plot_dive)
            {
                const BoxArray& ba = amrex::convert(boxArray(lev),IntVect::TheUnitVector());
                MultiFab dive(ba,DistributionMap(lev),1,0);
                ComputeDivE(dive, 0,
                            {Efield_aux[lev][0].get(), Efield_aux[lev][1].get(), Efield_aux[lev][2].get()},
                            WarpX::CellSize(lev));
                amrex::average_node_to_cellcenter(*mf[lev], dcomp, dive, 0, 1);
                if (lev == 0)
                {
                    varnames.push_back("divE");
                }
                dcomp += 1;
            }

            if (plot_finepatch)
            {
                PackPlotDataPtrs(srcmf, Efield_fp[lev]);
                amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
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
#if (AMREX_SPACEDIM == 2)
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

            if (plot_crsepatch)
            {
                // First the electric field
                if (lev == 0)
                {
                    mf[lev]->setVal(0.0, dcomp, AMREX_SPACEDIM, ngrow);
                }
                else
                {
                    std::array<std::unique_ptr<MultiFab>, 3> E = getInterpolatedE(lev);
                    PackPlotDataPtrs(srcmf, E);
                    amrex::average_edge_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
                    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
                    amrex::average_node_to_cellcenter(*mf[lev], dcomp+1, *E[1], 0, 1);
#endif
                }
                if (lev == 0)
                {
                    varnames.push_back("Ex_cp");
                    varnames.push_back("Ey_cp");
                    varnames.push_back("Ez_cp");
                }
                dcomp += 3;

                // now the magnetic field                
                if (lev == 0)
                {
                    mf[lev]->setVal(0.0, dcomp, AMREX_SPACEDIM, ngrow);
                }
                else
                {
                    std::array<std::unique_ptr<MultiFab>, 3> B = getInterpolatedB(lev);
                    PackPlotDataPtrs(srcmf, B);
                    amrex::average_face_to_cellcenter(*mf[lev], dcomp, srcmf);
#if (AMREX_SPACEDIM == 2)
                    MultiFab::Copy(*mf[lev], *mf[lev], dcomp+1, dcomp+2, 1, ngrow);
                    MultiFab::Copy(*mf[lev], *B[1], 0, dcomp+1, 1, ngrow);
#endif
                }
                if (lev == 0)
                {
                    varnames.push_back("Bx_cp");
                    varnames.push_back("By_cp");
                    varnames.push_back("Bz_cp");
                }
                dcomp += 3;
            }
            
            if (costs[0] != nullptr)
            {
                MultiFab::Copy(*mf[lev], *costs[lev], 0, dcomp, 1, 0);
                if (lev == 0)
                {
                    varnames.push_back("costs");
                }
                dcomp += 1;
            }

            BL_ASSERT(dcomp == ncomp);
	}

#if 0
        for (int lev = finest_level; lev > 0; --lev)
        {
            amrex::average_down(*mf[lev], *mf[lev-1], 0, ncomp, refRatio(lev-1));
        }
#endif

        Vector<std::string> rfs;
        if (plot_raw_fields) rfs.emplace_back("raw_fields"); // pre-build raw_fields/
	amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                       amrex::GetVecOfConstPtrs(mf),
                                       varnames, Geom(), t_new[0], istep, refRatio(),
                                       "HyperCLaw-V1.1",
                                       "Level_",
                                       "Cell",
                                       rfs);
    }

    if (plot_raw_fields)
    {
        const int nlevels = finestLevel()+1;
        for (int lev = 0; lev < nlevels; ++lev)
        {
            const std::string raw_plotfilename = plotfilename + "/raw_fields";
            // Plot auxilary patch
            if (plot_raw_fields_guards) {
                VisMF::Write(*Efield_aux[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_aux"));
                VisMF::Write(*Efield_aux[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_aux"));
                VisMF::Write(*Efield_aux[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_aux"));
                VisMF::Write(*Bfield_aux[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_aux"));
                VisMF::Write(*Bfield_aux[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_aux"));
                VisMF::Write(*Bfield_aux[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_aux"));
            } else {
                const DistributionMapping& dm = DistributionMap(lev);
                MultiFab Ex(Efield_aux[lev][0]->boxArray(), dm, 1, 0);
                MultiFab Ey(Efield_aux[lev][1]->boxArray(), dm, 1, 0);
                MultiFab Ez(Efield_aux[lev][2]->boxArray(), dm, 1, 0);
                MultiFab Bx(Bfield_aux[lev][0]->boxArray(), dm, 1, 0);
                MultiFab By(Bfield_aux[lev][1]->boxArray(), dm, 1, 0);
                MultiFab Bz(Bfield_aux[lev][2]->boxArray(), dm, 1, 0);
                MultiFab::Copy(Ex, *Efield_aux[lev][0], 0, 0, 1, 0);
                MultiFab::Copy(Ey, *Efield_aux[lev][1], 0, 0, 1, 0);
                MultiFab::Copy(Ez, *Efield_aux[lev][2], 0, 0, 1, 0);
                MultiFab::Copy(Bx, *Bfield_aux[lev][0], 0, 0, 1, 0);
                MultiFab::Copy(By, *Bfield_aux[lev][1], 0, 0, 1, 0);
                MultiFab::Copy(Bz, *Bfield_aux[lev][2], 0, 0, 1, 0);
                VisMF::Write(Ex, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_aux"));
                VisMF::Write(Ey, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_aux"));
                VisMF::Write(Ez, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_aux"));
                VisMF::Write(Bx, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_aux"));
                VisMF::Write(By, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_aux"));
                VisMF::Write(Bz, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_aux"));
            }

            // Plot fine patch
            if (plot_finepatch) {
            if (plot_raw_fields_guards) {
                VisMF::Write(*Efield_fp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_fp"));
                VisMF::Write(*Efield_fp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_fp"));
                VisMF::Write(*Efield_fp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_fp"));
                VisMF::Write(*Bfield_fp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_fp"));
                VisMF::Write(*Bfield_fp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_fp"));
                VisMF::Write(*Bfield_fp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_fp"));
                VisMF::Write(*current_fp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jx_fp"));
                VisMF::Write(*current_fp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jy_fp"));
                VisMF::Write(*current_fp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jz_fp"));
            } else {
                const DistributionMapping& dm = DistributionMap(lev);
                MultiFab Ex(Efield_fp[lev][0]->boxArray(), dm, 1, 0);
                MultiFab Ey(Efield_fp[lev][1]->boxArray(), dm, 1, 0);
                MultiFab Ez(Efield_fp[lev][2]->boxArray(), dm, 1, 0);
                MultiFab Bx(Bfield_fp[lev][0]->boxArray(), dm, 1, 0);
                MultiFab By(Bfield_fp[lev][1]->boxArray(), dm, 1, 0);
                MultiFab Bz(Bfield_fp[lev][2]->boxArray(), dm, 1, 0);
                MultiFab jx(current_fp[lev][0]->boxArray(), dm, 1, 0);
                MultiFab jy(current_fp[lev][1]->boxArray(), dm, 1, 0);
                MultiFab jz(current_fp[lev][2]->boxArray(), dm, 1, 0);
                MultiFab::Copy(Ex, *Efield_fp[lev][0], 0, 0, 1, 0);
                MultiFab::Copy(Ey, *Efield_fp[lev][1], 0, 0, 1, 0);
                MultiFab::Copy(Ez, *Efield_fp[lev][2], 0, 0, 1, 0);
                MultiFab::Copy(Bx, *Bfield_fp[lev][0], 0, 0, 1, 0);
                MultiFab::Copy(By, *Bfield_fp[lev][1], 0, 0, 1, 0);
                MultiFab::Copy(Bz, *Bfield_fp[lev][2], 0, 0, 1, 0);
                MultiFab::Copy(jx, *current_fp[lev][0], 0, 0, 1, 0);
                MultiFab::Copy(jy, *current_fp[lev][1], 0, 0, 1, 0);
                MultiFab::Copy(jz, *current_fp[lev][2], 0, 0, 1, 0);
                VisMF::Write(Ex, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_fp"));
                VisMF::Write(Ey, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_fp"));
                VisMF::Write(Ez, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_fp"));
                VisMF::Write(Bx, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_fp"));
                VisMF::Write(By, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_fp"));
                VisMF::Write(Bz, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_fp"));
                VisMF::Write(jx, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jx_fp"));
                VisMF::Write(jy, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jy_fp"));
                VisMF::Write(jz, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "jz_fp"));
            }
            }

            // Plot coarse patch
            if (plot_crsepatch)
            {
                if (lev == 0)
                {
                    const DistributionMapping& dm = DistributionMap(lev);
                    MultiFab Ex(Efield_aux[lev][0]->boxArray(), dm, 1, 0);
                    MultiFab Ey(Efield_aux[lev][1]->boxArray(), dm, 1, 0);
                    MultiFab Ez(Efield_aux[lev][2]->boxArray(), dm, 1, 0);
                    MultiFab Bx(Bfield_aux[lev][0]->boxArray(), dm, 1, 0);
                    MultiFab By(Bfield_aux[lev][1]->boxArray(), dm, 1, 0);
                    MultiFab Bz(Bfield_aux[lev][2]->boxArray(), dm, 1, 0);

                    Ex.setVal(0.0); Ey.setVal(0.0); Ez.setVal(0.0);
                    Bx.setVal(0.0); By.setVal(0.0); Bz.setVal(0.0);

                    VisMF::Write(Ex, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_cp"));
                    VisMF::Write(Ey, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_cp"));
                    VisMF::Write(Ez, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_cp"));
                    VisMF::Write(Bx, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_cp"));
                    VisMF::Write(By, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_cp"));
                    VisMF::Write(Bz, amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_cp"));              
                } else {

                    if (plot_raw_fields_guards) {
                        VisMF::Write(*Efield_cp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_cp"));
                        VisMF::Write(*Efield_cp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_cp"));
                        VisMF::Write(*Efield_cp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_cp"));
                        VisMF::Write(*Bfield_cp[lev][0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_cp"));
                        VisMF::Write(*Bfield_cp[lev][1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_cp"));
                        VisMF::Write(*Bfield_cp[lev][2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_cp"));
                    } else {
                        std::array<std::unique_ptr<MultiFab>, 3> E = getInterpolatedE(lev);
                        std::array<std::unique_ptr<MultiFab>, 3> B = getInterpolatedB(lev);

                        VisMF::Write(*E[0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ex_cp"));
                        VisMF::Write(*E[1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ey_cp"));
                        VisMF::Write(*E[2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Ez_cp"));
                        VisMF::Write(*B[0], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bx_cp"));
                        VisMF::Write(*B[1], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "By_cp"));
                        VisMF::Write(*B[2], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "Bz_cp"));
                    }
                }
            }
            
            if (F_fp[lev]) {
                VisMF::Write(*F_fp[lev], amrex::MultiFabFileFullPrefix(lev, raw_plotfilename, level_prefix, "F_fp"));
            }

        }
    }

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("weight");

    particle_varnames.push_back("momentum_x");
    particle_varnames.push_back("momentum_y");
    particle_varnames.push_back("momentum_z");

    particle_varnames.push_back("Ex");
    particle_varnames.push_back("Ey");
    particle_varnames.push_back("Ez");

    particle_varnames.push_back("Bx");
    particle_varnames.push_back("By");
    particle_varnames.push_back("Bz");

#ifdef WARPX_STORE_OLD_PARTICLE_ATTRIBS
    particle_varnames.push_back("xold");
    particle_varnames.push_back("yold");
    particle_varnames.push_back("zold");

    particle_varnames.push_back("uxold");
    particle_varnames.push_back("uyold");
    particle_varnames.push_back("uzold");    
#endif
    
    mypc->Checkpoint(plotfilename, "particle", true, particle_varnames);

    WriteJobInfo(plotfilename);

    WriteWarpXHeader(plotfilename);

    VisMF::SetHeaderVersion(current_version);
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
            for (int n = 0; n < AMREX_SPACEDIM; n++)
	    {
                jobInfoFile << geom[i].Domain().length(n) << " ";
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

std::array<std::unique_ptr<MultiFab>, 3> WarpX::getInterpolatedE(int lev) const
{

    const int ngrow = 0;
    
    std::array<std::unique_ptr<MultiFab>, 3> interpolated_E;
    for (int i = 0; i < 3; ++i) {
        interpolated_E[i].reset( new MultiFab(Efield_aux[lev][i]->boxArray(), dmap[lev], 1, ngrow) );
        interpolated_E[i]->setVal(0.0);
    }
    
    const int r_ratio = refRatio(lev-1)[0];
    const int use_limiter = 0;
#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,3> efab;
        for (MFIter mfi(*interpolated_E[0]); mfi.isValid(); ++mfi)
        {
            Box ccbx = mfi.fabbox();
            ccbx.enclosedCells();
            ccbx.coarsen(r_ratio).refine(r_ratio); // so that ccbx is coarsenable
            
            const FArrayBox& cxfab = (*Efield_cp[lev][0])[mfi];
            const FArrayBox& cyfab = (*Efield_cp[lev][1])[mfi];
            const FArrayBox& czfab = (*Efield_cp[lev][2])[mfi];
            
            efab[0].resize(amrex::convert(ccbx,Ex_nodal_flag));
            efab[1].resize(amrex::convert(ccbx,Ey_nodal_flag));
            efab[2].resize(amrex::convert(ccbx,Ez_nodal_flag));
            
#if (AMREX_SPACEDIM == 3)
            amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                BL_TO_FORTRAN_ANYD(efab[0]),
                                BL_TO_FORTRAN_ANYD(efab[1]),
                                BL_TO_FORTRAN_ANYD(efab[2]),
                                BL_TO_FORTRAN_ANYD(cxfab),
                                BL_TO_FORTRAN_ANYD(cyfab),
                                BL_TO_FORTRAN_ANYD(czfab),
                                &r_ratio,&use_limiter);
#else
            amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                BL_TO_FORTRAN_ANYD(efab[0]),
                                BL_TO_FORTRAN_ANYD(efab[2]),
                                BL_TO_FORTRAN_ANYD(cxfab),
                                BL_TO_FORTRAN_ANYD(czfab),
                                &r_ratio,&use_limiter);
            amrex_interp_nd_efield(ccbx.loVect(), ccbx.hiVect(),
                                   BL_TO_FORTRAN_ANYD(efab[1]),
                                   BL_TO_FORTRAN_ANYD(cyfab),
                                   &r_ratio);
#endif

            for (int i = 0; i < 3; ++i) {
                const Box& bx = (*interpolated_E[i])[mfi].box();
                (*interpolated_E[i])[mfi].plus(efab[i], bx, bx, 0, 0, 1);
            }
        }
    }

    return interpolated_E;
}

std::array<std::unique_ptr<MultiFab>, 3> WarpX::getInterpolatedB(int lev) const
{
    const int ngrow = 0;
    
    std::array<std::unique_ptr<MultiFab>, 3> interpolated_B;
    for (int i = 0; i < 3; ++i) {
        interpolated_B[i].reset( new MultiFab(Bfield_aux[lev][i]->boxArray(), dmap[lev], 1, ngrow) );
        interpolated_B[i]->setVal(0.0);
    }
    
    const Real* dx = Geom(lev-1).CellSize();
    const int r_ratio = refRatio(lev-1)[0];
    const int use_limiter = 0;
#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,3> bfab;
        for (MFIter mfi(*interpolated_B[0]); mfi.isValid(); ++mfi)
        {
            Box ccbx = mfi.fabbox();
            ccbx.enclosedCells();
            ccbx.coarsen(r_ratio).refine(r_ratio); // so that ccbx is coarsenable
            
            const FArrayBox& cxfab = (*Bfield_cp[lev][0])[mfi];
            const FArrayBox& cyfab = (*Bfield_cp[lev][1])[mfi];
            const FArrayBox& czfab = (*Bfield_cp[lev][2])[mfi];
            
            bfab[0].resize(amrex::convert(ccbx,Bx_nodal_flag));
            bfab[1].resize(amrex::convert(ccbx,By_nodal_flag));
            bfab[2].resize(amrex::convert(ccbx,Bz_nodal_flag));
            
#if (AMREX_SPACEDIM == 3)
            amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                         BL_TO_FORTRAN_ANYD(bfab[0]),
                                         BL_TO_FORTRAN_ANYD(bfab[1]),
                                         BL_TO_FORTRAN_ANYD(bfab[2]),
                                         BL_TO_FORTRAN_ANYD(cxfab),
                                         BL_TO_FORTRAN_ANYD(cyfab),
                                         BL_TO_FORTRAN_ANYD(czfab),
                                         dx, &r_ratio, &use_limiter);
#else
            amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                         BL_TO_FORTRAN_ANYD(bfab[0]),
                                         BL_TO_FORTRAN_ANYD(bfab[2]),
                                         BL_TO_FORTRAN_ANYD(cxfab),
                                         BL_TO_FORTRAN_ANYD(czfab),
                                         dx, &r_ratio, &use_limiter);
            amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                   BL_TO_FORTRAN_ANYD(bfab[1]),
                                   BL_TO_FORTRAN_ANYD(cyfab),
                                   &r_ratio, &use_limiter);
#endif

            for (int i = 0; i < 3; ++i) {
                const Box& bx = (*interpolated_B[i])[mfi].box();
                (*interpolated_B[i])[mfi].plus(bfab[i], bx, bx, 0, 0, 1);
            }
        }
    }
    return interpolated_B;
}

void
WarpX::PackPlotDataPtrs (Vector<const MultiFab*>& pmf,
                         const std::array<std::unique_ptr<MultiFab>,3>& data)
{
    BL_ASSERT(pmf.size() == AMREX_SPACEDIM);
#if (AMREX_SPACEDIM == 3)
    pmf[0] = data[0].get();
    pmf[1] = data[1].get();
    pmf[2] = data[2].get();
#elif (AMREX_SPACEDIM == 2)
    pmf[0] = data[0].get();
    pmf[1] = data[2].get();
#endif
}
