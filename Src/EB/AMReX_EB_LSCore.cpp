
#include <AMReX_Vector.H>
#include <AMReX_AmrCore.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <LSCore.H>
#include <LSCore_F.H>

//using namespace amrex;
namespace amrex {

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
LSCore::LSCore () {
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    level_set.resize(nlevs_max);
    ls_factory.resize(nlevs_max);

    istep.resize(nlevs_max, 0);

    // // periodic boundaries
    // int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    // int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs.setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs.setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }
}


LSCore::~LSCore () {}


// Initializes multilevel data
void LSCore::InitData () {
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
}


// Make a new level using provided BoxArray and DistributionMapping and fill
// with interpolated coarse level data. Overrides the pure virtual function in
// AmrCore
void LSCore::MakeNewLevelFromCoarse ( int lev, Real time, const BoxArray & ba,
                                      const DistributionMapping & dm) {
    const int ncomp  = level_set[lev - 1].nComp();
    const int nghost = level_set[lev - 1].nGrow();

    level_set[lev].define(ba, dm, ncomp, nghost);

    FillCoarsePatch(lev, time, level_set[lev], 0, ncomp);
}


// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void LSCore::RemakeLevel ( int lev, Real time, const BoxArray & ba,
                           const DistributionMapping & dm) {
    const int ncomp  = level_set[lev].nComp();
    const int nghost = level_set[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, level_set[lev]);
}


// Delete level data overrides the pure virtual function in AmrCore
void LSCore::ClearLevel (int lev) {
    level_set[lev].clear();
    ls_factory[lev].clear();
}


// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void LSCore::MakeNewLevelFromScratch (int lev, Real time, const BoxArray & ba,
                                      const DistributionMapping & dm) {
    const int ncomp  = 1;
    const int nghost = 1;

    level_set[lev].define(ba, dm, ncomp, nghost);
    // NOTE: use move semantics since LSFactory contains unique_ptrs
    ls_factory[lev] = std::move(LSFacory(0, ));

    const Real * dx = geom[lev].CellSize();

    // Number of particles
    int np = particle_data.size();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(level_set[lev], true); mfi.isValid(); ++mfi) {
        const Box &  tile_box = mfi.tilebox();
              auto & phi_tile = level_set[lev][mfi];

        fill_levelset_ib ( tile_box.loVect(), tile_box.hiVect(),
                           particle_data.dataPtr(), & np,
                           BL_TO_FORTRAN_3D(phi_tile),
                           dx                                     );

    }

    level_set[lev].FillBoundary(geom[lev].periodicity());

}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void LSCore::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) {
    static bool first = true;
    static Vector<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first) {
        first = false;
        // read in an array of "phierr", which is the tagging threshold in this
        // example, we tag values of "phi" which are greater than phierr for
        // that particular level in subroutine state_error, you could use more
        // elaborate tagging, such as more advanced logical expressions, or
        // gradients, etc.
        ParmParse pp("amr");
        int n = pp.countval("phierr");
        if (n > 0) {
            pp.getarr("phierr", phierr, 0, n);
        }
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real * dx      = geom[lev].CellSize();
    const Real * prob_lo = geom[lev].ProbLo();

    const MultiFab & state = level_set[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

        for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
        {
            const Box &    tilebox = mfi.tilebox();
                  TagBox & tagfab  = tags[mfi];

            // We cannot pass tagfab to Fortran because it is BaseFab<char>. So
            // we are going to get a temporary integer array. set itags
            // initially to 'untagged' everywhere we define itags over the
            // tilebox region
            tagfab.get_itags(itags, tilebox);

            // data pointer and index space
            int *       tptr = itags.dataPtr();
            const int * tlo  = tilebox.loVect();
            const int * thi  = tilebox.hiVect();

            // tag cells for refinement
            state_error ( tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
                          BL_TO_FORTRAN_3D(state[mfi]),
                          &tagval, &clearval,
                          AMREX_ARLIM_3D(tilebox.loVect()), AMREX_ARLIM_3D(tilebox.hiVect()),
                          AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &phierr[lev]);
            //
            // Now update the tags in the TagBox in the tilebox region to be
            // equal to itags
            //
            tagfab.tags_and_untags(itags, tilebox);
        }
    }
}

// read in some parameters from inputs file
void LSCore::ReadParameters () {


    /************************************************************************
     * Parse inputs                                                         *
     ***********************************************************************/


    {
        ParmParse pp("particles");

        int n_part;
        Vector<Real> pos, rad;

        pp.query("n_particles", n_part);
        pp.getarr("positions", pos, 0, 3 * n_part);
        pp.getarr("radii", rad, 0, n_part);


        // Convert input arrays to vectors of type particle_info
        int i_part = 0;
        for(int i = 0; i + 2 < pos.size(); i += 3) {
            particle_info particle;

            particle.pos = RealVect(pos[i], pos[i+1], pos[i+2]);
            particle.radius = rad[i_part];
            particle_data.push_back(particle);
            i_part ++;
        }
    }
}


// Set covered coarse cells to be the average of overlying fine cells
void LSCore::AverageDown () {
    for (int lev = finest_level-1; lev >= 0; --lev) {
        amrex::average_down(level_set[lev+1], level_set[lev],
                            geom[lev+1], geom[lev],
                            0, level_set[lev].nComp(), refRatio(lev));
    }
}


// More flexible version of AverageDown() that lets you average down across
// multiple levels
void LSCore::AverageDownTo (int crse_lev) {
    amrex::average_down(level_set[crse_lev+1], level_set[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, level_set[crse_lev].nComp(), refRatio(crse_lev));
}


// Compute a new multifab by coping in phi from valid region and filling ghost
// cells works for single level and 2-level cases (fill fine grid ghost by
// interpolating from coarse)
void LSCore::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp) {
    if (lev == 0) {

        PhysBCFunct physbc(geom[lev], bcs, BndryFunctBase(phifill));

        // NOTE: if source MultiFab vector as size = 1 => no interpolation
        amrex::FillPatchSingleLevel(mf, time, {& level_set[0]}, {0.}, 0, icomp, ncomp,
                                    geom[lev], physbc);

    } else {

        PhysBCFunct cphysbc(geom[lev-1], bcs, BndryFunctBase(phifill));
        PhysBCFunct fphysbc(geom[lev  ], bcs, BndryFunctBase(phifill));

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, {& level_set[lev - 1]}, {0.}, {& level_set[lev]}, {0.},
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, fphysbc, refRatio(lev-1),
                                  mapper, bcs);

    }
}


// Fill an entire multifab by interpolating from the coarser level. This comes
// into play when a new level of refinement appears
void LSCore::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp) {
    BL_ASSERT(lev > 0);

    PhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
    PhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

    Interpolater * mapper = & cell_cons_interp;

    amrex::InterpFromCoarseLevel(mf, time, level_set[lev - 1], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                 cphysbc, fphysbc, refRatio(lev-1),
                                 mapper, bcs);
}


// get plotfile name
std::string LSCore::PlotFileName (int lev) const {
    return amrex::Concatenate(plot_file, lev, 5);
}


// put together an array of multifabs for writing
Vector<const MultiFab*> LSCore::PlotFileMF () const {
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(&level_set[i]);
    }
    return r;
}


// set plotfile variable names
Vector<std::string> LSCore::PlotFileVarNames () const {
    return {"phi"};
}


// write plotfile to disk
void LSCore::WritePlotFile () const {
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                   Geom(), 0., istep, refRatio());
}


void LSCore::WriteCheckpointFile () const {

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g.,
    //                     finest_level, t_new, etc.) and also the BoxArrays at
    //                     each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at
    //                     each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (ParallelDescriptor::IOProcessor()) {

        std::string HeaderFileName(checkpointname + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for AmrCoreAdv\n";

        // write out finest_level
        HeaderFile << finest_level << "\n";

        // write out array of istep
        for (int i = 0; i < istep.size(); ++i) {
            HeaderFile << istep[i] << " ";
        }
        HeaderFile << "\n";

        // write the BoxArray at each level
        for (int lev = 0; lev <= finest_level; ++lev) {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Write(level_set[lev],
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
    }

}


void LSCore::ReadCheckpointFile () {

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = 1;
        int nghost = 0;
        level_set[lev].define(grids[lev], dmap[lev], ncomp, nghost);
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(level_set[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}


// utility to skip to next line in Header
void LSCore::GotoNextLine (std::istream & is) {
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}



}
