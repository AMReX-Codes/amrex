
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

#include <AMReX_EB_LSCoreBase.H>

//using namespace amrex;
namespace amrex {

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
LSCoreBase::LSCoreBase () {
    //NOTE: Geometry on all levels has been defined already.

    ReadParameters();
    InitLSCoreBase();
}


LSCoreBase::LSCoreBase(const RealBox * rb, int max_level_in, const Vector<int> & n_cell_in, int coord)
    : AmrCore(rb, max_level_in, n_cell_in, coord)
{
    //NOTE: Geometry on all levels has been defined already.

    ReadParameters();
    InitLSCoreBase();
}


LSCoreBase::~LSCoreBase () {}


void LSCoreBase::InitLSCoreBase() {

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    level_set.resize(nlevs_max);
    level_set_valid.resize(nlevs_max);

    ls_factory.resize(nlevs_max);
    eb_levels.resize(nlevs_max);
    rebuild_eb.resize(nlevs_max, 1); // At first rebuild eb on each level

    bcs.resize(1);

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
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }

}

void LSCoreBase::LoadTagLevels () {
    // read in an array of "phierr", which is the tagging threshold in this
    // example, we tag values of "phi" which are greater than phierr for that
    // particular level in subroutine state_error, you could use more elaborate
    // tagging, such as more advanced logical expressions, or gradients, etc.
    ParmParse pp("amr");
    int n = pp.countval("phierr");
    if (n > 0) {
        pp.getarr("phierr", phierr, 0, n);
    }
}

void LSCoreBase::SetTagLevels (const Vector<Real> & m_phierr) {
    phierr = m_phierr;
}


// Initializes multilevel data
void LSCoreBase::Init () {
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;

        // This tells the AmrMesh class not to iterate when creating the 
        //    initial grid hierarchy
        SetIterateToFalse();

        // This tells the Cluster routine to use the new chopping
        // routine which rejects cuts if they don't improve the efficiency
        SetUseNewChop();

        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0)
            WriteCheckpointFile();

    } else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0)
        WritePlotFile();
}


void LSCoreBase::InitData (bool a_use_phierr) {
    use_phierr = a_use_phierr;
    if (use_phierr)
        LoadTagLevels();
    Init();
}


void LSCoreBase::InitData (const Vector<Real> & m_phierr) {
    SetTagLevels(m_phierr);
    Init();
}


// Make a new level using provided BoxArray and DistributionMapping and fill
// with interpolated coarse level data. Overrides the pure virtual function in
// AmrCore
void LSCoreBase::MakeNewLevelFromCoarse ( int lev, Real time, const BoxArray & ba,
                                          const DistributionMapping & dm) {
    const int ncomp  = level_set[lev - 1].nComp();
    const int nghost = level_set[lev - 1].nGrow();

    BoxArray ba_nd = amrex::convert(ba, IntVect{1, 1, 1});
    level_set[lev].define(ba_nd, dm, ncomp, nghost);

    FillCoarsePatch(lev, time, level_set[lev], 0, ncomp);

    // At this point, we consider _everywhere_ as valid. This is maintained for
    // legacy reasons. TODO: There might be a better way of doing things.
    level_set_valid[lev].define(ba_nd, dm, ncomp, nghost);
    level_set_valid[lev].setVal(1);
}


// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data. Overrides the pure virtual function
// in AmrCore
void LSCoreBase::RemakeLevel ( int lev, Real time, const BoxArray & ba,
                               const DistributionMapping & dm) {
    const int ncomp  = level_set[lev].nComp();
    const int nghost = level_set[lev].nGrow();

    BoxArray ba_nd = amrex::convert(ba, IntVect{1, 1, 1});
    MultiFab new_state(ba_nd, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, level_set[lev]);

    // At this point, we consider _everywhere_ as valid. This is maintained for
    // legacy reasons. TODO: There might be a better way of doing things.
    level_set_valid[lev].define(ba_nd, dm, ncomp, nghost);
    level_set_valid[lev].setVal(1);
}


void LSCoreBase::UpdateGrids (int lev, const BoxArray & ba, const DistributionMapping & dm){

    bool ba_changed = ( ba != grids[lev] );
    bool dm_changed = ( dm != dmap[lev] );

    if (! (ba_changed || dm_changed))
        return;


    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    BoxArray ba_nd = amrex::convert(ba, IntVect{AMREX_D_DECL(1, 1, 1)});

    MultiFab ls_regrid = MFUtil::duplicate<MultiFab, MFUtil::SymmetricGhost> (ba_nd, dm, level_set[lev]);
    iMultiFab valid_regrid = MFUtil::duplicate<iMultiFab, MFUtil::SymmetricGhost>
        (ba_nd, dm, level_set_valid[lev]);

    std::swap(ls_regrid, level_set[lev]);
    std::swap(valid_regrid, level_set_valid[lev]);
}



// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void LSCoreBase::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) {

    if (use_phierr && (lev >= phierr.size())) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real * dx      = geom[lev].CellSize();
    const Real * prob_lo = geom[lev].ProbLo();

    MultiFab volfrac(grids[lev], dmap[lev], 1, 1);
    const MultiFab * state = & level_set[lev];

    if (! use_phierr) {
        eb_levels[lev]->fillVolFrac( volfrac, geom[lev]);
        state = & volfrac;
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

        for (MFIter mfi(* state, true); mfi.isValid(); ++mfi) {
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
            if (use_phierr) {
                amrex_eb_levelset_error ( tptr, AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
                                          BL_TO_FORTRAN_3D((* state)[mfi]),
                                          & tagval, & clearval,
                                          BL_TO_FORTRAN_BOX(tilebox),
                                          AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo),
                                          & time, & phierr[lev]);

            } else {

                Real tol = 0.000001; // TODO: Fix magic numbers (after talking with ANN)
                amrex_eb_volfrac_error ( tptr, AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
                                         BL_TO_FORTRAN_3D((*state)[mfi]),
                                         & tagval, & clearval, & tol,
                                         BL_TO_FORTRAN_BOX(tilebox),
                                         AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo));
            }

            //___________________________________________________________________
            // Now update the tags in the TagBox in the tilebox region to be
            // equal to itags
            //
            if (! use_phierr)
                tagfab.buffer(8, 0); // TODO: Fix magic numbers (after talking with ANN)
            tagfab.tags_and_untags(itags, tilebox);
        }
    }
}


// Read some parameters from inputs file
void LSCoreBase::ReadParameters () {

    /************************************************************************
     * Parse inputs                                                         *
     ***********************************************************************/

    ParmParse pp("eb_amr");
    pp.query("eb_pad", eb_pad);
    pp.query("max_eb_pad", max_eb_pad);

}


// Set covered coarse cells to be the average of overlying fine cells
void LSCoreBase::AverageDown () {
    for (int lev = finest_level-1; lev >= 0; lev--) {

        amrex::average_down(level_set[lev + 1], level_set[lev],
                            0, level_set[lev].nComp(), refRatio(lev));
    }
}


// More flexible version of AverageDown() that lets you average down across
// multiple levels
void LSCoreBase::AverageDownTo (int crse_lev) {

    amrex::average_down(level_set[crse_lev+1], level_set[crse_lev],
                        0, level_set[crse_lev].nComp(), refRatio(crse_lev));

}


// Compute a new multifab by coping in phi from valid region and filling ghost
// cells works for single level and 2-level cases (fill fine grid ghost by
// interpolating from coarse)
void LSCoreBase::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp) {
    if (lev == 0) {

        BndryFuncArray bfunc(amrex_eb_phifill);
        PhysBCFunct<BndryFuncArray> physbc(geom[lev], bcs, bfunc);

        // NOTE: if source MultiFab vector as size = 1 => no interpolation
        amrex::FillPatchSingleLevel(mf, time, {& level_set[0]}, {0.}, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);

    } else {

        BndryFuncArray bfunc(amrex_eb_phifill);
        PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1], bcs, bfunc);
        PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ], bcs, bfunc);

        Interpolater * mapper = & node_bilinear_interp;


        amrex::FillPatchTwoLevels(mf, time, {& level_set[lev - 1]}, {0.}, {& level_set[lev]}, {0.},
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);

    }
}


// Fill an entire multifab by interpolating from the coarser level. This comes
// into play when a new level of refinement appears
void LSCoreBase::FillCoarsePatch (int lev, Real time, MultiFab & mf, int icomp, int ncomp) {
    BL_ASSERT(lev > 0);

    BndryFuncArray bfunc(amrex_eb_phifill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1], bcs, bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ], bcs, bfunc);

    Interpolater * mapper = & node_bilinear_interp;

    amrex::InterpFromCoarseLevel(mf, time, level_set[lev - 1], 0, icomp, ncomp,
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev - 1), mapper, bcs, 0);

}


// Constructs a box over which to look for EB facets. The Box size grows based
// on the coarse-level level-set value. But it never grows larger than
// max_eb_pad.
Box LSCoreBase::EBSearchBox(const Box & tilebox, const FArrayBox & ls_crse,
                            const Geometry & geom_fine, bool & bail) {

    // Infinities don't work well with std::max, so just bail and construct the
    // maximum box.
    if (ls_crse.contains_inf()){
        IntVect n_grow(AMREX_D_DECL(max_eb_pad, max_eb_pad, max_eb_pad));
        Box bx = amrex::convert(ls_crse.box(), IntVect{0, 0, 0});
        bx.grow(n_grow);

        bail = true;
        return bx;
    }

    // Something's gone wrong :( ... so just bail and construct the maximum box.
    if (ls_crse.contains_nan()){
        IntVect n_grow(AMREX_D_DECL(max_eb_pad, max_eb_pad, max_eb_pad));
        Box bx = amrex::convert(ls_crse.box(), IntVect{0, 0, 0});
        bx.grow(n_grow);

        bail = true;
        return bx;
    }


    Real max_ls = std::max(ls_crse.max(), std::abs(ls_crse.min()));

    IntVect n_grow(AMREX_D_DECL(geom_fine.InvCellSize(0)*max_ls,
                                geom_fine.InvCellSize(1)*max_ls,
                                geom_fine.InvCellSize(2)*max_ls));

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        if (n_grow[i] > max_eb_pad) {
            n_grow[i] = max_eb_pad;
            bail = true;
        }

    Box bx = amrex::convert(tilebox, IntVect{AMREX_D_DECL(0, 0, 0)});
    bx.grow(n_grow);

    return bx;
}



// Get plotfile name
std::string LSCoreBase::PlotFileName (int lev) const {
    // return amrex::Concatenate(plot_file, lev, 5);
    return plot_file;
}


// Put together an array of multifabs for writing
Vector<MultiFab> LSCoreBase::PlotFileMF () const {
    Vector<MultiFab> r(max_level + 1);
    for (int i = 0; i < max_level + 1; i++) {
        const int ncomp  = level_set[i].nComp();
        const int nghost = level_set[i].nGrow();
        r[i].define(grids[i], dmap[i], ncomp, nghost);

        amrex::average_node_to_cellcenter(r[i], 0, level_set[i], 0, 1);
    }
    return std::move(r);
}


// Set plotfile variable names
Vector<std::string> LSCoreBase::PlotFileVarNames () const {
    return {"level-set"};
}


// Write plotfile to disk
void LSCoreBase::WritePlotFile () const {
    // Get plotfile name
    const std::string & plotfilename = PlotFileName(0);

    // Generate cell-centered data to put into plotfile
    const Vector<MultiFab> mf_plt = PlotFileMF();
    Vector<const MultiFab*> mf_ptr;
    for (const MultiFab & mf : mf_plt)
        mf_ptr.push_back(& mf);

    // Get variable names
    const auto & varnames = PlotFileVarNames();

    // Keep user informed
    amrex::Print() << "Writing ";
    for (const std::string & str_name : varnames)
        amrex::Print() << str_name << " ";
    amrex::Print() << "plotfile: " << plotfilename << "\n";

    // Save plot file
    Vector<int> istep(max_level + 1, 0);
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, mf_ptr, varnames,
                                   Geom(), 0., istep, refRatio());
}


void LSCoreBase::WriteCheckpointFile () const {

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g.,
    //                     finest_level, t_new, etc.) and also the BoxArrays at
    //                     each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at
    //                     each level of refinement

    // checkpoint file name, e.g., chk00010
    // const std::string & checkpointname = amrex::Concatenate(chk_file,istep[0]);
    const std::string & checkpointname = chk_file;

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level + 1;

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
        Vector<int> istep(max_level + 1, 0);
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


void LSCoreBase::ReadCheckpointFile () {

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
void LSCoreBase::GotoNextLine (std::istream & is) {
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}



}
