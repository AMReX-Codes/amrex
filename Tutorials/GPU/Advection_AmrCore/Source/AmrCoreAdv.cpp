
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Kernels_3d.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <AmrCoreAdv.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreAdv::AmrCoreAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

/*
    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
*/

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
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

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge = new amrex::AmrMeshInSituBridge;
    insitu_bridge->initialize();
#endif
}

AmrCoreAdv::~AmrCoreAdv ()
{
#ifdef BL_USE_SENSEI_INSITU
    delete insitu_bridge;
#endif
}

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration);

	cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

	// sync up time
	for (lev = 0; lev <= finest_level; ++lev) {
	    t_new[lev] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();
	}

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef BL_USE_SENSEI_INSITU
        insitu_bridge->update(step, cur_time,
            static_cast<amrex::AmrMesh*>(this), {&phi_new}, {{"phi"}});
#endif

#ifdef BL_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

	if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
	WritePlotFile();
    }

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge->finalize();
#endif
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
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

// Make a new level using provided BoxArray and DistributionMapping and 
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				    const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1].nComp();
    const int nghost = phi_new[lev-1].nGrow();
    
    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, phi_new[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and 
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
			 const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev].nComp();
    const int nghost = phi_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }    
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 0;

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    Real cur_time = t_new[lev];
    MultiFab& state = phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        Array4<Real> fab = state[mfi].array();
        GeometryData geomData = geom[lev].data();
        const Box& box = mfi.validbox();

        amrex::launch(box,
        [=] AMREX_GPU_DEVICE (Box const& tbx)
        {
            initdata(tbx, fab, geomData);
        });
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Vector<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
	first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
	ParmParse pp("adv");
	int n = pp.countval("phierr");
	if (n > 0) {
	    pp.getarr("phierr", phierr, 0, n);
	}
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;
	
	for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box& tilebox  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];
	    
	    // We cannot pass tagfab to Fortran because it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
            // set itags initially to 'untagged' everywhere
            // we define itags over the tilebox region
	    tagfab.get_itags(itags, tilebox);
	    
            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebox.loVect();
	    const int*  thi     = tilebox.hiVect();
/*
            state_error(tagfab, tilebox, state[mfi], &tagval, &clearval,
                        geomdata, &time, &phierr[lev]);
*/
            // tag cells for refinement
	    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
			BL_TO_FORTRAN_3D(state[mfi]),
			&tagval, &clearval, 
			AMREX_ARLIM_3D(tilebox.loVect()), AMREX_ARLIM_3D(tilebox.hiVect()), 
			AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &phierr[lev]);
	    //
	    // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
	    //
	    tagfab.tags_and_untags(itags, tilebox);
	}
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("regrid_int", regrid_int);
	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);
	pp.query("chk_file", chk_file);
	pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
    }

    {
	ParmParse pp("adv");
	
	pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(phi_new[lev+1], phi_new[lev],
                            geom[lev+1], geom[lev],
                            0, phi_new[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(phi_new[crse_lev+1], phi_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, phi_new[crse_lev].nComp(), refRatio(crse_lev));
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
	Vector<MultiFab*> smf;
	Vector<Real> stime;
	GetData(0, time, smf, stime);

        BndryFuncArray bfunc(phifill);
	PhysBCFunct<BndryFuncArray> physbc(geom[lev],bcs,bfunc);
	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
				     geom[lev], physbc, 0);
    }
    else
    {
	Vector<MultiFab*> cmf, fmf;
	Vector<Real> ctime, ftime;
	GetData(lev-1, time, cmf, ctime);
	GetData(lev  , time, fmf, ftime);

        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

	Interpolater* mapper = &cell_cons_interp;

	amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
				   0, icomp, ncomp, geom[lev-1], geom[lev],
				   cphysbc, 0, fphysbc, 0, refRatio(lev-1),
				   mapper, bcs, 0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    
    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    BndryFuncArray bfunc(phifill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

    Interpolater* mapper = &cell_cons_interp;

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, 0, fphysbc, 0, refRatio(lev-1),
				 mapper, bcs, 0);
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
	data.push_back(&phi_new[lev]);
	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
	data.push_back(&phi_old[lev]);
	datatime.push_back(t_old[lev]);
    }
    else
    {
	data.push_back(&phi_old[lev]);
	data.push_back(&phi_new[lev]);
	datatime.push_back(t_old[lev]);
	datatime.push_back(t_new[lev]);
    }
}


// advance a level by dt
// includes a recursive call for finer levels
void
AmrCoreAdv::timeStep (int lev, Real time, int iteration)
{
    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if 
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) 
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
		int old_finest = finest_level; 
		regrid(lev, time);

                // mark that we have regridded this level already
		for (int k = lev; k <= finest_level; ++k) {
		    last_regrid_step[k] = istep[k];
		}

                // if there are newly created levels, set the time step
		for (int k = old_finest+1; k <= finest_level; ++k) {
		    dt[k] = dt[k-1] / MaxRefRatio(k-1);
		}
	    }
	}
    }

    if (Verbose()) {
	amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
	amrex::Print() << "ADVANCE with time = " << t_new[lev] 
                       << " dt = " << dt[lev] << std::endl;
    }

    // advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose())
    {
	amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
	for (int i = 1; i <= nsubsteps[lev+1]; ++i)
	{
	    timeStep(lev+1, time+(i-1)*dt[lev+1], i);
	}

	if (do_reflux)
	{
            // update lev based on coarse-fine flux mismatch
	    flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);
	}

	AverageDownTo(lev); // average lev+1 down to lev
    }
    
}

// advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::Advance (int lev, Real time, Real dt_lev, int iteration, int ncycle)
{
    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    MultiFab& S_new = phi_new[lev];

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const auto dx = geom[lev].CellSizeArray();
    GpuArray<Real, AMREX_SPACEDIM> dtdx;
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        dtdx[i] = dt_lev/(dx[i]);
    }

    const Real* prob_lo = geom[lev].ProbLo();

    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux)
    {
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    BoxArray ba = grids[lev];
	    ba.surroundingNodes(i);
	    fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
	}
    }

    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());


/*
    // Allocate fabs for fluxes and Godunov velocities. (Kept for reference).
    for (int i = 0; i < BL_SPACEDIM ; i++) {
	const Box& bxtmp = amrex::surroundingNodes(bx,i);
	flux[i].resize(bxtmp,S_new.nComp());
	uface[i].resize(amrex::grow(bxtmp,1),1);
    }
*/

    // Build temporary multiFabs to work on.
    Array<MultiFab, AMREX_SPACEDIM> fluxcalc;
    Array<MultiFab, AMREX_SPACEDIM> facevel;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        BoxArray ba = amrex::convert(S_new.boxArray(), IntVect::TheDimensionVector(idim));

        fluxcalc[idim].define (ba,         S_new.DistributionMap(), S_new.nComp(), 0);
        facevel [idim].define (ba.grow(1), S_new.DistributionMap(),             1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{

        // ======== GET FACE VELOCITY =========
            GpuArray<Box, AMREX_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);,
                         nbx[1] = mfi.nodaltilebox(1);,
                         nbx[2] = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0),1);,
                         const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1),1);,
                         const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2),1););

            GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( facevel[0].array(mfi),
                                                                      facevel[1].array(mfi),
                                                                      facevel[2].array(mfi)) };

            const Box& psibox = Box(IntVect(AMREX_D_DECL(std::min(ngbxx.smallEnd(0)-1, ngbxy.smallEnd(0)-1),
                                                         std::min(ngbxx.smallEnd(1)-1, ngbxy.smallEnd(0)-1),
                                                         0)),
                                    IntVect(AMREX_D_DECL(std::max(ngbxx.bigEnd(0),   ngbxy.bigEnd(0)+1),
                                                         std::max(ngbxx.bigEnd(1)+1, ngbxy.bigEnd(1)),
                                                         0)));

            AsyncFab psifab(psibox, 1);
            Array4<Real> psi = psifab.array();
            GeometryData geomdata = geom[lev].data();
            auto prob_lo = geom[lev].ProbLoArray();
            auto dx = geom[lev].CellSizeArray();

            amrex::launch(psibox,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                get_face_velocity_psi(tbx, ctr_time,
                                      psi, geomdata); 
            });

            AMREX_D_TERM(
                         amrex::ParallelFor(ngbxx,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_x(i, j, k, vel[0], psi, prob_lo, dx); 
                         });,

                         amrex::ParallelFor(ngbxy,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_y(i, j, k, vel[1], psi, prob_lo, dx);
                         });,

                         amrex::ParallelFor(ngbxz,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_z(i, j, k, vel[2], psi, prob_lo, dx);
                         });
                        );

            psifab.clear();

        // ======== FLUX CALC AND UPDATE =========

	    const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> statein  = Sborder.array(mfi);
            Array4<Real> stateout = S_new.array(mfi);

            GpuArray<Array4<Real>, AMREX_SPACEDIM> flux{ AMREX_D_DECL(fluxcalc[0].array(mfi),
                                                                      fluxcalc[1].array(mfi),
                                                                      fluxcalc[2].array(mfi)) };

            AMREX_D_TERM(const Box& dqbxx = amrex::grow(bx, IntVect{2, 1, 1});,
                         const Box& dqbxy = amrex::grow(bx, IntVect{1, 2, 1});,
                         const Box& dqbxz = amrex::grow(bx, IntVect{1, 1, 2}););

            AsyncFab slope2fab (amrex::grow(bx, 2), 1);
            Array4<Real> slope2 = slope2fab.array();
            AsyncFab slope4fab (amrex::grow(bx, 1), 1);
            Array4<Real> slope4 = slope4fab.array();

            // compute longitudinal fluxes
            // ===========================

            // x -------------------------
            AsyncFab phixfab (gbx, 1);
            Array4<Real> phix = phixfab.array();

            amrex::launch(dqbxx,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                slopex2(tbx, statein, slope2);
            });

            amrex::launch(gbx,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                slopex4(tbx, statein, slope2, slope4);
            });

            amrex::ParallelFor(amrex::growLo(gbx, 0, -1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_x(i, j, k, statein, vel[0], phix, slope4, dtdx); 
            });

            // y -------------------------
            AsyncFab phiyfab (gbx, 1);
            Array4<Real> phiy = phiyfab.array();

            amrex::launch(dqbxy,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                slopey2(tbx, statein, slope2);
            });

            amrex::launch(gbx,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                slopey4(tbx, statein, slope2, slope4);
            });

            amrex::ParallelFor(amrex::growLo(gbx, 1, -1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_y(i, j, k, statein, vel[1], phiy, slope4, dtdx); 
            });

            // z -------------------------
            AsyncFab phizfab (gbx, 1);
            Array4<Real> phiz = phizfab.array();

            amrex::launch(dqbxz,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                slopez2(tbx, statein, slope2);
            });

            amrex::launch(gbx,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                slopez4(tbx, statein, slope2, slope4);
            });

            amrex::ParallelFor(amrex::growLo(gbx, 2, -1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_z(i, j, k, statein, vel[2], phiz, slope4, dtdx); 
            });

            slope2fab.clear();
            slope4fab.clear();

            // compute transverse fluxes
            // ===========================

            AMREX_D_TERM(const Box& gbxx = amrex::grow(bx, 0, 1);,
                         const Box& gbxy = amrex::grow(bx, 1, 1);,
                         const Box& gbxz = amrex::grow(bx, 2, 1););

            // xy & xz --------------------
            AsyncFab phix_yfab (gbx, 1);
            AsyncFab phix_zfab (gbx, 1);
            Array4<Real> phix_y = phix_yfab.array();
            Array4<Real> phix_z = phix_zfab.array();

            amrex::ParallelFor(amrex::growHi(gbxz, 0, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_xy(i, j, k, 
                        AMREX_D_DECL(vel[0], vel[1], vel[2]),
                        AMREX_D_DECL(phix, phiy, phiz),
                        phix_y, dtdx);
            }); 

            amrex::ParallelFor(amrex::growHi(gbxy, 0, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_xz(i, j, k,
                        AMREX_D_DECL(vel[0], vel[1], vel[2]),
                        AMREX_D_DECL(phix, phiy, phiz),
                        phix_z, dtdx);
            }); 

            // yz & yz --------------------
            AsyncFab phiy_xfab (gbx, 1);
            AsyncFab phiy_zfab (gbx, 1);
            Array4<Real> phiy_x = phiy_xfab.array();
            Array4<Real> phiy_z = phiy_zfab.array();

            amrex::ParallelFor(amrex::growHi(gbxz, 1, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_yx(i, j, k,
                        AMREX_D_DECL(vel[0], vel[1], vel[2]),
                        AMREX_D_DECL(phix, phiy, phiz),
                        phiy_x, dtdx);
            }); 

            amrex::ParallelFor(amrex::growHi(gbxx, 1, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_yz(i, j, k,
                        AMREX_D_DECL(vel[0], vel[1], vel[2]),
                        AMREX_D_DECL(phix, phiy, phiz),
                        phiy_z, dtdx);
            }); 

            // zx & zy --------------------
            AsyncFab phiz_xfab (gbx, 1);
            AsyncFab phiz_yfab (gbx, 1);
            Array4<Real> phiz_x = phiz_xfab.array();
            Array4<Real> phiz_y = phiz_yfab.array();

            amrex::ParallelFor(amrex::growHi(gbxy, 2, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_zx(i, j, k, 
                        AMREX_D_DECL(vel[0], vel[1], vel[2]),
                        AMREX_D_DECL(phix, phiy, phiz),
                        phiz_x, dtdx);
            }); 

            amrex::ParallelFor(amrex::growHi(gbxx, 2, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                flux_zy(i, j, k,
                        AMREX_D_DECL(vel[0], vel[1], vel[2]),
                        AMREX_D_DECL(phix, phiy, phiz),
                        phiz_y, dtdx);
            }); 

            // final edge states 
            // ===========================
            amrex::ParallelFor(amrex::growHi(bx, 0, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                combine_flux_x(i, j, k,
                               vel[0], vel[1], vel[2],
                               phix, phiy_z, phiz_y,
                               flux[0], dtdx);
            });

            amrex::ParallelFor(amrex::growHi(bx, 1, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                combine_flux_y(i, j, k,
                               vel[0], vel[1], vel[2],
                               phiy, phix_z, phiz_x,
                               flux[1], dtdx);
            });

            amrex::ParallelFor(amrex::growHi(bx, 2, 1),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                combine_flux_z(i, j, k,
                               vel[0], vel[1], vel[2],
                               phiz, phix_y, phiy_x,
                               flux[2], dtdx);
            });

            // Flux has been updated. These temporaries are no longer needed.
            phixfab.clear();
            phiyfab.clear();
            phizfab.clear();
            phix_yfab.clear();
            phix_zfab.clear();
            phiy_xfab.clear();
            phiy_zfab.clear();
            phiz_xfab.clear();
            phiz_yfab.clear();

            // compute new state (stateout) and scale fluxes based on face area.
            // ===========================

            // Do a conservative update 
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                conservative(i, j, k,
                             statein, stateout,
                             AMREX_D_DECL(flux[0], flux[1], flux[2]),
                             dtdx);
            });

            // Scale by face area in order to correctly reflux
            AMREX_D_TERM(
                         amrex::ParallelFor(amrex::growHi(bx, 0, 1),
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             flux_scale_x(i, j, k, flux[0], dt_lev, dx);
                         });,
 
                         amrex::ParallelFor(amrex::growHi(bx, 1, 1),
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             flux_scale_y(i, j, k, flux[1], dt_lev, dx);
                         });,

                         amrex::ParallelFor(amrex::growHi(bx, 2, 1),
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             flux_scale_z(i, j, k, flux[2], dt_lev, dx);
                         });
                        );

            GpuArray<Array4<Real>, AMREX_SPACEDIM> fluxout{ AMREX_D_DECL(fluxes[0].array(mfi),
                                                                         fluxes[1].array(mfi),
                                                                         fluxes[2].array(mfi)) };
          
            if (do_reflux) {
                for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
                    amrex::ParallelFor(nbx[idim],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        fluxout[idim](i,j,k) = flux[idim](i,j,k);
                    });
                }
            }
        }
    }

    // ======== CFL CHECK, MOVED OUTSIDE MFITER LOOP =========

    AMREX_D_TERM(Real umax = facevel[0].norm0(0,0,false);,
                 Real vmax = facevel[1].norm0(0,0,false);,
                 Real wmax = facevel[2].norm0(0,0,false););

    if (AMREX_D_TERM(umax*dt_lev > dx[0], ||
                     vmax*dt_lev > dx[1], ||
                     wmax*dt_lev > dx[2]))
    {
        amrex::Print() << "umax = " << umax << ", vmax = " << vmax << ", wmax = " << wmax 
                       << ", dt = " << ctr_time << " dx = " << dx[1] << " " << dx[2] << " " << dx[3] << std::endl;
        amrex::Abort("CFL violation. use smaller adv.cfl.");
    }

    // ======== END OF GPU EDIT, (FOR NOW) =========

    // increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // In this example we are solving phi_t = -div(+F)
    // The fluxes contain, e.g., F_{i+1/2,j} = (phi*u)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (do_reflux) { 
	if (flux_reg[lev+1]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
	        // update the lev+1/lev flux register (index lev+1)   
	        flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
	    }	    
	}
	if (flux_reg[lev]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
	        // update the lev/lev-1 flux register (index lev) 
		flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
	    }
	}
    }
}

// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
	dt_tmp[lev] = EstTimeStep(lev, true);
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
	dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
	n_factor *= nsubsteps[lev];
	dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
	dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
	dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, bool local) const
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
//    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];
    const MultiFab& S_new = phi_new[lev];

    Array<MultiFab,AMREX_SPACEDIM> facevel;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        facevel[idim].define(amrex::convert(S_new.boxArray(), IntVect::TheDimensionVector(idim)),
                             S_new.DistributionMap(), 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est) if (Gpu::notInLaunchRegion())
#endif
    {
        // Calculate face velocities.
	for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
            AMREX_D_TERM(const Box& nbxx = mfi.nodaltilebox(0);,
                         const Box& nbxy = mfi.nodaltilebox(1);,
                         const Box& nbxz = mfi.nodaltilebox(2););

            GpuArray<Array4<Real>, AMREX_SPACEDIM> vel { AMREX_D_DECL(facevel[0].array(mfi),
                                                                      facevel[1].array(mfi),
                                                                      facevel[2].array(mfi)) };

            const Box& psibox = Box(IntVect(AMREX_D_DECL(std::min(nbxx.smallEnd(0)-1, nbxy.smallEnd(0)-1),
                                                         std::min(nbxx.smallEnd(1)-1, nbxy.smallEnd(0)-1),
                                                         0)),
                                    IntVect(AMREX_D_DECL(std::max(nbxx.bigEnd(0),     nbxy.bigEnd(0)+1),
                                                         std::max(nbxx.bigEnd(1)+1,   nbxy.bigEnd(1)),
                                                         0)));

            AsyncFab psifab(psibox, 1);
            Array4<Real> psi = psifab.array();
            GeometryData geomdata = geom[lev].data();
            auto prob_lo = geom[lev].ProbLoArray();
            auto dx = geom[lev].CellSizeArray();

            amrex::launch(psibox, 
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                get_face_velocity_psi(tbx, cur_time, psi, geomdata); 
            });

            AMREX_D_TERM(
                         amrex::ParallelFor(nbxx,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_x(i, j, k, vel[0], psi, prob_lo, dx); 
                         });,

                         amrex::ParallelFor(nbxy,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_y(i, j, k, vel[1], psi, prob_lo, dx);
                         });,

                         amrex::ParallelFor(nbxz,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_z(i, j, k, vel[2], psi, prob_lo, dx);
                         });
                        );

            psifab.clear();
	}
    }

    for (int i=0; i<BL_SPACEDIM; ++i)
    {
        Real est = facevel[i].norm0(0,0,true);
        dt_est = std::min(dt_est, dx[i]/est);
    }

    // Currently, this never happens (function called with local = true).
    // Reduction occurs outside this function.
    if (!local) {
	ParallelDescriptor::ReduceRealMin(dt_est);
    }

    dt_est *= cfl;

    return dt_est;
}

// get plotfile name
std::string
AmrCoreAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreAdv::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(&phi_new[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreAdv::PlotFileVarNames () const
{
    return {"phi"};
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				   Geom(), t_new[0], istep, refRatio());
}

void
AmrCoreAdv::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

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

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
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
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

}


void
AmrCoreAdv::ReadCheckpointFile ()
{

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

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
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
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}

// utility to skip to next line in Header
void
AmrCoreAdv::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
