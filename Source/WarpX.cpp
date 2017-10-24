
#include <limits>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_ParmParse.H>

#include <AMReX_MGT_Solver.H>
#include <AMReX_stencil_types.H>

#include <AMReX_MultiFabUtil.H>

#include <WarpX.H>
#include <WarpX_f.H>
#include <WarpXConst.H>

using namespace amrex;

Vector<Real> WarpX::B_external(3, 0.0);

Real WarpX::gamma_boost = 1.;
Real WarpX::beta_boost = 0.;
Vector<Real> WarpX::boost_direction(3, 0.0);

long WarpX::current_deposition_algo = 3;
long WarpX::charge_deposition_algo = 0;
long WarpX::field_gathering_algo = 1;
long WarpX::particle_pusher_algo = 0;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

bool WarpX::use_laser         = false;
bool WarpX::use_filter        = false;
bool WarpX::serialize_ics     = false;

bool WarpX::do_boosted_frame_diagnostic = false;
int  WarpX::num_snapshots_lab = std::numeric_limits<int>::lowest();
Real WarpX::dt_snapshots_lab  = std::numeric_limits<Real>::lowest();

#if (BL_SPACEDIM == 3)
IntVect WarpX::Bx_nodal_flag(1,0,0);
IntVect WarpX::By_nodal_flag(0,1,0);
IntVect WarpX::Bz_nodal_flag(0,0,1);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::Bx_nodal_flag(1,0);  // x is the first dimension to AMReX
IntVect WarpX::By_nodal_flag(0,0);  // y is the missing dimension to 2D AMReX
IntVect WarpX::Bz_nodal_flag(0,1);  // z is the second dimension to 2D AMReX
#endif

#if (BL_SPACEDIM == 3)
IntVect WarpX::Ex_nodal_flag(0,1,1);
IntVect WarpX::Ey_nodal_flag(1,0,1);
IntVect WarpX::Ez_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::Ex_nodal_flag(0,1);  // x is the first dimension to AMReX
IntVect WarpX::Ey_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
IntVect WarpX::Ez_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

#if (BL_SPACEDIM == 3)
IntVect WarpX::jx_nodal_flag(0,1,1);
IntVect WarpX::jy_nodal_flag(1,0,1);
IntVect WarpX::jz_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::jx_nodal_flag(0,1);  // x is the first dimension to AMReX
IntVect WarpX::jy_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
IntVect WarpX::jz_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

WarpX* WarpX::m_instance = nullptr;

WarpX&
WarpX::GetInstance ()
{
    if (!m_instance) {
	m_instance = new WarpX();
    }
    return *m_instance;
}

void
WarpX::ResetInstance ()
{
    delete m_instance;
    m_instance = nullptr;
}

WarpX::WarpX ()
{
    m_instance = this;

    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
#if 0
    // no subcycling yet
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }
#endif

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    // Particle Container
    mypc = std::unique_ptr<MultiParticleContainer> (new MultiParticleContainer(this));

    Efield_aux.resize(nlevs_max);
    Bfield_aux.resize(nlevs_max);

    F_fp.resize(nlevs_max);
    rho_fp.resize(nlevs_max);
    current_fp.resize(nlevs_max);
    Efield_fp.resize(nlevs_max);
    Bfield_fp.resize(nlevs_max);

    F_cp.resize(nlevs_max);
    rho_cp.resize(nlevs_max);
    current_cp.resize(nlevs_max);
    Efield_cp.resize(nlevs_max);
    Bfield_cp.resize(nlevs_max);

    pml.resize(nlevs_max);

    masks.resize(nlevs_max);
    gather_masks.resize(nlevs_max);

    costs.resize(nlevs_max);
}

WarpX::~WarpX ()
{
}

void
WarpX::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("check_file", check_file);
	pp.query("check_int", check_int);

	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);

	pp.query("restart", restart_chkfile);
    }

    {
	ParmParse pp("warpx");

	pp.query("cfl", cfl);
	pp.query("verbose", verbose);
	pp.query("regrid_int", regrid_int);

        // Boosted-frame parameters
        pp.query("gamma_boost", gamma_boost);
        beta_boost = std::sqrt(1.-1./pow(gamma_boost,2));
        pp.queryarr("boost_direction", boost_direction);
        if( gamma_boost > 1. ){
            // Read and normalize the boost direction
            Real s = 1.0/std::sqrt(boost_direction[0]*boost_direction[0] +
                                   boost_direction[1]*boost_direction[1] +
                                   boost_direction[2]*boost_direction[2]);
	    boost_direction = { boost_direction[0]*s,
                                boost_direction[1]*s,
                                boost_direction[2]*s };
        }
        
        pp.queryarr("B_external", B_external);

	pp.query("do_moving_window", do_moving_window);
	if (do_moving_window)
	{
	    std::string s;
	    pp.get("moving_window_dir", s);
	    if (s == "x" || s == "X") {
		moving_window_dir = 0;
	    }
#if (BL_SPACEDIM == 3)
	    else if (s == "y" || s == "Y") {
		moving_window_dir = 1;
	    }
#endif
	    else if (s == "z" || s == "Z") {
		moving_window_dir = BL_SPACEDIM-1;
	    }
	    else {
		const std::string msg = "Unknown moving_window_dir: "+s;
		amrex::Abort(msg.c_str());
	    }

	    moving_window_x = geom[0].ProbLo(moving_window_dir);

	    pp.get("moving_window_v", moving_window_v);
	    moving_window_v *= PhysConst::c;
	}

	pp.query("do_plasma_injection", do_plasma_injection);
	if (do_plasma_injection) {
	  pp.get("num_injected_species", num_injected_species);
	  injected_plasma_species.resize(num_injected_species);
	  pp.getarr("injected_plasma_species", injected_plasma_species,
		 0, num_injected_species);
	}

        pp.query("do_boosted_frame_diagnostic", do_boosted_frame_diagnostic);       
        if (do_boosted_frame_diagnostic) {
            
            if (gamma_boost <= 1.0) {
                amrex::Abort("gamma_boost must be > 1 to use the boosted frame diagnostic.");
            }
            
            if ( (std::abs(boost_direction[0] - 0.0) > 1.0e-12) or 
                 (std::abs(boost_direction[1] - 0.0) > 1.0e-12) or 
                 (std::abs(boost_direction[2] - 1.0) > 1.0e012)) {
                amrex::Abort("The boosted frame diagnostic currently only works if the boost is in the z direction.");
            }

            pp.get("num_snapshots_lab", num_snapshots_lab);
            pp.get("dt_snapshots_lab", dt_snapshots_lab);
            pp.get("gamma_boost", gamma_boost);
            if (not do_moving_window) {
                amrex::Abort("The moving window shoudd be on if using the boosted frame diagnostic");
            }
            std::string s;
	    pp.get("moving_window_dir", s);
	    if ( not (s == "z" || s == "Z")) {
                amrex::Abort("The boosted frame diagnostic currently only works if the window is moving in the z direction.");
            }                     
        }
        
        pp.query("do_electrostatic", do_electrostatic);
        pp.query("n_buffer", n_buffer);
        pp.query("const_dt", const_dt);

	pp.query("use_laser", use_laser);
	pp.query("use_filter", use_filter);
	pp.query("serialize_ics", serialize_ics);
        pp.query("do_dive_cleaning", do_dive_cleaning);

        pp.query("do_pml", do_pml);
        pp.query("pml_ncell", pml_ncell);
        pp.query("pml_delta", pml_delta);
        pp.query("pml_type", pml_type);

        pp.query("plot_raw_fields", plot_raw_fields);
        pp.query("plot_raw_fields_guards", plot_raw_fields_guards);
        if (ParallelDescriptor::NProcs() == 1) {
            plot_proc_number = false;
        }
        pp.query("plot_part_per_cell", plot_part_per_cell);
        pp.query("plot_part_per_grid", plot_part_per_grid);
        pp.query("plot_part_per_proc", plot_part_per_proc);
        pp.query("plot_proc_number"  , plot_proc_number);
        pp.query("plot_dive"         , plot_dive);
        pp.query("plot_divb"         , plot_divb);

        if (maxLevel() > 0) {
            pp.query("plot_finepatch", plot_finepatch);
            pp.query("plot_crsepatch", plot_crsepatch);
        }

        if (maxLevel() > 0) {
            Vector<Real> lo, hi;
            pp.getarr("fine_tag_lo", lo);
            pp.getarr("fine_tag_hi", hi);
            fine_tag_lo = RealVect{lo};
            fine_tag_hi = RealVect{hi};
        }

        pp.query("load_balance_int", load_balance_int);
    }

    {
	ParmParse pp("interpolation");
	pp.query("nox", nox);
	pp.query("noy", noy);
	pp.query("noz", noz);
	if (nox != noy || nox != noz) {
	    amrex::Abort("warpx.nox, noy and noz must be equal");
	}
	if (nox < 1) {
	    amrex::Abort("warpx.nox must >= 1");
	}
    }
    
    {
	ParmParse pp("algo");
	pp.query("current_deposition", current_deposition_algo);
	pp.query("charge_deposition", charge_deposition_algo);
	pp.query("field_gathering", field_gathering_algo);
	pp.query("particle_pusher", particle_pusher_algo);
    }
}

// This is a virtual function.
void
WarpX::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                const DistributionMapping& new_dmap)
{
    AllocLevelData(lev, new_grids, new_dmap);
    InitLevelData(lev, time);
}

void
WarpX::ClearLevel (int lev)
{
    for (int i = 0; i < 3; ++i) {
	Efield_aux[lev][i].reset();
	Bfield_aux[lev][i].reset();

	current_fp[lev][i].reset();
	Efield_fp [lev][i].reset();
	Bfield_fp [lev][i].reset();

	current_cp[lev][i].reset();
	Efield_cp [lev][i].reset();
	Bfield_cp [lev][i].reset();
    }

    F_fp  [lev].reset();
    rho_fp[lev].reset();
    F_cp  [lev].reset();
    rho_cp[lev].reset();

    costs[lev].reset();
}

void
WarpX::AllocLevelData (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    // WarpX assumes the same number of guard cells for Ex, Ey, Ez, Bx, By, Bz
    int ngE   = (WarpX::nox % 2) ? WarpX::nox+1 : WarpX::nox;  // Always even number
    int ngJ = ngE;
    int ngRho = ngE;

    //
    // The fine patch
    //
    Bfield_fp[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,1,ngE));
    Bfield_fp[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,1,ngE));
    Bfield_fp[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,1,ngE));

    Efield_fp[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,1,ngE));
    Efield_fp[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,1,ngE));
    Efield_fp[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,1,ngE));

    current_fp[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,1,ngJ));
    current_fp[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,1,ngJ));
    current_fp[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,1,ngJ));

    if (do_dive_cleaning)
    {
        F_fp[lev].reset  (new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,1, 0));
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,1,ngRho));
    }

    //
    // The Aux patch (i.e., the full solution)
    //
    if (lev == 0)
    {
        Bfield_aux[lev][0].reset(new MultiFab(*Bfield_fp[lev][0], amrex::make_alias, 0, 1));
        Bfield_aux[lev][1].reset(new MultiFab(*Bfield_fp[lev][1], amrex::make_alias, 0, 1));
        Bfield_aux[lev][2].reset(new MultiFab(*Bfield_fp[lev][2], amrex::make_alias, 0, 1));

        Efield_aux[lev][0].reset(new MultiFab(*Efield_fp[lev][0], amrex::make_alias, 0, 1));
        Efield_aux[lev][1].reset(new MultiFab(*Efield_fp[lev][1], amrex::make_alias, 0, 1));
        Efield_aux[lev][2].reset(new MultiFab(*Efield_fp[lev][2], amrex::make_alias, 0, 1));
    }
    else
    {
        Bfield_aux[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,1,ngE));
        Bfield_aux[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,1,ngE));
        Bfield_aux[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,1,ngE));

        Efield_aux[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,1,ngE));
        Efield_aux[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,1,ngE));
        Efield_aux[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,1,ngE));
    }

    //
    // The coarse patch
    //
    if (lev > 0)
    {
        BoxArray cba = ba;
        cba.coarsen(refRatio(lev-1));

        // Create the MultiFabs for B
        Bfield_cp[lev][0].reset( new MultiFab(amrex::convert(cba,Bx_nodal_flag),dm,1,ngE));
        Bfield_cp[lev][1].reset( new MultiFab(amrex::convert(cba,By_nodal_flag),dm,1,ngE));
        Bfield_cp[lev][2].reset( new MultiFab(amrex::convert(cba,Bz_nodal_flag),dm,1,ngE));

        // Create the MultiFabs for E
        Efield_cp[lev][0].reset( new MultiFab(amrex::convert(cba,Ex_nodal_flag),dm,1,ngE));
        Efield_cp[lev][1].reset( new MultiFab(amrex::convert(cba,Ey_nodal_flag),dm,1,ngE));
        Efield_cp[lev][2].reset( new MultiFab(amrex::convert(cba,Ez_nodal_flag),dm,1,ngE));

        // Create the MultiFabs for the current
        current_cp[lev][0].reset( new MultiFab(amrex::convert(cba,jx_nodal_flag),dm,1,ngJ));
        current_cp[lev][1].reset( new MultiFab(amrex::convert(cba,jy_nodal_flag),dm,1,ngJ));
        current_cp[lev][2].reset( new MultiFab(amrex::convert(cba,jz_nodal_flag),dm,1,ngJ));

        if (do_dive_cleaning)
        {
            F_cp[lev].reset  (new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,1, 0));
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,1,ngRho));
        }
    }

    if (load_balance_int > 0) {
        costs[lev].reset(new MultiFab(ba, dm, 1, 0));
    }
}

std::array<Real,3>
WarpX::CellSize (int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const Real* dx = gm.CellSize();
#if (BL_SPACEDIM == 3)
    return { dx[0], dx[1], dx[2] };
#elif (BL_SPACEDIM == 2)
    return { dx[0], 1.0, dx[1] };
#else
    static_assert(BL_SPACEDIM != 1, "1D is not supported");
#endif
}

std::array<Real,3>
WarpX::LowerCorner(const Box& bx, int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const RealBox grid_box{bx, gm.CellSize(), gm.ProbLo()};
    const Real* xyzmin = grid_box.lo();
#if (BL_SPACEDIM == 3)
    return { xyzmin[0], xyzmin[1], xyzmin[2] };
#elif (BL_SPACEDIM == 2)
    return { xyzmin[0], -1.e100, xyzmin[1] };
#endif
}

std::array<Real,3>
WarpX::UpperCorner(const Box& bx, int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const RealBox grid_box{bx, gm.CellSize(), gm.ProbLo()};
    const Real* xyzmax = grid_box.hi();
#if (BL_SPACEDIM == 3)
    return { xyzmax[0], xyzmax[1], xyzmax[2] };
#elif (BL_SPACEDIM == 2)
    return { xyzmax[0], 1.e100, xyzmax[1] };
#endif
}

void WarpX::zeroOutBoundary(amrex::MultiFab& input_data,
                            amrex::MultiFab& bndry_data,
                            const FabArray<BaseFab<int> >& mask) const {
    bndry_data.setVal(0.0, 1);
    for (MFIter mfi(input_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        WRPX_ZERO_OUT_BNDRY(bx.loVect(), bx.hiVect(),
                            input_data[mfi].dataPtr(),
                            bndry_data[mfi].dataPtr(),
                            mask[mfi].dataPtr());
    }
    bndry_data.FillBoundary();
}

void WarpX::sumFineToCrseNodal(const amrex::MultiFab& fine,
                               amrex::MultiFab& crse,
                               const amrex::Geometry& cgeom,
                               const amrex::IntVect& ratio) {
    const BoxArray& fine_BA = fine.boxArray();
    const DistributionMapping& fine_dm = fine.DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(ratio);

    MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, 1, 0);
    coarsened_fine_data.setVal(0.0);

    for (MFIter mfi(coarsened_fine_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Box& crse_box = coarsened_fine_data[mfi].box();
        const Box& fine_box = fine[mfi].box();
        WRPX_SUM_FINE_TO_CRSE_NODAL(bx.loVect(), bx.hiVect(), ratio.getVect(),
                                    coarsened_fine_data[mfi].dataPtr(), crse_box.loVect(), crse_box.hiVect(),
                                    fine[mfi].dataPtr(), fine_box.loVect(), fine_box.hiVect());
    }

    crse.copy(coarsened_fine_data, cgeom.periodicity(), FabArrayBase::ADD);
}

void
WarpX::fixRHSForSolve(Vector<std::unique_ptr<MultiFab> >& rhs,
                      const Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks) const {
    int num_levels = rhs.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        MultiFab& fine_rhs = *rhs[lev];
        const FabArray<BaseFab<int> >& mask = *masks[lev];
        const BoxArray& fine_ba = fine_rhs.boxArray();
        const DistributionMapping& fine_dm = fine_rhs.DistributionMap();
        MultiFab fine_bndry_data(fine_ba, fine_dm, 1, 1);
        zeroOutBoundary(fine_rhs, fine_bndry_data, mask);
    }
}

void WarpX::getLevelMasks(Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks,
                          const int nnodes) {
    int num_levels = grids.size();
    BL_ASSERT(num_levels == dmap.size());

    int covered = 0;
    int notcovered = 1;
    int physbnd = 1;
    int interior = 0;

    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();

        FabArray<BaseFab<int> > tmp_mask(nba, dmap[lev], 1, nnodes);
        tmp_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
                           covered, notcovered, physbnd, interior);
        masks[lev].reset(new FabArray<BaseFab<int> >(nba, dmap[lev], 1, 0));
        for (MFIter mfi(tmp_mask); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            WRPX_BUILD_MASK(bx.loVect(), bx.hiVect(),
                            tmp_mask[mfi].dataPtr(), (*masks[lev])[mfi].dataPtr(), &nnodes);
        }
    }
}


void WarpX::computePhi(const Vector<std::unique_ptr<MultiFab> >& rho,
                             Vector<std::unique_ptr<MultiFab> >& phi) const {


    int num_levels = rho.size();
    Vector<std::unique_ptr<MultiFab> > rhs(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        phi[lev]->setVal(0.0, 2);
        rhs[lev].reset(new MultiFab(rho[lev]->boxArray(), dmap[lev], 1, 0));
        MultiFab::Copy(*rhs[lev], *rho[lev], 0, 0, 1, 0);
        rhs[lev]->mult(-1.0/PhysConst::ep0, 0);
    }

    fixRHSForSolve(rhs, masks);

    bool nodal = true;
    bool have_rhcc = false;
    int  nc = 0;
    int Ncomp = 1;
    int stencil = ND_CROSS_STENCIL;
    int verbose = 0;
    Vector<int> mg_bc(2*BL_SPACEDIM, 1); // this means Dirichlet
    Real rel_tol = 1.0e-14;
    Real abs_tol = 1.0e-14;

    Vector<Geometry>            level_geom(1);
    Vector<BoxArray>            level_grids(1);
    Vector<DistributionMapping> level_dm(1);
    Vector<MultiFab*>           level_phi(1);
    Vector<MultiFab*>           level_rhs(1);

    for (int lev = 0; lev < num_levels; ++lev) {
        level_phi[0]   = phi[lev].get();
        level_rhs[0]   = rhs[lev].get();
        level_geom[0]  = geom[lev];
        level_grids[0] = grids[lev];
        level_dm[0]    = dmap[lev];

        MGT_Solver solver(level_geom, mg_bc.dataPtr(), level_grids,
                          level_dm, nodal,
                          stencil, have_rhcc, nc, Ncomp, verbose);

        solver.set_nodal_const_coefficients(1.0);

        solver.solve_nodal(level_phi, level_rhs, rel_tol, abs_tol);

        if (lev < num_levels-1) {

            NoOpPhysBC cphysbc, fphysbc;
#if BL_SPACEDIM == 3
            int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
            int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
#else
            int lo_bc[] = {INT_DIR, INT_DIR};
            int hi_bc[] = {INT_DIR, INT_DIR};
#endif
            Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
            NodeBilinear mapper;

            amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0, *phi[lev],
                                         0, 0, 1, geom[lev], geom[lev+1],
                                         cphysbc, fphysbc,
                                         IntVect(D_DECL(2, 2, 2)), &mapper, bcs);
        }
    }

    for (int lev = 0; lev < num_levels; ++lev) {
        const Geometry& gm = geom[lev];
        phi[lev]->FillBoundary(gm.periodicity());
    }
}

void WarpX::computeE(Vector<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                     const Vector<std::unique_ptr<MultiFab> >& phi) const {

    const int num_levels = E.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = GetInstance().Geom(lev);
        const Real* dx = gm.CellSize();
        for (MFIter mfi(*phi[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            WRPX_COMPUTE_E_NODAL(bx.loVect(), bx.hiVect(),
                                 (*phi[lev] )[mfi].dataPtr(),
                                 (*E[lev][0])[mfi].dataPtr(),
                                 (*E[lev][1])[mfi].dataPtr(),
#if BL_SPACEDIM == 3
                                 (*E[lev][2])[mfi].dataPtr(),
#endif
                                 dx);
        }

        E[lev][0]->FillBoundary(gm.periodicity());
        E[lev][1]->FillBoundary(gm.periodicity());
#if BL_SPACEDIM == 3
        E[lev][2]->FillBoundary(gm.periodicity());
#endif
    }
}

void
WarpX::ComputeDivB (MultiFab& divB, int dcomp,
                    const std::array<const MultiFab*, 3>& B,
                    const std::array<Real,3>& dx)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(divB, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        WRPX_COMPUTE_DIVB(bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_N_ANYD(divB[mfi],dcomp),
                           BL_TO_FORTRAN_ANYD((*B[0])[mfi]),
                           BL_TO_FORTRAN_ANYD((*B[1])[mfi]),
                           BL_TO_FORTRAN_ANYD((*B[2])[mfi]),
                           dx.data());
    }
}

void
WarpX::ComputeDivE (MultiFab& divE, int dcomp,
                    const std::array<const MultiFab*, 3>& E,
                    const std::array<Real,3>& dx)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(divE, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        WRPX_COMPUTE_DIVE(bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_N_ANYD(divE[mfi],dcomp),
                           BL_TO_FORTRAN_ANYD((*E[0])[mfi]),
                           BL_TO_FORTRAN_ANYD((*E[1])[mfi]),
                           BL_TO_FORTRAN_ANYD((*E[2])[mfi]),
                           dx.data());
    }
}
