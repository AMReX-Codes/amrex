
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <WarpX.H>
#include <WarpX_f.H>
#include <BilinearFilter.H>
#include <NCIGodfreyFilter.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

using namespace amrex;

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

    if (restart_chkfile.empty())
    {
        ComputeDt();
        InitFromScratch();
    }
    else
    {
	InitFromCheckpoint();
        if (is_synchronized) {
            ComputeDt();
        }
	PostRestart();
    }

    ComputePMLFactors();

    if (WarpX::use_fdtd_nci_corr) {
        WarpX::InitNCICorrector();
    }

    if (WarpX::use_filter) {
        WarpX::InitFilter();
    }

    BuildBufferMasks();

    InitDiagnostics();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nGrids Summary:\n";
        printGridSummary(std::cout, 0, finestLevel());
    }

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge = new amrex::AmrMeshInSituBridge;
    insitu_bridge->setEnabled(insitu_int > 0 ? 1 : 0);
    insitu_bridge->setConfig(insitu_config);
    insitu_bridge->setPinMesh(insitu_pin_mesh);
    if (insitu_bridge->initialize())
    {
        amrex::ErrorStream()
            << "WarpX::InitData : Failed to initialize the in situ bridge."
            << std::endl;

        amrex::Abort();
    }
    insitu_bridge->setFrequency(1);
#endif

    if (restart_chkfile.empty())
    {
        if (plot_int > 0)
            WritePlotFile();

        if (check_int > 0)
            WriteCheckPointFile();

        if ((insitu_int > 0) && (insitu_start == 0))
            UpdateInSitu();
    }
}

void
WarpX::InitDiagnostics () {
    if (do_boosted_frame_diagnostic) {
        const Real* current_lo = geom[0].ProbLo();
        const Real* current_hi = geom[0].ProbHi();
        Real dt_boost = dt[0];

	// Find the positions of the lab-frame box that corresponds to the boosted-frame box at t=0
	Real zmin_lab = current_lo[moving_window_dir]/( (1.+beta_boost)*gamma_boost );
	Real zmax_lab = current_hi[moving_window_dir]/( (1.+beta_boost)*gamma_boost );

        myBFD.reset(new BoostedFrameDiagnostic(zmin_lab,
					       zmax_lab,
                                               moving_window_v, dt_snapshots_lab,
                                               num_snapshots_lab, gamma_boost,
                                               t_new[0], dt_boost,
                                               moving_window_dir, geom[0]));
    }
}

void
WarpX::InitFromScratch ()
{
    const Real time = 0.0;

    AmrCore::InitFromScratch(time);  // This will call MakeNewLevelFromScratch

    mypc->AllocData();
    mypc->InitData();

#ifdef USE_OPENBC_POISSON
    InitOpenbc();
#endif

    InitPML();

#ifdef WARPX_DO_ELECTROSTATIC
    if (do_electrostatic) {
        getLevelMasks(masks);

        // the plus one is to convert from num_cells to num_nodes
        getLevelMasks(gather_masks, n_buffer + 1);
    }
#endif // WARPX_DO_ELECTROSTATIC
}

void
WarpX::InitPML ()
{
    if (do_pml)
    {
        amrex::IntVect do_pml_Lo_corrected = do_pml_Lo;

#ifdef WARPX_DIM_RZ
        do_pml_Lo_corrected[0] = 0; // no PML at r=0, in cylindrical geometry
#endif
        pml[0].reset(new PML(boxArray(0), DistributionMap(0), &Geom(0), nullptr,
                             pml_ncell, pml_delta, 0,
#ifdef WARPX_USE_PSATD
                             dt[0], nox_fft, noy_fft, noz_fft, do_nodal,
#endif
                             do_dive_cleaning, do_moving_window,
                             pml_has_particles, do_pml_in_domain,
                             do_pml_Lo_corrected, do_pml_Hi));
        for (int lev = 1; lev <= finest_level; ++lev)
        {
            amrex::IntVect do_pml_Lo_MR = amrex::IntVect::TheUnitVector();
#ifdef WARPX_DIM_RZ
            //In cylindrical geometry, if the edge of the patch is at r=0, do not add PML
            if ((max_level > 0) && (fine_tag_lo[0]==0.)) {
                do_pml_Lo_MR[0] = 0;
            }
#endif
            pml[lev].reset(new PML(boxArray(lev), DistributionMap(lev),
                                   &Geom(lev), &Geom(lev-1),
                                   pml_ncell, pml_delta, refRatio(lev-1)[0],
#ifdef WARPX_USE_PSATD
                                   dt[lev], nox_fft, noy_fft, noz_fft, do_nodal,
#endif
                                   do_dive_cleaning, do_moving_window,
                                   pml_has_particles, do_pml_in_domain,
                                   do_pml_Lo_MR, amrex::IntVect::TheUnitVector()));
        }
    }
}

void
WarpX::ComputePMLFactors ()
{
    if (do_pml)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            pml[lev]->ComputePMLFactors(dt[lev]);
        }
    }
}

void
WarpX::InitNCICorrector ()
{
    if (WarpX::use_fdtd_nci_corr)
    {
        for (int lev = 0; lev <= max_level; ++lev)
        {
            const Geometry& gm = Geom(lev);
            const Real* dx = gm.CellSize();
            amrex::Real dz, cdtodz;
            if (AMREX_SPACEDIM == 3){
                dz = dx[2];
            }else{
                dz = dx[1];
            }
            cdtodz = PhysConst::c * dt[lev] / dz;

            // Initialize Godfrey filters
            // Same filter for fields Ex, Ey and Bz
            nci_godfrey_filter_exeybz[lev].reset( new NCIGodfreyFilter(godfrey_coeff_set::Ex_Ey_Bz, cdtodz, WarpX::l_lower_order_in_v) );
            // Same filter for fields Bx, By and Ez
            nci_godfrey_filter_bxbyez[lev].reset( new NCIGodfreyFilter(godfrey_coeff_set::Bx_By_Ez, cdtodz, WarpX::l_lower_order_in_v) );
            // Compute Godfrey filters stencils
            nci_godfrey_filter_exeybz[lev]->ComputeStencils();
            nci_godfrey_filter_bxbyez[lev]->ComputeStencils();
        }
    }
}

void
WarpX::InitFilter (){
    if (WarpX::use_filter){
        WarpX::bilinear_filter.npass_each_dir = WarpX::filter_npass_each_dir;
        WarpX::bilinear_filter.ComputeStencils();
    }
}

void
WarpX::PostRestart ()
{
#ifdef WARPX_USE_PSATD
    amrex::Abort("WarpX::PostRestart: TODO for PSATD");
#endif
    mypc->PostRestart();
}

#ifdef USE_OPENBC_POISSON
void
WarpX::InitOpenbc ()
{
#ifndef BL_USE_MPI
    static_assert(false, "must use MPI");
#endif

    static_assert(AMREX_SPACEDIM == 3, "Openbc is 3D only");
    BL_ASSERT(finestLevel() == 0);

    const int lev = 0;

    const Geometry& gm = Geom(lev);
    const Box& gbox = gm.Domain();
    int lohi[6];
    warpx_openbc_decompose(gbox.loVect(), gbox.hiVect(), lohi, lohi+3);

    int nprocs = ParallelDescriptor::NProcs();
    int myproc = ParallelDescriptor::MyProc();
    Vector<int> alllohi(6*nprocs,100000);

    MPI_Allgather(lohi, 6, MPI_INT, alllohi.data(), 6, MPI_INT, ParallelDescriptor::Communicator());

    BoxList bl{IndexType::TheNodeType()};
    for (int i = 0; i < nprocs; ++i)
    {
	bl.push_back(Box(IntVect(alllohi[6*i  ],alllohi[6*i+1],alllohi[6*i+2]),
			 IntVect(alllohi[6*i+3],alllohi[6*i+4],alllohi[6*i+5]),
			 IndexType::TheNodeType()));
    }
    BoxArray ba{bl};

    Vector<int> iprocmap(nprocs+1);
    std::iota(iprocmap.begin(), iprocmap.end(), 0);
    iprocmap.back() = myproc;

    DistributionMapping dm{iprocmap};

    MultiFab rho_openbc(ba, dm, 1, 0);
    MultiFab phi_openbc(ba, dm, 1, 0);

    bool local = true;
    const std::unique_ptr<MultiFab>& rho = mypc->GetChargeDensity(lev, local);

    rho_openbc.setVal(0.0);
    rho_openbc.copy(*rho, 0, 0, 1, rho->nGrow(), 0, gm.periodicity(), FabArrayBase::ADD);

    const Real* dx = gm.CellSize();

    warpx_openbc_potential(rho_openbc[myproc].dataPtr(), phi_openbc[myproc].dataPtr(), dx);

    BoxArray nba = boxArray(lev);
    nba.surroundingNodes();
    MultiFab phi(nba, DistributionMap(lev), 1, 0);
    phi.copy(phi_openbc, gm.periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.validbox();
	warpx_compute_E(bx.loVect(), bx.hiVect(),
			BL_TO_FORTRAN_3D(phi[mfi]),
			BL_TO_FORTRAN_3D((*Efield[lev][0])[mfi]),
			BL_TO_FORTRAN_3D((*Efield[lev][1])[mfi]),
			BL_TO_FORTRAN_3D((*Efield[lev][2])[mfi]),
			dx);
    }
}
#endif

void
WarpX::InitLevelData (int lev, Real time)
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

    if (F_fp[lev]) {
        F_fp[lev]->setVal(0.0);
    }

    if (rho_fp[lev]) {
        rho_fp[lev]->setVal(0.0);
    }

    if (F_cp[lev]) {
        F_cp[lev]->setVal(0.0);
    }

    if (rho_cp[lev]) {
        rho_cp[lev]->setVal(0.0);
    }

    if (costs[lev]) {
        costs[lev]->setVal(0.0);
    }
}

#ifdef WARPX_USE_PSATD_HYBRID

void
WarpX::InitLevelDataFFT (int lev, Real time)
{

    Efield_fp_fft[lev][0]->setVal(0.0);
    Efield_fp_fft[lev][1]->setVal(0.0);
    Efield_fp_fft[lev][2]->setVal(0.0);
    Bfield_fp_fft[lev][0]->setVal(0.0);
    Bfield_fp_fft[lev][1]->setVal(0.0);
    Bfield_fp_fft[lev][2]->setVal(0.0);
    current_fp_fft[lev][0]->setVal(0.0);
    current_fp_fft[lev][1]->setVal(0.0);
    current_fp_fft[lev][2]->setVal(0.0);
    rho_fp_fft[lev]->setVal(0.0);

    if (lev > 0)
    {
        Efield_cp_fft[lev][0]->setVal(0.0);
        Efield_cp_fft[lev][1]->setVal(0.0);
        Efield_cp_fft[lev][2]->setVal(0.0);
        Bfield_cp_fft[lev][0]->setVal(0.0);
        Bfield_cp_fft[lev][1]->setVal(0.0);
        Bfield_cp_fft[lev][2]->setVal(0.0);
        current_cp_fft[lev][0]->setVal(0.0);
        current_cp_fft[lev][1]->setVal(0.0);
        current_cp_fft[lev][2]->setVal(0.0);
        rho_cp_fft[lev]->setVal(0.0);
    }

}

#endif
