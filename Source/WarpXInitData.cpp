
#include <numeric>

#include <AMReX_ParallelDescriptor.H>

#include <WarpX.H>
#include <WarpX_f.H>

using namespace amrex;

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

    if (restart_chkfile.empty())
    {
	InitFromScratch();
        ComputeDt();
    }
    else
    {
	InitFromCheckpoint();
	PostRestart();
        if (is_synchronized) {
            ComputeDt();
        }
    }

    ComputePMLFactors();

    BuildBufferMasks();

    InitDiagnostics();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nGrids Summary:\n";
        printGridSummary(std::cout, 0, finestLevel());
    }

    if (restart_chkfile.empty())
    {
	if (plot_int > 0) {
            WritePlotFile();
	}
	if (check_int > 0) {
	    WriteCheckPointFile();
	}
    }
}

void
WarpX::InitDiagnostics () {
    if (do_boosted_frame_diagnostic) {
        const Real* current_lo = geom[0].ProbLo();
        const Real* current_hi = geom[0].ProbHi();
        Real dt_boost = dt[0];
        
        myBFD.reset(new BoostedFrameDiagnostic(current_lo[moving_window_dir], 
                                               current_hi[moving_window_dir],
                                               moving_window_v, dt_snapshots_lab,
                                               num_snapshots_lab, gamma_boost, dt_boost, 
                                               moving_window_dir));
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

    if (do_electrostatic) {
        getLevelMasks(masks);
        
        // the plus one is to convert from num_cells to num_nodes
        getLevelMasks(gather_masks, n_buffer + 1);
    }
}

void
WarpX::InitPML ()
{
    if (do_pml)
    {
        pml[0].reset(new PML(boxArray(0), DistributionMap(0), &Geom(0), nullptr,
                             pml_ncell, pml_delta, 0, do_dive_cleaning, do_moving_window));
        for (int lev = 1; lev <= finest_level; ++lev)
        {
            pml[lev].reset(new PML(boxArray(lev), DistributionMap(lev),
                                   &Geom(lev), &Geom(lev-1),
                                   pml_ncell, pml_delta, refRatio(lev-1)[0], do_dive_cleaning,
                                   do_moving_window));
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
            pml[lev]->ComputePMLFactors(dt[lev],pml_type);
        }
    }
}

void
WarpX::PostRestart ()
{
    mypc->PostRestart();
}

#ifdef USE_OPENBC_POISSON
void
WarpX::InitOpenbc ()
{
#ifndef BL_USE_MPI
    static_assert(false, "must use MPI");
#endif

    static_assert(BL_SPACEDIM == 3, "Openbc is 3D only");
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
