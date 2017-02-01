
#include <numeric>

#include <ParallelDescriptor.H>

#include <WarpX.H>
#include <WarpX_f.H>

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

    if (restart_chkfile.empty())
    {
	InitFromScratch();

	if (plot_int > 0) {
	    WritePlotFile();
	}
	if (check_int > 0) {
	    WriteCheckPointFile();
	}
    }
    else
    {
	InitFromCheckpoint();
	PostRestart();
    }
}

void
WarpX::InitFromScratch ()
{
    BL_ASSERT(max_level == 0);

    const Real time = 0.0;

    // define coarse level BoxArray and DistributionMap
    {
	finest_level = 0;

	t_new[0] = time;
	t_old[0] = time - 1.e200;
    
	const BoxArray& ba = MakeBaseGrids();
	DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	MakeNewLevel(0, ba, dm);

	InitLevelData(0);
    }

    // if max_level > 0, define fine levels

    mypc->AllocData();
    mypc->InitData();

#ifdef USE_OPENBC_POISSON
    InitOpenbc();
#endif
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
    Array<int> alllohi(6*nprocs,100000);

    MPI_Allgather(lohi, 6, MPI_INT, alllohi.data(), 6, MPI_INT, ParallelDescriptor::Communicator());
    
    BoxList bl{IndexType::TheNodeType()};
    for (int i = 0; i < nprocs; ++i)
    {
	bl.push_back(Box(IntVect(alllohi[6*i  ],alllohi[6*i+1],alllohi[6*i+2]),
			 IntVect(alllohi[6*i+3],alllohi[6*i+4],alllohi[6*i+5]),
			 IndexType::TheNodeType()));
    }
    BoxArray ba{bl};

    Array<int> iprocmap(nprocs+1);
    std::iota(iprocmap.begin(), iprocmap.end(), 0);
    iprocmap.back() = myproc;

    DistributionMapping dm{iprocmap};

    MultiFab rho_openbc(ba, 1, 0, dm);
    MultiFab phi_openbc(ba, 1, 0, dm);

    bool local = true;
    const std::unique_ptr<MultiFab>& rho = mypc->GetChargeDensity(lev, local);

    rho_openbc.setVal(0.0);
    rho_openbc.copy(*rho, 0, 0, 1, rho->nGrow(), 0, gm.periodicity(), FabArrayBase::ADD);

    const Real* dx = gm.CellSize();
    
    warpx_openbc_potential(rho_openbc[myproc].dataPtr(), phi_openbc[myproc].dataPtr(), dx);

    BoxArray nba = boxArray(lev);
    nba.surroundingNodes();
    MultiFab phi(nba, 1, 0);
    phi.copy(phi_openbc, gm.periodicity());

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
