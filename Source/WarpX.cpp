
#include <ParmParse.H>

#include <WarpX.H>

long WarpX::current_deposition_algo = 3;
long WarpX::charge_deposition_algo = 0;
long WarpX::field_gathering_algo = 1;
long WarpX::particle_pusher_algo = 0;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

#if (BL_SPACEDIM == 3)
IntVect WarpX::Bx_nodal_flag(1,0,0);
IntVect WarpX::By_nodal_flag(0,1,0);
IntVect WarpX::Bz_nodal_flag(0,0,1);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::Bx_nodal_flag(1,0);  // x is the first dimension to BoxLib
IntVect WarpX::By_nodal_flag(0,0);  // y is the missing dimension to 2D BoxLib
IntVect WarpX::Bz_nodal_flag(0,1);  // z is the second dimension to 2D BoxLib
#endif

#if (BL_SPACEDIM == 3)
IntVect WarpX::Ex_nodal_flag(1,0,0);
IntVect WarpX::Ey_nodal_flag(0,1,0);
IntVect WarpX::Ez_nodal_flag(0,0,1);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::Ex_nodal_flag(0,1);  // x is the first dimension to BoxLib
IntVect WarpX::Ey_nodal_flag(1,1);  // y is the missing dimension to 2D BoxLib
IntVect WarpX::Ez_nodal_flag(1,0);  // z is the second dimension to 2D BoxLib
#endif

#if (BL_SPACEDIM == 3)
IntVect WarpX::jx_nodal_flag(1,0,0);
IntVect WarpX::jy_nodal_flag(0,1,0);
IntVect WarpX::jz_nodal_flag(0,0,1);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::jx_nodal_flag(0,1);  // x is the first dimension to BoxLib
IntVect WarpX::jy_nodal_flag(1,1);  // y is the missing dimension to 2D BoxLib
IntVect WarpX::jz_nodal_flag(1,0);  // z is the second dimension to 2D BoxLib
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
    ReadParameters();

    if (max_level != 0) {
	BoxLib::Abort("WaprX: max_level must be zero");
    }

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    // Particle Container
    mypc = std::unique_ptr<MyParticleContainer> (new MyParticleContainer(this));

    current.resize(nlevs_max);
    Efield.resize(nlevs_max);
    Bfield.resize(nlevs_max);
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
    }

    {
	ParmParse pp("interpolation");
	pp.query("nox", nox);
	pp.query("noy", noy);
	pp.query("noz", noz);  
    }

    {
	ParmParse pp("algo");
	pp.query("current_deposition", current_deposition_algo);
	pp.query("charge_deposition", charge_deposition_algo);
	pp.query("field_gathering", field_gathering_algo);
	pp.query("particle_pusher", particle_pusher_algo);
    }
}

void
WarpX::MakeNewLevel (int lev, Real time,
		      const BoxArray& new_grids, const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // PICSAR assumes all fields are nodal plus one ghost cell.
    const IntVect& nodalflag = IntVect::TheUnitVector();
    const int ng = 1;

    current[lev].resize(3);
    Efield [lev].resize(3);
    Bfield [lev].resize(3);
    for (int i = 0; i < 3; ++i) {
	current[lev][i].reset(new MultiFab(grids[lev],1,ng,dmap[lev],Fab_allocate,nodalflag));
	Efield [lev][i].reset(new MultiFab(grids[lev],1,ng,dmap[lev],Fab_allocate,nodalflag));
	Bfield [lev][i].reset(new MultiFab(grids[lev],1,ng,dmap[lev],Fab_allocate,nodalflag));
    }
}

void 
WarpX::FillBoundary (MultiFab& mf, const Geometry& geom, const IntVect& nodalflag)
{
    const IndexType correct_typ(nodalflag);
    BoxArray ba = mf.boxArray();
    ba.convert(correct_typ);

    MultiFab tmpmf(ba, mf.nComp(), mf.nGrow(), mf.DistributionMap());

    const IndexType& mf_typ = mf.boxArray().ixType();

    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi) {
	FArrayBox& tmpfab = tmpmf[mfi];
	tmpfab.SetBoxType(mf_typ);
	tmpfab.copy(mf[mfi]);
	tmpfab.SetBoxType(correct_typ);
    }

    tmpmf.FillBoundary(geom.periodicity());

    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi) {
	FArrayBox& tmpfab = tmpmf[mfi];
	tmpfab.SetBoxType(mf_typ);
	mf[mfi].copy(tmpmf[mfi]);
	tmpfab.SetBoxType(correct_typ);
    }
}

void
WarpX::Copy(MultiFab& dstmf, int dcomp, int ncomp, const MultiFab& srcmf, int scomp)
{
    const IndexType& dst_typ = dstmf.boxArray().ixType();
    const IndexType& src_typ = srcmf.boxArray().ixType();

    for (MFIter mfi(dstmf); mfi.isValid(); ++mfi)
    {
	FArrayBox& dfab = dstmf[mfi];
	const FArrayBox& sfab = srcmf[mfi];
	dfab.SetBoxType(src_typ);
	dfab.copy(sfab, scomp, dcomp, ncomp);
	dfab.SetBoxType(dst_typ);
    }
}

