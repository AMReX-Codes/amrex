
#include <limits>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>

#include <AMReX_ParmParse.H>

#include <WarpX.H>
#include <WarpXConst.H>

using namespace amrex;

long WarpX::current_deposition_algo = 3;
long WarpX::charge_deposition_algo = 0;
long WarpX::field_gathering_algo = 1;
long WarpX::particle_pusher_algo = 0;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

bool WarpX::use_laser         = false;

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

    if (max_level != 0) {
	amrex::Abort("WaprX: max_level must be zero");
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
    mypc = std::unique_ptr<MultiParticleContainer> (new MultiParticleContainer(this));

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
	pp.query("regrid_int", regrid_int);

        // PML
        if (Geometry::isAllPeriodic()) {
            do_pml = 0;  // no PML for all periodic boundaries
        } else {
            pp.query("do_pml", do_pml);
            pp.query("pml_ncell", pml_ncell);
        }

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

	pp.query("use_laser", use_laser);

        pp.query("plot_raw_fields", plot_raw_fields);
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
    InitLevelData(time);
}

void
WarpX::ClearLevel (int lev)
{
    for (int i = 0; i < 3; ++i) {
	current[lev][i].reset();
	Efield [lev][i].reset();
	Bfield [lev][i].reset();
    }
}

void
WarpX::AllocLevelData (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    const int ng = WarpX::nox;  // need to update this

    // Create the MultiFabs for B
    Bfield[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,1,ng));
    Bfield[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,1,ng));
    Bfield[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,1,ng));

    // Create the MultiFabs for E
    Efield[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,1,ng));
    Efield[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,1,ng));
    Efield[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,1,ng));

    // Create the MultiFabs for the current
    current[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,1,ng));
    current[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,1,ng));
    current[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,1,ng));
}

void
WarpX::InitPML ()
{
    if (!do_pml) return;

    const Geometry& gm0 = Geom(0);
    const Box& domainbox = gm0.Domain();
    Box grownbox = domainbox;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if (!Geometry::isPeriodic(i)) {
            grownbox.grow(i,pml_ncell);
        }
    }

    int block_size = maxGridSize(0);
    while (block_size < pml_ncell) {
        block_size += blockingFactor(0);
    }

    Array<IntVect> shift;
    {
        int len[3] = {0,0,0};
        int jmp[3] = {1,1,1};

        for (int i = 0; i < BL_SPACEDIM; ++i) {
            if (!Geometry::isPeriodic(i)) {
                len[i] = jmp[i] = domainbox.length(i);
            }
        }

        for (int i = -len[0]; i <= len[0]; i += jmp[0]) {
        for (int j = -len[1]; j <= len[1]; j += jmp[1]) {
        for (int k = -len[2]; k <= len[2]; k += jmp[2]) {
            if (i != 0 || j != 0 || k!= 0) {
                shift.push_back(IntVect(D_DECL(i,j,k)));
            }
        }
        }
        }
    }

    BoxList bl;
    for (const IntVect& iv : shift) {
        Box bbx = domainbox;
        bbx.shift(iv);
        bbx &= grownbox;
        BoxList bltmp(bbx);
        bltmp.maxSize(block_size);
        bl.catenate(bltmp);
    }
    pml_ba.define(bl);
    
    // xxxxx for now let's just use a simple distributionmapping
    pml_dm.RoundRobinProcessorMap(pml_ba.size(), ParallelDescriptor::NProcs());

    const int ng = Efield[0][0]->nGrow();
    pml_E[0].reset(new MultiFab(amrex::convert(pml_ba,Ex_nodal_flag),pml_dm,2,ng));
    pml_E[1].reset(new MultiFab(amrex::convert(pml_ba,Ey_nodal_flag),pml_dm,2,ng));
    pml_E[2].reset(new MultiFab(amrex::convert(pml_ba,Ez_nodal_flag),pml_dm,2,ng));
    pml_B[0].reset(new MultiFab(amrex::convert(pml_ba,Bx_nodal_flag),pml_dm,2,ng));
    pml_B[1].reset(new MultiFab(amrex::convert(pml_ba,By_nodal_flag),pml_dm,2,ng));
    pml_B[2].reset(new MultiFab(amrex::convert(pml_ba,Bz_nodal_flag),pml_dm,2,ng));

    pml_E[0]->setVal(0.0);
    pml_E[1]->setVal(0.0);
    pml_E[2]->setVal(0.0);
    pml_B[0]->setVal(0.0);
    pml_B[1]->setVal(0.0);
    pml_B[2]->setVal(0.0);

#if (BL_SPACEDIM == 3)
    {
        const Real* dx = gm0.CellSize();
        const int* dlo = domainbox.loVect();
        const int* dhi = domainbox.hiVect();
        const int* glo = grownbox.loVect();
        const int* ghi = grownbox.hiVect();
        const IntVect& sz = grownbox.size();
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            pml_sigma     [idim].resize(sz[idim]+1, 0.0);
            pml_sigma_star[idim].resize(sz[idim]  , 0.0);
            pml_sigma     [idim].m_lo = glo[idim];
            pml_sigma     [idim].m_hi = ghi[idim]+1;
            pml_sigma_star[idim].m_lo = glo[idim];
            pml_sigma_star[idim].m_hi = ghi[idim];

            const Real fac = 4.0*PhysConst::c/(dx[idim]*static_cast<Real>(pml_ncell*pml_ncell));

            for (int ind = glo[idim]; ind < dlo[idim]; ++ind) {
                Real offset = static_cast<Real>(dlo[idim] - ind);
                pml_sigma[idim][ind-glo[idim]] = fac*(offset*offset);
            }
            for (int ind = dhi[idim]+2; ind < ghi[idim]+2; ++ind) {
                Real offset = static_cast<Real>(ind - (dhi[idim]+1));
                pml_sigma[idim][ind-glo[idim]] = fac*(offset*offset);
            }

            for (int icc = glo[idim]; icc < dlo[idim]; ++icc) {
                Real offset = static_cast<Real>(dlo[idim] - icc) - 0.5;
                pml_sigma_star[idim][icc-glo[idim]] = fac*(offset*offset);
            }
            for (int icc = dhi[idim]+1; icc < ghi[idim]+1; ++icc) {
                Real offset = static_cast<Real>(icc - dhi[idim]) - 0.5;
                pml_sigma_star[idim][icc-glo[idim]] = fac*(offset*offset);                
            }
        }
    }
#else
    amrex::Abort("InitPML: 2d not supported yet");
#endif

    const int lev = 0;
    ComputePMLFactors(lev, dt[lev]);
}

void
WarpX::shiftMF(MultiFab& mf, const Geometry& geom, int num_shift, int dir)
{
    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    const int nc = mf.nComp();
    const int ng = std::max(mf.nGrow(), std::abs(num_shift));
    MultiFab tmpmf(ba, dm, nc, ng);
    MultiFab::Copy(tmpmf, mf, 0, 0, nc, ng);
    tmpmf.FillBoundary(geom.periodicity());

    // Zero out the region that the window moved into
    const IndexType& typ = ba.ixType();
    const Box& domainBox = geom.Domain();
    Box adjBox;
    if (num_shift > 0) {
        adjBox = adjCellHi(domainBox, dir, ng);
    } else {
        adjBox = adjCellLo(domainBox, dir, ng);
    }
    adjBox = amrex::convert(adjBox, typ);

    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        if (idim == dir and typ.nodeCentered(dir)) {
            if (num_shift > 0) {
                adjBox.growLo(idim, -1);
            } else {
                adjBox.growHi(idim, -1);
            }
        } else if (idim != dir) {
            adjBox.growLo(idim, ng);
            adjBox.growHi(idim, ng);
        }
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi )
    {
        FArrayBox& srcfab = tmpmf[mfi];

        Box outbox = mfi.fabbox();
        outbox &= adjBox;
        if (outbox.ok()) {  // outbox is the region that the window moved into
            srcfab.setVal(0.0, outbox, 0, nc);
        }

        FArrayBox& dstfab = mf[mfi];
        dstfab.setVal(0.0);

        Box dstBox = dstfab.box();

        if (num_shift > 0) {
            dstBox.growHi(dir, -num_shift);
        } else {
            dstBox.growLo(dir,  num_shift);
        }

        dstfab.copy(srcfab, amrex::shift(dstBox,dir,num_shift), 0, dstBox, 0, nc);
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


#if (BL_SPACEDIM == 3)

Box
WarpX::getIndexBox(const RealBox& real_box) const
{
  BL_ASSERT(max_level == 0);

  IntVect slice_lo, slice_hi;

  D_TERM(slice_lo[0]=std::floor((real_box.lo(0) - geom[0].ProbLo(0))/geom[0].CellSize(0));,
	 slice_lo[1]=std::floor((real_box.lo(1) - geom[0].ProbLo(1))/geom[0].CellSize(1));,
	 slice_lo[2]=std::floor((real_box.lo(2) - geom[0].ProbLo(2))/geom[0].CellSize(2)););

  D_TERM(slice_hi[0]=std::floor((real_box.hi(0) - geom[0].ProbLo(0))/geom[0].CellSize(0));,
	 slice_hi[1]=std::floor((real_box.hi(1) - geom[0].ProbLo(1))/geom[0].CellSize(1));,
	 slice_hi[2]=std::floor((real_box.hi(2) - geom[0].ProbLo(2))/geom[0].CellSize(2)););

  return Box(slice_lo, slice_hi) & geom[0].Domain();
}

void
WarpX::fillSlice(Real z_coord) const
{
  BL_ASSERT(max_level == 0);

  // Get our slice and convert to index space
  RealBox real_slice = geom[0].ProbDomain();
  real_slice.setLo(2, z_coord);
  real_slice.setHi(2, z_coord);
  Box slice_box = getIndexBox(real_slice);

  // define the multifab that stores slice
  BoxArray ba = Efield[0][0]->boxArray();
  const IndexType correct_typ(IntVect::TheZeroVector());
  ba.convert(correct_typ);
  const DistributionMapping& dm = dmap[0];
  std::vector< std::pair<int, Box> > isects;
  ba.intersections(slice_box, isects, false, 0);
  Array<Box> boxes;
  Array<int> procs;
  for (int i = 0; i < isects.size(); ++i) {
    procs.push_back(dm[isects[i].first]);
    boxes.push_back(isects[i].second);
  }
  procs.push_back(ParallelDescriptor::MyProc());
  BoxArray slice_ba(&boxes[0], boxes.size());
  DistributionMapping slice_dmap(procs);
  MultiFab slice(slice_ba, slice_dmap, 6, 0);

  const MultiFab* mfs[6];
  mfs[0] = Efield[0][0].get();
  mfs[1] = Efield[0][1].get();
  mfs[2] = Efield[0][2].get();
  mfs[3] = Bfield[0][0].get();
  mfs[4] = Bfield[0][1].get();
  mfs[5] = Bfield[0][2].get();

  IntVect flags[6];
  flags[0] = WarpX::Ex_nodal_flag;
  flags[1] = WarpX::Ey_nodal_flag;
  flags[2] = WarpX::Ez_nodal_flag;
  flags[3] = WarpX::Bx_nodal_flag;
  flags[4] = WarpX::By_nodal_flag;
  flags[5] = WarpX::Bz_nodal_flag;

  // Fill the slice with sampled data
  const Real* dx      = geom[0].CellSize();
  const Geometry& g   = geom[0];
  for (MFIter mfi(slice); mfi.isValid(); ++mfi) {
    int slice_gid = mfi.index();
    Box grid = slice_ba[slice_gid];
    int xlo = grid.smallEnd(0), ylo = grid.smallEnd(1), zlo = grid.smallEnd(2);
    int xhi = grid.bigEnd(0), yhi = grid.bigEnd(1), zhi = grid.bigEnd(2);

    IntVect iv = grid.smallEnd();
    ba.intersections(Box(iv, iv), isects, true, 0);
    BL_ASSERT(!isects.empty());
    int full_gid = isects[0].first;

    for (int k = zlo; k <= zhi; k++) {
      for (int j = ylo; j <= yhi; j++) {
	for (int i = xlo; i <= xhi; i++) {
	  for (int comp = 0; comp < 6; comp++) {
	    Real x = g.ProbLo(0) + i*dx[0];
	    Real y = g.ProbLo(1) + j*dx[1];
	    Real z = z_coord;

	    D_TERM(iv[0]=std::floor((x - g.ProbLo(0) + 0.5*g.CellSize(0)*flags[comp][0])/g.CellSize(0));,
		   iv[1]=std::floor((y - g.ProbLo(1) + 0.5*g.CellSize(1)*flags[comp][1])/g.CellSize(1));,
		   iv[2]=std::floor((z - g.ProbLo(2) + 0.5*g.CellSize(2)*flags[comp][2])/g.CellSize(2)););

	    slice[slice_gid](IntVect(i, j, k), comp) = (*mfs[comp])[full_gid](iv);
	  }
	}
      }
    }
  }
}

void WarpX::sampleAtPoints(const Array<Real>& x,
			   const Array<Real>& y,
			   const Array<Real>& z,
			   Array<Array<Real> >& result) const {
  BL_ASSERT((x.size() == y.size()) and (y.size() == z.size()));
  BL_ASSERT(max_level == 0);

  const MultiFab* mfs[6];
  mfs[0] = Efield[0][0].get();
  mfs[1] = Efield[0][1].get();
  mfs[2] = Efield[0][2].get();
  mfs[3] = Bfield[0][0].get();
  mfs[4] = Bfield[0][1].get();
  mfs[5] = Bfield[0][2].get();

  IntVect flags[6];
  flags[0] = WarpX::Ex_nodal_flag;
  flags[1] = WarpX::Ey_nodal_flag;
  flags[2] = WarpX::Ez_nodal_flag;
  flags[3] = WarpX::Bx_nodal_flag;
  flags[4] = WarpX::By_nodal_flag;
  flags[5] = WarpX::Bz_nodal_flag;

  const unsigned npoints = x.size();
  const int ncomps = 6;
  result.resize(ncomps);
  for (int i = 0; i < ncomps; i++) {
    result[i].resize(npoints, 0);
  }

  BoxArray ba = Efield[0][0]->boxArray();
  const IndexType correct_typ(IntVect::TheZeroVector());
  ba.convert(correct_typ);
  std::vector< std::pair<int, Box> > isects;

  IntVect iv;
  const Geometry& g   = geom[0];
  for (int i = 0; i < npoints; ++i) {
    for (int comp = 0; comp < ncomps; comp++) {
      D_TERM(iv[0]=std::floor((x[i] - g.ProbLo(0) + 0.5*g.CellSize(0)*flags[comp][0])/g.CellSize(0));,
	     iv[1]=std::floor((y[i] - g.ProbLo(1) + 0.5*g.CellSize(1)*flags[comp][1])/g.CellSize(1));,
	     iv[2]=std::floor((z[i] - g.ProbLo(2) + 0.5*g.CellSize(2)*flags[comp][2])/g.CellSize(2)););

      ba.intersections(Box(iv, iv), isects, true, 0);
      const int grid = isects[0].first;
      const int who = dmap[0][grid];
      if (who == ParallelDescriptor::MyProc()) {
	result[comp][i] = (*mfs[comp])[grid](iv);
      }
    }
  }

  for (int i = 0; i < ncomps; i++) {
    ParallelDescriptor::ReduceRealSum(result[i].dataPtr(), result[i].size());
  }
}
#endif
