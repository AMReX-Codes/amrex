
#include <limits>

#include <MultiParticleContainer.H>
#include <WarpXParticleContainer.H>
#include <AMReX_AmrParGDB.H>
#include <WarpX_f.H>
#include <WarpX.H>

// Import low-level single-particle kernels
#include <GetAndSetPosition.H>
#include <UpdatePosition.H>

using namespace amrex;

int WarpXParticleContainer::do_not_push = 0;

WarpXParIter::WarpXParIter (ContainerType& pc, int level)
    : ParIter(pc, level, MFItInfo().SetDynamic(WarpX::do_dynamic_scheduling))
{
}

#if (AMREX_SPACEDIM == 2)
void
WarpXParIter::GetPosition (Cuda::ManagedDeviceVector<Real>& x, Cuda::ManagedDeviceVector<Real>& y, Cuda::ManagedDeviceVector<Real>& z) const
{
    amrex::ParIter<0,0,PIdx::nattribs>::GetPosition(x, z);
#ifdef WARPX_RZ
    const auto& attribs = GetAttribs();
    const auto& theta = attribs[PIdx::theta];
    y.resize(x.size());
    for (unsigned int i=0 ; i < x.size() ; i++) {
        // The x stored in the particles is actually the radius
        y[i] = x[i]*std::sin(theta[i]);
        x[i] = x[i]*std::cos(theta[i]);
    }
#else
    y.resize(x.size(), std::numeric_limits<Real>::quiet_NaN());
#endif
}

void
WarpXParIter::SetPosition (const Cuda::ManagedDeviceVector<Real>& x, const Cuda::ManagedDeviceVector<Real>& y, const Cuda::ManagedDeviceVector<Real>& z)
{
#ifdef WARPX_RZ
    auto& attribs = GetAttribs();
    auto& theta = attribs[PIdx::theta];
    Cuda::DeviceVector<Real> r(x.size());
    for (unsigned int i=0 ; i < x.size() ; i++) {
        theta[i] = std::atan2(y[i], x[i]);
        r[i] = std::sqrt(x[i]*x[i] + y[i]*y[i]);
    }
    amrex::ParIter<0,0,PIdx::nattribs>::SetPosition(r, z);
#else
    amrex::ParIter<0,0,PIdx::nattribs>::SetPosition(x, z);
#endif
}
#endif

WarpXParticleContainer::WarpXParticleContainer (AmrCore* amr_core, int ispecies)
    : ParticleContainer<0,0,PIdx::nattribs>(amr_core->GetParGDB())
    , species_id(ispecies)
{
    for (unsigned int i = PIdx::Ex; i <= PIdx::Bz; ++i) {
        communicate_real_comp[i] = false; // Don't need to communicate E and B.
    }
    SetParticleSize();
    ReadParameters();

    // build up the map of string names to particle component numbers
    particle_comps["w"]  = PIdx::w;
    particle_comps["ux"] = PIdx::ux;
    particle_comps["uy"] = PIdx::uy;
    particle_comps["uz"] = PIdx::uz;
    particle_comps["Ex"] = PIdx::Ex;
    particle_comps["Ey"] = PIdx::Ey;
    particle_comps["Ez"] = PIdx::Ez;
    particle_comps["Bx"] = PIdx::Bx;
    particle_comps["By"] = PIdx::By;
    particle_comps["Bz"] = PIdx::Bz;
#ifdef WARPX_RZ
    particle_comps["theta"] = PIdx::theta;
#endif

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        particle_comps["xold"]  = PIdx::nattribs;
        particle_comps["yold"]  = PIdx::nattribs+1;
        particle_comps["zold"]  = PIdx::nattribs+2;
        particle_comps["uxold"] = PIdx::nattribs+3;
        particle_comps["uyold"] = PIdx::nattribs+4;
        particle_comps["uzold"] = PIdx::nattribs+5;

    }

    // Initialize temporary local arrays for charge/current deposition
    int num_threads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    #pragma omp single
    num_threads = omp_get_num_threads();
    #endif
    local_rho.resize(num_threads);
    local_jx.resize(num_threads);
    local_jy.resize(num_threads);
    local_jz.resize(num_threads);
    m_xp.resize(num_threads);
    m_yp.resize(num_threads);
    m_zp.resize(num_threads);
    m_giv.resize(num_threads);
}

void
WarpXParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
	ParmParse pp("particles");

#ifdef AMREX_USE_GPU
        do_tiling = false; // By default, tiling is off on GPU
#else
        do_tiling = true;
#endif
        pp.query("do_tiling",  do_tiling);
        pp.query("do_not_push", do_not_push);

	initialized = true;
    }
}

void
WarpXParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because grids have not
    // been built when constructor was called.
    reserveData();
    resizeData();
}

void
WarpXParticleContainer::AddOneParticle (int lev, int grid, int tile,
                                        Real x, Real y, Real z,
                                        std::array<Real,PIdx::nattribs>& attribs)
{
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];
    AddOneParticle(particle_tile, x, y, z, attribs);
}

void
WarpXParticleContainer::AddOneParticle (ParticleTileType& particle_tile,
                                        Real x, Real y, Real z,
                                        std::array<Real,PIdx::nattribs>& attribs)
{
    ParticleType p;
    p.id()  = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();
#if (AMREX_SPACEDIM == 3)
    p.pos(0) = x;
    p.pos(1) = y;
    p.pos(2) = z;
#elif (AMREX_SPACEDIM == 2)
#ifdef WARPX_RZ
    attribs[PIdx::theta] = std::atan2(y, x);
    x = std::sqrt(x*x + y*y);
#endif
    p.pos(0) = x;
    p.pos(1) = z;
#endif

    particle_tile.push_back(p);
    particle_tile.push_back_real(attribs);
}

void
WarpXParticleContainer::AddNParticles (int lev,
                                       int n, const Real* x, const Real* y, const Real* z,
				       const Real* vx, const Real* vy, const Real* vz,
				       int nattr, const Real* attr, int uniqueparticles, int id)
{
    BL_ASSERT(nattr == 1);
    const Real* weight = attr;

    int ibegin, iend;
    if (uniqueparticles) {
	ibegin = 0;
	iend = n;
    } else {
	int myproc = ParallelDescriptor::MyProc();
	int nprocs = ParallelDescriptor::NProcs();
	int navg = n/nprocs;
	int nleft = n - navg * nprocs;
	if (myproc < nleft) {
	    ibegin = myproc*(navg+1);
	    iend = ibegin + navg+1;
	} else {
	    ibegin = myproc*navg + nleft;
	    iend = ibegin + navg;
	}
    }

    //  Add to grid 0 and tile 0
    // Redistribute() will move them to proper places.
    std::pair<int,int> key {0,0};
    auto& particle_tile = GetParticles(lev)[key];

    std::size_t np = iend-ibegin;

#ifdef WARPX_RZ
    Vector<Real> theta(np);
#endif

    for (int i = ibegin; i < iend; ++i)
    {
        ParticleType p;
	if (id==-1)
	{
	    p.id() = ParticleType::NextID();
	} else {
	    p.id() = id;
	}
        p.cpu() = ParallelDescriptor::MyProc();
#if (AMREX_SPACEDIM == 3)
        p.pos(0) = x[i];
        p.pos(1) = y[i];
        p.pos(2) = z[i];
#elif (AMREX_SPACEDIM == 2)
#ifdef WARPX_RZ
        theta[i-ibegin] = std::atan2(y[i], x[i]);
        p.pos(0) = std::sqrt(x[i]*x[i] + y[i]*y[i]);
#else
        p.pos(0) = x[i];
#endif
        p.pos(1) = z[i];
#endif

        if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
        {
            auto& ptile = DefineAndReturnParticleTile(0, 0, 0);
            ptile.push_back_real(particle_comps["xold"], x[i]);
            ptile.push_back_real(particle_comps["yold"], y[i]);
            ptile.push_back_real(particle_comps["zold"], z[i]);
        }

        particle_tile.push_back(p);
    }

    if (np > 0)
    {
        particle_tile.push_back_real(PIdx::w , weight + ibegin, weight + iend);
        particle_tile.push_back_real(PIdx::ux,     vx + ibegin,     vx + iend);
        particle_tile.push_back_real(PIdx::uy,     vy + ibegin,     vy + iend);
        particle_tile.push_back_real(PIdx::uz,     vz + ibegin,     vz + iend);

        if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
        {
            auto& ptile = DefineAndReturnParticleTile(0, 0, 0);
            ptile.push_back_real(particle_comps["uxold"], vx + ibegin, vx + iend);
            ptile.push_back_real(particle_comps["uyold"], vy + ibegin, vy + iend);
            ptile.push_back_real(particle_comps["uzold"], vz + ibegin, vz + iend);
        }

        for (int comp = PIdx::uz+1; comp < PIdx::nattribs; ++comp)
        {
#ifdef WARPX_RZ
            if (comp == PIdx::theta) {
                particle_tile.push_back_real(comp, theta.front(), theta.back());
            }
            else {
                particle_tile.push_back_real(comp, np, 0.0);
            }
#else
            particle_tile.push_back_real(comp, np, 0.0);
#endif
        }
    }

    Redistribute();
}


void
WarpXParticleContainer::DepositCurrent(WarpXParIter& pti,
                                       RealVector& wp, RealVector& uxp,
                                       RealVector& uyp, RealVector& uzp,
                                       MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                       MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                       const long np_current, const long np,
                                       int thread_num, int lev, Real dt )
{
  Real *jx_ptr, *jy_ptr, *jz_ptr;
  const std::array<Real,3>& xyzmin_tile = WarpX::LowerCorner(pti.tilebox(), lev);
  const std::array<Real,3>& dx = WarpX::CellSize(lev);
  const std::array<Real,3>& cdx = WarpX::CellSize(std::max(lev-1,0));
  const std::array<Real, 3>& xyzmin = xyzmin_tile;
  const long lvect = 8;

  BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr_cd);
  BL_PROFILE_VAR_NS("PPC::Evolve::Accumulate", blp_accumulate);

  Box tbx = convert(pti.tilebox(), WarpX::jx_nodal_flag);
  Box tby = convert(pti.tilebox(), WarpX::jy_nodal_flag);
  Box tbz = convert(pti.tilebox(), WarpX::jz_nodal_flag);

  // WarpX assumes the same number of guard cells for Jx, Jy, Jz
  long ngJ = jx.nGrow();

  bool j_is_nodal = jx.is_nodal() and jy.is_nodal() and jz.is_nodal();

  // Deposit charge for particles that are not in the current buffers
  if (np_current > 0)
  {
#ifdef AMREX_USE_GPU
      jx_ptr = jx[pti].dataPtr();
      jy_ptr = jy[pti].dataPtr();
      jz_ptr = jz[pti].dataPtr();

      auto jxntot = jx[pti].length();
      auto jyntot = jy[pti].length();
      auto jzntot = jz[pti].length();
#else
      tbx.grow(ngJ);
      tby.grow(ngJ);
      tbz.grow(ngJ);

      local_jx[thread_num].resize(tbx);
      local_jy[thread_num].resize(tby);
      local_jz[thread_num].resize(tbz);

      jx_ptr = local_jx[thread_num].dataPtr();
      jy_ptr = local_jy[thread_num].dataPtr();
      jz_ptr = local_jz[thread_num].dataPtr();

      local_jx[thread_num].setVal(0.0);
      local_jy[thread_num].setVal(0.0);
      local_jz[thread_num].setVal(0.0);

      auto jxntot = local_jx[thread_num].length();
      auto jyntot = local_jy[thread_num].length();
      auto jzntot = local_jz[thread_num].length();
#endif

      BL_PROFILE_VAR_START(blp_pxr_cd);
      if (j_is_nodal) {
          const Real* p_wp = wp.dataPtr();
          const Real* p_gaminv = m_giv[thread_num].dataPtr();
          const Real* p_uxp = uxp.dataPtr();
          const Real* p_uyp = uyp.dataPtr();
          const Real* p_uzp = uzp.dataPtr();
          AsyncArray<Real> wptmp_aa(np_current);
          Real* const wptmp = wptmp_aa.data();
          const Box& tile_box = pti.tilebox();
#if (AMREX_SPACEDIM == 3)
          const long nx = tile_box.length(0);
          const long ny = tile_box.length(1);
          const long nz = tile_box.length(2);
#else
          const long nx = tile_box.length(0);
          const long ny = 0;
          const long nz = tile_box.length(1);
#endif
          amrex::ParallelFor (np_current, [=] AMREX_GPU_DEVICE (long ip) {
                  wptmp[ip] = p_wp[ip] * p_gaminv[ip] * p_uxp[ip];
          });
          warpx_charge_deposition(jx_ptr, &np_current,
                                  m_xp[thread_num].dataPtr(),
                                  m_yp[thread_num].dataPtr(),
                                  m_zp[thread_num].dataPtr(),
                                  wptmp,
                                  &this->charge,
                                  &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                  &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                  &ngJ, &ngJ, &ngJ,
                                  &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                  &lvect, &WarpX::current_deposition_algo);
          amrex::ParallelFor (np_current, [=] AMREX_GPU_DEVICE (long ip) {
                  wptmp[ip] = p_wp[ip] * p_gaminv[ip] * p_uyp[ip];
          });
          warpx_charge_deposition(jy_ptr, &np_current,
                                  m_xp[thread_num].dataPtr(),
                                  m_yp[thread_num].dataPtr(),
                                  m_zp[thread_num].dataPtr(),
                                  wptmp,
                                  &this->charge,
                                  &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                  &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                  &ngJ, &ngJ, &ngJ,
                                  &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                  &lvect, &WarpX::current_deposition_algo);
          amrex::ParallelFor (np_current, [=] AMREX_GPU_DEVICE (long ip) {
                  wptmp[ip] = p_wp[ip] * p_gaminv[ip] * p_uzp[ip];
          });
          warpx_charge_deposition(jz_ptr, &np_current,
                                  m_xp[thread_num].dataPtr(),
                                  m_yp[thread_num].dataPtr(),
                                  m_zp[thread_num].dataPtr(),
                                  wptmp,
                                  &this->charge,
                                  &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                  &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                  &ngJ, &ngJ, &ngJ,
                                  &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                  &lvect, &WarpX::current_deposition_algo);
      } else {
          warpx_current_deposition(
                               jx_ptr, &ngJ, jxntot.getVect(),
                               jy_ptr, &ngJ, jyntot.getVect(),
                               jz_ptr, &ngJ, jzntot.getVect(),
                               &np_current,
                               m_xp[thread_num].dataPtr(),
                               m_yp[thread_num].dataPtr(),
                               m_zp[thread_num].dataPtr(),
                               uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
                               m_giv[thread_num].dataPtr(),
                               wp.dataPtr(), &this->charge,
                               &xyzmin[0], &xyzmin[1], &xyzmin[2],
                               &dt, &dx[0], &dx[1], &dx[2],
                               &WarpX::nox,&WarpX::noy,&WarpX::noz,
                               &lvect,&WarpX::current_deposition_algo);

#ifdef WARPX_RZ
         warpx_current_deposition_rz_volume_scaling(
                                  jx_ptr, &ngJ, jxntot.getVect(),
                                  jy_ptr, &ngJ, jyntot.getVect(),
                                  jz_ptr, &ngJ, jzntot.getVect(),
                                  &xyzmin[0], &dx[0]);
#endif
      }

      BL_PROFILE_VAR_STOP(blp_pxr_cd);

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_START(blp_accumulate);

      jx[pti].atomicAdd(local_jx[thread_num], tbx, tbx, 0, 0, 1);
      jy[pti].atomicAdd(local_jy[thread_num], tby, tby, 0, 0, 1);
      jz[pti].atomicAdd(local_jz[thread_num], tbz, tbz, 0, 0, 1);

      BL_PROFILE_VAR_STOP(blp_accumulate);
#endif
  }

  // Deposit charge for particles that are in the current buffers
  if (np_current < np)
  {
      const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
      const Box& ctilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
      const std::array<Real,3>& cxyzmin_tile = WarpX::LowerCorner(ctilebox, lev-1);

#ifdef AMREX_USE_GPU
      jx_ptr = (*cjx)[pti].dataPtr();
      jy_ptr = (*cjy)[pti].dataPtr();
      jz_ptr = (*cjz)[pti].dataPtr();

      auto jxntot = jx[pti].length();
      auto jyntot = jy[pti].length();
      auto jzntot = jz[pti].length();
#else

      tbx = amrex::convert(ctilebox, WarpX::jx_nodal_flag);
      tby = amrex::convert(ctilebox, WarpX::jy_nodal_flag);
      tbz = amrex::convert(ctilebox, WarpX::jz_nodal_flag);
      tbx.grow(ngJ);
      tby.grow(ngJ);
      tbz.grow(ngJ);

      local_jx[thread_num].resize(tbx);
      local_jy[thread_num].resize(tby);
      local_jz[thread_num].resize(tbz);

      jx_ptr = local_jx[thread_num].dataPtr();
      jy_ptr = local_jy[thread_num].dataPtr();
      jz_ptr = local_jz[thread_num].dataPtr();

      local_jx[thread_num].setVal(0.0);
      local_jy[thread_num].setVal(0.0);
      local_jz[thread_num].setVal(0.0);

      auto jxntot = local_jx[thread_num].length();
      auto jyntot = local_jy[thread_num].length();
      auto jzntot = local_jz[thread_num].length();
#endif

      long ncrse = np - np_current;
      BL_PROFILE_VAR_START(blp_pxr_cd);
      if (j_is_nodal) {
          const Real* p_wp = wp.dataPtr() + np_current;
          const Real* p_gaminv = m_giv[thread_num].dataPtr() + np_current;
          const Real* p_uxp = uxp.dataPtr() + np_current;
          const Real* p_uyp = uyp.dataPtr() + np_current;
          const Real* p_uzp = uzp.dataPtr() + np_current;
          AsyncArray<Real> wptmp_aa(ncrse);
          Real* const wptmp = wptmp_aa.data();
          const Box& tile_box = pti.tilebox();
#if (AMREX_SPACEDIM == 3)
          const long nx = tile_box.length(0);
          const long ny = tile_box.length(1);
          const long nz = tile_box.length(2);
#else
          const long nx = tile_box.length(0);
          const long ny = 0;
          const long nz = tile_box.length(1);
#endif
          amrex::ParallelFor (ncrse, [=] AMREX_GPU_DEVICE (long ip) {
                  wptmp[ip] = p_wp[ip] * p_gaminv[ip] * p_uxp[ip];
          });
          warpx_charge_deposition(jx_ptr, &ncrse,
                                  m_xp[thread_num].dataPtr() +np_current,
                                  m_yp[thread_num].dataPtr() +np_current,
                                  m_zp[thread_num].dataPtr() +np_current,
                                  wptmp,
                                  &this->charge,
                                  &cxyzmin_tile[0], &cxyzmin_tile[1], &cxyzmin_tile[2],
                                  &cdx[0], &cdx[1], &cdx[2], &nx, &ny, &nz,
                                  &ngJ, &ngJ, &ngJ,
                                  &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                  &lvect, &WarpX::current_deposition_algo);
          amrex::ParallelFor (ncrse, [=] AMREX_GPU_DEVICE (long ip) {
                  wptmp[ip] = p_wp[ip] * p_gaminv[ip] * p_uyp[ip];
          });
          warpx_charge_deposition(jy_ptr, &ncrse,
                                  m_xp[thread_num].dataPtr() +np_current,
                                  m_yp[thread_num].dataPtr() +np_current,
                                  m_zp[thread_num].dataPtr() +np_current,
                                  wptmp,
                                  &this->charge,
                                  &cxyzmin_tile[0], &cxyzmin_tile[1], &cxyzmin_tile[2],
                                  &cdx[0], &cdx[1], &cdx[2], &nx, &ny, &nz,
                                  &ngJ, &ngJ, &ngJ,
                                  &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                  &lvect, &WarpX::current_deposition_algo);
          amrex::ParallelFor (ncrse, [=] AMREX_GPU_DEVICE (long ip) {
                  wptmp[ip] = p_wp[ip] * p_gaminv[ip] * p_uzp[ip];
          });
          warpx_charge_deposition(jz_ptr, &ncrse,
                                  m_xp[thread_num].dataPtr() +np_current,
                                  m_yp[thread_num].dataPtr() +np_current,
                                  m_zp[thread_num].dataPtr() +np_current,
                                  wptmp,
                                  &this->charge,
                                  &cxyzmin_tile[0], &cxyzmin_tile[1], &cxyzmin_tile[2],
                                  &cdx[0], &cdx[1], &cdx[2], &nx, &ny, &nz,
                                  &ngJ, &ngJ, &ngJ,
                                  &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                  &lvect, &WarpX::current_deposition_algo);
      } else {
          warpx_current_deposition(
                               jx_ptr, &ngJ, jxntot.getVect(),
                               jy_ptr, &ngJ, jyntot.getVect(),
                               jz_ptr, &ngJ, jzntot.getVect(),
                               &ncrse,
                               m_xp[thread_num].dataPtr() +np_current,
                               m_yp[thread_num].dataPtr() +np_current,
                               m_zp[thread_num].dataPtr() +np_current,
                               uxp.dataPtr()+np_current,
                               uyp.dataPtr()+np_current,
                               uzp.dataPtr()+np_current,
                               m_giv[thread_num].dataPtr()+np_current,
                               wp.dataPtr()+np_current, &this->charge,
                               &cxyzmin_tile[0], &cxyzmin_tile[1], &cxyzmin_tile[2],
                               &dt, &cdx[0], &cdx[1], &cdx[2],
                               &WarpX::nox,&WarpX::noy,&WarpX::noz,
                               &lvect,&WarpX::current_deposition_algo);
#ifdef WARPX_RZ
         warpx_current_deposition_rz_volume_scaling(
                                  jx_ptr, &ngJ, jxntot.getVect(),
                                  jy_ptr, &ngJ, jyntot.getVect(),
                                  jz_ptr, &ngJ, jzntot.getVect(),
                                  &xyzmin[0], &dx[0]);
#endif
      }

      BL_PROFILE_VAR_STOP(blp_pxr_cd);

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_START(blp_accumulate);

      (*cjx)[pti].atomicAdd(local_jx[thread_num], tbx, tbx, 0, 0, 1);
      (*cjy)[pti].atomicAdd(local_jy[thread_num], tby, tby, 0, 0, 1);
      (*cjz)[pti].atomicAdd(local_jz[thread_num], tbz, tbz, 0, 0, 1);

      BL_PROFILE_VAR_STOP(blp_accumulate);
#endif
    }
};


void
WarpXParticleContainer::DepositCharge ( WarpXParIter& pti, RealVector& wp,
                                        MultiFab* rhomf, MultiFab* crhomf, int icomp,
                                        const long np_current,
                                        const long np, int thread_num, int lev )
{

  BL_PROFILE_VAR_NS("PICSAR::ChargeDeposition", blp_pxr_chd);
  BL_PROFILE_VAR_NS("PPC::Evolve::Accumulate", blp_accumulate);

  const std::array<Real,3>& xyzmin_tile = WarpX::LowerCorner(pti.tilebox(), lev);
  const long lvect = 8;

  long ngRho = rhomf->nGrow();
  Real* data_ptr;
  Box tile_box = convert(pti.tilebox(), IntVect::TheUnitVector());

  const std::array<Real,3>& dx = WarpX::CellSize(lev);
  const std::array<Real,3>& cdx = WarpX::CellSize(std::max(lev-1,0));

  // Deposit charge for particles that are not in the current buffers
  if (np_current > 0)
  {
      const std::array<Real, 3>& xyzmin = xyzmin_tile;

#ifdef AMREX_USE_GPU
      data_ptr = (*rhomf)[pti].dataPtr(icomp);
      auto rholen = (*rhomf)[pti].length();
#else
      tile_box.grow(ngRho);
      local_rho[thread_num].resize(tile_box);

      data_ptr = local_rho[thread_num].dataPtr();
      auto rholen = local_rho[thread_num].length();

      local_rho[thread_num].setVal(0.0);
#endif

#if (AMREX_SPACEDIM == 3)
      const long nx = rholen[0]-1-2*ngRho;
      const long ny = rholen[1]-1-2*ngRho;
      const long nz = rholen[2]-1-2*ngRho;
#else
      const long nx = rholen[0]-1-2*ngRho;
      const long ny = 0;
      const long nz = rholen[1]-1-2*ngRho;
#endif
      BL_PROFILE_VAR_START(blp_pxr_chd);
      warpx_charge_deposition(data_ptr, &np_current,
                              m_xp[thread_num].dataPtr(),
                              m_yp[thread_num].dataPtr(),
                              m_zp[thread_num].dataPtr(),
                              wp.dataPtr(),
                              &this->charge,
                              &xyzmin[0], &xyzmin[1], &xyzmin[2],
                              &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                              &ngRho, &ngRho, &ngRho,
                              &WarpX::nox,&WarpX::noy,&WarpX::noz,
                              &lvect, &WarpX::charge_deposition_algo);
#ifdef WARPX_RZ
      warpx_charge_deposition_rz_volume_scaling(
                               data_ptr, &ngRho, rholen.getVect(),
                               &xyzmin[0], &dx[0]);
#endif
      BL_PROFILE_VAR_STOP(blp_pxr_chd);

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_START(blp_accumulate);

      (*rhomf)[pti].atomicAdd(local_rho[thread_num], tile_box, tile_box, 0, icomp, 1);

      BL_PROFILE_VAR_STOP(blp_accumulate);
#endif
  }

  // Deposit charge for particles that are in the current buffers
  if (np_current < np)
  {
      const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
      const Box& ctilebox = amrex::coarsen(pti.tilebox(), ref_ratio);
      const std::array<Real,3>& cxyzmin_tile = WarpX::LowerCorner(ctilebox, lev-1);

#ifdef AMREX_USE_GPU
      data_ptr = (*crhomf)[pti].dataPtr();
      auto rholen = (*crhomf)[pti].length();
#else
      tile_box = amrex::convert(ctilebox, IntVect::TheUnitVector());
      tile_box.grow(ngRho);
      local_rho[thread_num].resize(tile_box);

      data_ptr = local_rho[thread_num].dataPtr();
      auto rholen = local_rho[thread_num].length();

      local_rho[thread_num].setVal(0.0);
#endif

#if (AMREX_SPACEDIM == 3)
      const long nx = rholen[0]-1-2*ngRho;
      const long ny = rholen[1]-1-2*ngRho;
      const long nz = rholen[2]-1-2*ngRho;
#else
      const long nx = rholen[0]-1-2*ngRho;
      const long ny = 0;
      const long nz = rholen[1]-1-2*ngRho;
#endif

      long ncrse = np - np_current;
      BL_PROFILE_VAR_START(blp_pxr_chd);
      warpx_charge_deposition(data_ptr, &ncrse,
                              m_xp[thread_num].dataPtr() + np_current,
                              m_yp[thread_num].dataPtr() + np_current,
                              m_zp[thread_num].dataPtr() + np_current,
                              wp.dataPtr() + np_current,
                              &this->charge,
                              &cxyzmin_tile[0], &cxyzmin_tile[1], &cxyzmin_tile[2],
                              &cdx[0], &cdx[1], &cdx[2], &nx, &ny, &nz,
                              &ngRho, &ngRho, &ngRho,
                              &WarpX::nox,&WarpX::noy,&WarpX::noz,
                              &lvect, &WarpX::charge_deposition_algo);
#ifdef WARPX_RZ
      warpx_charge_deposition_rz_volume_scaling(
                               data_ptr, &ngRho, rholen.getVect(),
                               &cxyzmin_tile[0], &cdx[0]);
#endif
      BL_PROFILE_VAR_STOP(blp_pxr_chd);

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_START(blp_accumulate);

      (*crhomf)[pti].atomicAdd(local_rho[thread_num], tile_box, tile_box, 0, icomp, 1);

      BL_PROFILE_VAR_STOP(blp_accumulate);
#endif
    }
};

void
WarpXParticleContainer::DepositCharge (Vector<std::unique_ptr<MultiFab> >& rho, bool local)
{

    int num_levels = rho.size();
    int finest_level = num_levels - 1;

    // each level deposits it's own particles
    const int ng = rho[0]->nGrow();
    for (int lev = 0; lev < num_levels; ++lev) {

        rho[lev]->setVal(0.0, ng);

        const auto& gm = m_gdb->Geom(lev);
        const auto& ba = m_gdb->ParticleBoxArray(lev);

        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();
        BoxArray nba = ba;
        nba.surroundingNodes();

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Box& box = nba[pti];

            auto& wp = pti.GetAttribs(PIdx::w);
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();

            FArrayBox& rhofab = (*rho[lev])[pti];

            WRPX_DEPOSIT_CIC(particles.dataPtr(), nstride, np,
                             wp.dataPtr(), &this->charge,
                             rhofab.dataPtr(), box.loVect(), box.hiVect(),
                             plo, dx, &ng);
        }

        if (!local) rho[lev]->SumBoundary(gm.periodicity());
    }

    // now we average down fine to crse
    std::unique_ptr<MultiFab> crse;
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        const BoxArray& fine_BA = rho[lev+1]->boxArray();
        const DistributionMapping& fine_dm = rho[lev+1]->DistributionMap();
        BoxArray coarsened_fine_BA = fine_BA;
        coarsened_fine_BA.coarsen(m_gdb->refRatio(lev));

        MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, 1, 0);
        coarsened_fine_data.setVal(0.0);

        IntVect ratio(AMREX_D_DECL(2, 2, 2));  // FIXME

        for (MFIter mfi(coarsened_fine_data); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const Box& crse_box = coarsened_fine_data[mfi].box();
            const Box& fine_box = (*rho[lev+1])[mfi].box();
            WRPX_SUM_FINE_TO_CRSE_NODAL(bx.loVect(), bx.hiVect(), ratio.getVect(),
                                        coarsened_fine_data[mfi].dataPtr(), crse_box.loVect(), crse_box.hiVect(),
                                        (*rho[lev+1])[mfi].dataPtr(), fine_box.loVect(), fine_box.hiVect());
        }

        rho[lev]->copy(coarsened_fine_data, m_gdb->Geom(lev).periodicity(), FabArrayBase::ADD);
    }
}

std::unique_ptr<MultiFab>
WarpXParticleContainer::GetChargeDensity (int lev, bool local)
{
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;
    nba.surroundingNodes();

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    const int ng = WarpX::nox;

    auto rho = std::unique_ptr<MultiFab>(new MultiFab(nba,dm,1,ng));
    rho->setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
        Cuda::ManagedDeviceVector<Real> xp, yp, zp;
#ifdef _OPENMP
        FArrayBox rho_loc;
#endif

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& wp = pti.GetAttribs(PIdx::w);

            const long np  = pti.numParticles();

            pti.GetPosition(xp, yp, zp);

            // Data on the grid
            Real* data_ptr;
            FArrayBox& rhofab = (*rho)[pti];
#ifdef _OPENMP
            const std::array<Real,3>& xyzmin_tile = WarpX::LowerCorner(pti.tilebox(), lev);
            Box tile_box = convert(pti.tilebox(), IntVect::TheUnitVector());
            const std::array<Real, 3>& xyzmin = xyzmin_tile;
            tile_box.grow(ng);
            rho_loc.resize(tile_box);
            rho_loc = 0.0;
            data_ptr = rho_loc.dataPtr();
            auto rholen = rho_loc.length();
#else
            const Box& box = pti.validbox();
            const std::array<Real,3>& xyzmin_grid = WarpX::LowerCorner(box, lev);
            const std::array<Real, 3>& xyzmin = xyzmin_grid;
            data_ptr = rhofab.dataPtr();
            auto rholen = rhofab.length();
#endif

#if (AMREX_SPACEDIM == 3)
            const long nx = rholen[0]-1-2*ng;
            const long ny = rholen[1]-1-2*ng;
            const long nz = rholen[2]-1-2*ng;
#else
            const long nx = rholen[0]-1-2*ng;
            const long ny = 0;
            const long nz = rholen[1]-1-2*ng;
#endif

            long nxg = ng;
            long nyg = ng;
            long nzg = ng;
            long lvect = 8;

            warpx_charge_deposition(data_ptr,
                                    &np,
                                    xp.dataPtr(),
                                    yp.dataPtr(),
                                    zp.dataPtr(), wp.dataPtr(),
                                    &this->charge, &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                    &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                    &nxg, &nyg, &nzg, &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                    &lvect, &WarpX::charge_deposition_algo);
#ifdef WARPX_RZ
            long ngRho = WarpX::nox;
            warpx_charge_deposition_rz_volume_scaling(
                                     data_ptr, &ngRho, rholen.getVect(),
                                     &xyzmin[0], &dx[0]);
#endif

#ifdef _OPENMP
            rhofab.atomicAdd(rho_loc);
        }
#endif
    }

    if (!local) rho->SumBoundary(gm.periodicity());

    return rho;
}

Real WarpXParticleContainer::sumParticleCharge(bool local) {

    amrex::Real total_charge = 0.0;

    for (int lev = 0; lev < finestLevel(); ++lev)
    {

#ifdef _OPENMP
#pragma omp parallel reduction(+:total_charge)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& wp = pti.GetAttribs(PIdx::w);
            for (unsigned long i = 0; i < wp.size(); i++) {
                total_charge += wp[i];
            }
        }
    }

    if (!local) ParallelDescriptor::ReduceRealSum(total_charge);
    total_charge *= this->charge;
    return total_charge;
}

std::array<Real, 3> WarpXParticleContainer::meanParticleVelocity(bool local) {

    amrex::Real vx_total = 0.0;
    amrex::Real vy_total = 0.0;
    amrex::Real vz_total = 0.0;

    long np_total = 0;

    amrex::Real inv_clight_sq = 1.0/PhysConst::c/PhysConst::c;

    for (int lev = 0; lev <= finestLevel(); ++lev) {

#ifdef _OPENMP
#pragma omp parallel reduction(+:vx_total, vy_total, vz_total, np_total)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ux = pti.GetAttribs(PIdx::ux);
            auto& uy = pti.GetAttribs(PIdx::uy);
            auto& uz = pti.GetAttribs(PIdx::uz);

            np_total += pti.numParticles();

            for (unsigned long i = 0; i < ux.size(); i++) {
                Real usq = (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i])*inv_clight_sq;
                Real gaminv = 1.0/std::sqrt(1.0 + usq);
                vx_total += ux[i]*gaminv;
                vy_total += uy[i]*gaminv;
                vz_total += uz[i]*gaminv;
            }
        }
    }

    if (!local) {
        ParallelDescriptor::ReduceRealSum(vx_total);
        ParallelDescriptor::ReduceRealSum(vy_total);
        ParallelDescriptor::ReduceRealSum(vz_total);
        ParallelDescriptor::ReduceLongSum(np_total);
    }

    std::array<Real, 3> mean_v;
    if (np_total > 0) {
        mean_v[0] = vx_total / np_total;
        mean_v[1] = vy_total / np_total;
        mean_v[2] = vz_total / np_total;
    }

    return mean_v;
}

Real WarpXParticleContainer::maxParticleVelocity(bool local) {

    amrex::Real max_v = 0.0;

    for (int lev = 0; lev <= finestLevel(); ++lev)
    {

#ifdef _OPENMP
#pragma omp parallel reduction(max:max_v)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ux = pti.GetAttribs(PIdx::ux);
            auto& uy = pti.GetAttribs(PIdx::uy);
            auto& uz = pti.GetAttribs(PIdx::uz);
            for (unsigned long i = 0; i < ux.size(); i++) {
                max_v = std::max(max_v, sqrt(ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]));
            }
        }
    }

    if (!local) ParallelDescriptor::ReduceRealMax(max_v);
    return max_v;
}

void
WarpXParticleContainer::PushXES (Real dt)
{
    BL_PROFILE("WPC::PushXES()");

    int num_levels = finestLevel() + 1;

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];

            WRPX_PUSH_LEAPFROG_POSITIONS(particles.dataPtr(), nstride, np,
                                         uxp.dataPtr(), uyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                                         uzp.dataPtr(),
#endif
                                         &dt,
                                         prob_domain.lo(), prob_domain.hi());
        }
    }
}

void
WarpXParticleContainer::PushX (Real dt)
{
    for (int lev = 0; lev <= finestLevel(); ++lev) {
        PushX(lev, dt);
    }
}

void
WarpXParticleContainer::PushX (int lev, Real dt)
{
    BL_PROFILE("WPC::PushX()");

    if (do_not_push) return;

    MultiFab* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            Real wt = amrex::second();

            //
            // Particle Push
            //
            // Extract pointers to particle position and momenta, for this particle tile
            // - positions are stored as an array of struct, in `ParticleType`
            ParticleType * AMREX_RESTRICT pstructs = &(pti.GetArrayOfStructs()[0]);
            // - momenta are stored as a struct of array, in `attribs`
            auto& attribs = pti.GetAttribs();
            Real* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            Real* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            Real* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
#ifdef WARPX_RZ
            Real* AMREX_RESTRICT theta = attribs[PIdx::theta].dataPtr();
#endif
            // Loop over the particles and update their position
            amrex::ParallelFor( pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    ParticleType& p = pstructs[i]; // Particle object that gets updated
                    Real x, y, z; // Temporary variables
#ifndef WARPX_RZ
                    GetPosition( x, y, z, p ); // Initialize x, y, z
                    UpdatePosition( x, y, z, ux[i], uy[i], uz[i], dt);
                    SetPosition( p, x, y, z ); // Update the object p
#else
                    // For WARPX_RZ, the particles are still pushed in 3D Cartesian
                    GetCartesianPositionFromCylindrical( x, y, z, p, theta[i] );
                    UpdatePosition( x, y, z, ux[i], uy[i], uz[i], dt);
                    SetCylindricalPositionFromCartesian( p, theta[i], x, y, z );
#endif
                }
            );

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (amrex::second() - wt) / tbx.d_numPts();
                FArrayBox* costfab = cost->fabPtr(pti);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, work_box,
                {
                    costfab->plus(wt, work_box);
                });
            }
        }
    }
}

// This function is called in Redistribute, just after locate
void
WarpXParticleContainer::particlePostLocate(ParticleType& p,
                                           const ParticleLocData& pld,
                                           const int lev)
{
    // Tag particle if goes to higher level.
    // It will be split later in the loop
    if (pld.m_lev == lev+1
        and p.m_idata.id != NoSplitParticleID
        and p.m_idata.id >= 0)
    {
        p.m_idata.id = DoSplitParticleID;
    }
    // For the moment, do not do anything if particles goes
    // to lower level.
    if (pld.m_lev == lev-1){
    }
}
