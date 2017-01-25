
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <ParticleContainer.H>
#include <ParticleIterator.H>


//
// xxxxx need to make this work in 2D!
//

namespace
{
    Array<Real> CrossProduct (const Array<Real>& a, const Array<Real>& b)
    {
	return { a[1]*b[2]-a[2]*b[1],  a[2]*b[0]-a[0]*b[2],  a[0]*b[1]-a[1]*b[0] };
    }
}

LaserParticleContainer::LaserParticleContainer (AmrCore* amr_core, int ispecies)
    : WarpXParticleContainer(amr_core, ispecies)
{
    charge = 1.0;
    mass = std::numeric_limits<Real>::max();

    if (WarpX::use_laser)
    {
	ParmParse pp("laser");

	std::string laser_type_s;
	pp.get("profile", laser_type_s);
	std::transform(laser_type_s.begin(), laser_type_s.end(), laser_type_s.begin(), ::tolower);
	if (laser_type_s == "gaussian") {
	    profile = laser_t::Gaussian;
	} else {
	    BoxLib::Abort("Unknown laser type");
	}

	pp.getarr("position", position);
	pp.getarr("direction", nvec);
	pp.getarr("polarization", p_X);
	pp.query("pusher_algo", pusher_algo);

	pp.get("e_max", e_max);
	pp.get("profile_waist", profile_waist);
	pp.get("profile_duration", profile_duration);
	pp.get("wavelength", wavelength);

	// Plane normal
	Real s = 1.0/std::sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
	nvec = { nvec[0]*s, nvec[1]*s, nvec[2]*s };

	// The first polarization vector
	s = 1.0/std::sqrt(p_X[0]*p_X[0] + p_X[1]*p_X[1] + p_X[2]*p_X[2]);
	p_X = { p_X[0]*s, p_X[1]*s, p_X[2]*s };

	Real dp = std::inner_product(nvec.begin(), nvec.end(), p_X.begin(), 0.0);
	if (std::abs(dp) > 1.e-14) {
	    BoxLib::Abort("Laser plane vector is not perpendicular to the main polarization vector");
	}

	p_Y = CrossProduct(nvec, p_X);   // The second polarization vector

#if BL_SPACEDIM == 3
	u_X = p_X;
	u_Y = p_Y;
#else
	u_X = CrossProduct({0., 1., 0.}, nvec);
	u_Y = {0., 1., 0.};
#endif
    }
}

void
LaserParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because GDB was not
    // ready in constructor.
    m_particles.resize(GDB().finestLevel()+1);
}

void
LaserParticleContainer::InitData ()
{
    m_particles.resize(GDB().finestLevel()+1);

    const int lev = 0;

    const Geometry& geom = GDB().Geom(lev);
    const Real* dx  = geom.CellSize();
    const RealBox& prob_domain = geom.ProbDomain();

    // spacing of laser particles in the laser plane
    const Real eps = dx[0]*1.e-50;
    const Real S_X = std::min(std::min(dx[0]/(std::abs(p_X[0])+eps),
				       dx[1]/(std::abs(p_X[1])+eps)),
			               dx[2]/(std::abs(p_X[2])+eps));
    const Real S_Y = std::min(std::min(dx[0]/(std::abs(p_Y[0])+eps),
				       dx[1]/(std::abs(p_Y[1])+eps)),
                                       dx[2]/(std::abs(p_Y[2])+eps));

    const Real particle_weight = ComputeWeight(S_X, S_Y);

    // Give integer coordinates in the laser plane, return the real coordinates in the "lab" frame
    auto Transform = [&](int i, int j) -> Array<Real>
    {
	return { position[0] + (S_X*i)*p_X[0] + (S_Y*j)*p_Y[0],
		 position[1] + (S_X*i)*p_X[1] + (S_Y*j)*p_Y[1],
		 position[2] + (S_X*i)*p_X[2] + (S_Y*j)*p_Y[2] };
    };

    // Given the "lab" frame coordinates, return the real coordinates in the laser plane coordinates
    auto InverseTransform = [&](const Array<Real>& pos) -> Array<Real>
    {
	return { p_X[0]*pos[0] + p_X[1]*pos[1] + p_X[2]*pos[2],
		 p_Y[0]*pos[0] + p_Y[1]*pos[1] + p_Y[2]*pos[2],
		 nvec[0]*pos[0] + nvec[1]*pos[1] + nvec[2]*pos[2] };
    };

    Array<int> plane_lo(2, std::numeric_limits<int>::max());
    Array<int> plane_hi(2, std::numeric_limits<int>::min());
    {
	auto compute_min_max = [&](Real x, Real y, Real z)
	{
	    const Array<Real> pos_plane = InverseTransform({x, y, z});
	    int i = pos_plane[0]/S_X;
	    int j = pos_plane[1]/S_Y;
	    plane_lo[0] = std::min(plane_lo[0], i);
	    plane_lo[1] = std::min(plane_lo[1], j);
	    plane_hi[0] = std::max(plane_hi[0], i+1);
	    plane_hi[1] = std::max(plane_hi[1], j+1);
	};

	const Real* prob_lo = prob_domain.lo();
	const Real* prob_hi = prob_domain.hi();
	compute_min_max(prob_lo[0], prob_lo[1], prob_lo[2]);
	compute_min_max(prob_hi[0], prob_lo[1], prob_lo[2]);
	compute_min_max(prob_lo[0], prob_hi[1], prob_lo[2]);
	compute_min_max(prob_hi[0], prob_hi[1], prob_lo[2]);
	compute_min_max(prob_lo[0], prob_lo[1], prob_hi[2]);
	compute_min_max(prob_hi[0], prob_lo[1], prob_hi[2]);
	compute_min_max(prob_lo[0], prob_hi[1], prob_hi[2]);
	compute_min_max(prob_hi[0], prob_hi[1], prob_hi[2]);
    }

    const int nprocs = ParallelDescriptor::NProcs();
    const int myproc = ParallelDescriptor::MyProc();

    const Box plane_box {IntVect(D_DECL(plane_lo[0],plane_lo[1],0)),
                         IntVect(D_DECL(plane_hi[0],plane_hi[1],0))};
    BoxArray plane_ba {plane_box};
    {
	IntVect chunk(plane_box.size());
	const int min_size = 8;
	while (plane_ba.size() < nprocs && chunk[0] > min_size && chunk[1] > min_size)
	{
	    for (int j = 1; j >= 0 ; j--)
	    {
		chunk[j] /= 2;

		if (plane_ba.size() < nprocs) {
		    plane_ba.maxSize(chunk);
		}
	    }
	}
    }

    Array<Real> particle_x, particle_y, particle_z, particle_w;

    const DistributionMapping plane_dm {plane_ba, nprocs};
    const Array<int>& procmap = plane_dm.ProcessorMap();
    for (int i = 0, n = plane_ba.size(); i < n; ++i)
    {
	if (procmap[i] == myproc)
	{
	    const Box& bx = plane_ba[i];
	    for (IntVect cell = bx.smallEnd(); cell <= bx.bigEnd(); bx.next(cell))
	    {
		const Array<Real>& pos = Transform(cell[0], cell[1]);
		if (prob_domain.contains(pos.data()))
		{
		    for (int k = 0; k<2; ++k) {
			particle_x.push_back(pos[0]);
			particle_y.push_back(pos[1]);
			particle_z.push_back(pos[2]);
		    }
		    particle_w.push_back( particle_weight);
		    particle_w.push_back(-particle_weight);
		}
	    }
	}
    }
    const int np = particle_z.size();
    Array<Real> particle_ux(np, 0.0);
    Array<Real> particle_uy(np, 0.0);
    Array<Real> particle_uz(np, 0.0);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Adding laser particles\n";
    }
    AddNParticles(np, particle_x.data(), particle_y.data(), particle_z.data(),
		  particle_ux.data(), particle_uy.data(), particle_uz.data(),
		  1, particle_w.data(), 1);
}

void
LaserParticleContainer::Evolve (int lev,
				const MultiFab&, const MultiFab&, const MultiFab&,
				const MultiFab&, const MultiFab&, const MultiFab&,
				MultiFab& jx, MultiFab& jy, MultiFab& jz, Real t, Real dt)
{
    BL_PROFILE("Laser::Evolve()");
    BL_PROFILE_VAR_NS("Laser::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::LaserParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::LaserCurrentDepo", blp_pxr_cd);

    const Geometry& gm  = GDB().Geom(lev);
    const BoxArray& ba  = jx.boxArray();

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif

#if (BL_SPACEDIM == 3)
    long ngx_j  = jx.nGrow();
    long ngy_j  = ngx_j;
    long ngz_j  = ngx_j;
#elif (BL_SPACEDIM == 2)
    long ngx_j  = jx.nGrow();;
    long ngy_j  = 0;
    long ngz_j  = ngx_j;
#endif

    BL_ASSERT(OnSameGrids(lev,jx));

    {
	Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, giv, plane_Xp, plane_Yp;

	PartIterInfo info {lev, do_tiling, tile_size};
	for (PartIter pti(*this, info); pti.isValid(); ++pti)
	{
	    const int  gid = pti.index();
	    const Box& tbx = pti.tilebox();
	    const Box& vbx = pti.validbox();
	    const long np  = pti.numParticles();

	    // Data on the grid
	    FArrayBox& jxfab = jx[gid];
	    FArrayBox& jyfab = jy[gid];
	    FArrayBox& jzfab = jz[gid];

	    xp.resize(np);
	    yp.resize(np);
	    zp.resize(np);
	    wp.resize(np);
	    uxp.resize(np);
	    uyp.resize(np);
	    uzp.resize(np);
	    giv.resize(np);

	    //
	    // copy data from particle container to temp arrays
	    //
	    BL_PROFILE_VAR_START(blp_copy);
	    pti.foreach([&](int i, ParticleType& p) {
#if (BL_SPACEDIM == 3)
		    xp[i] = p.m_pos[0];
		    yp[i] = p.m_pos[1];
		    zp[i] = p.m_pos[2];
#elif (BL_SPACEDIM == 2)
		    xp[i] = p.m_pos[0];
		    yp[i] = std::numeric_limits<Real>::quiet_NaN();
		    zp[i] = p.m_pos[1];
#endif
		    wp[i]  = p.m_data[PIdx::w];
		    uxp[i] = p.m_data[PIdx::ux];
		    uyp[i] = p.m_data[PIdx::uy];
		    uzp[i] = p.m_data[PIdx::uz];
		});
	    BL_PROFILE_VAR_STOP(blp_copy);


	    const Box& box = BoxLib::enclosedCells(ba[gid]);
	    BL_ASSERT(box == vbx);
#if (BL_SPACEDIM == 3)
	    long nx = box.length(0);
	    long ny = box.length(1);
	    long nz = box.length(2);
#elif (BL_SPACEDIM == 2)
	    long nx = box.length(0);
	    long ny = 0;
	    long nz = box.length(1);
#endif
	    RealBox grid_box = RealBox( box, gm.CellSize(), gm.ProbLo() );
#if (BL_SPACEDIM == 3)
	    const Real* xyzmin = grid_box.lo();
#elif (BL_SPACEDIM == 2)
	    Real xyzmin[3] = { grid_box.lo(0), std::numeric_limits<Real>::quiet_NaN(), grid_box.lo(1) };
#endif

	    //
	    // Particle Push
	    //
      // Find the coordinates of the particles in the emission plane
      pti.foreach([&](int i, ParticleType& p) {
#if (BL_SPACEDIM == 3)
          plane_Xp[i] = u_X[0]*(xp[i] - position[0])
                      + u_X[1]*(yp[i] - position[1])
                      + u_X[2]*(zp[i] - position[2]);
          plane_Yp[i] = u_Y[0]*(xp[i] - position[0])
                      + u_Y[1]*(yp[i] - position[1])
                      + u_Y[2]*(zp[i] - position[2]);
#elif (BL_SPACEDIM == 2)
          plane_Xp[i] = u_X[0]*(xp[i] - position[0])
                      + u_X[2]*(zp[i] - position[2]);
          plane_Yp[i] = 0;
#endif
		  });
      // Calculate the laser amplitude to be emitted,
      // at the position of the emission plane
//      if (profile == laser_t::Gaussian) {
//        warpx_gaussian_pulse( plane_Xp, plane_Yp, t,
//                              profile_waist, profile_duration );
//      }
      // Calculate the corresponding momentum for the particles
//      pti.foreach([&](int i, ParticleType& p) {
//          v_over_c =
//          gamma = 1./( 1 - v_over_c**2 )
//
//		  });


	    BL_PROFILE_VAR_START(blp_pxr_pp);
	    warpx_laser_pusher(&np, xp.data(), yp.data(), zp.data(),
			       uxp.data(), uyp.data(), uzp.data(), giv.data(),
			       &this->charge, &dt,
			       &pusher_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_pp);

	    //
	    // Current Deposition
	    // xxxxx this part needs to be thread safe if we have OpenMP over tiles
	    //
	    long lvect = 8;
	    BL_PROFILE_VAR_START(blp_pxr_cd);
	    warpx_current_deposition(jxfab.dataPtr(), jyfab.dataPtr(), jzfab.dataPtr(),
				     &np, xp.data(), yp.data(), zp.data(),
				     uxp.data(), uyp.data(), uzp.data(),
				     giv.data(), wp.data(), &this->charge,
				     &xyzmin[0], &xyzmin[1], &xyzmin[2],
				     &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ngx_j, &ngy_j, &ngz_j,
                                     &WarpX::nox,&WarpX::noy,&WarpX::noz,
				     &lvect,&WarpX::current_deposition_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_cd);

	    //
	    // copy particle data back
	    //
	    BL_PROFILE_VAR_START(blp_copy);
	    pti.foreach([&](int i, ParticleType& p) {
                    BL_ASSERT(p.m_id > 0);
#if (BL_SPACEDIM == 3)
		    p.m_pos[0] = xp[i];
		    p.m_pos[1] = yp[i];
		    p.m_pos[2] = zp[i];
#elif (BL_SPACEDIM == 2)
		    p.m_pos[0] = xp[i];
		    p.m_pos[1] = zp[i];
#endif
		    p.m_data[PIdx::ux] = uxp[i];
		    p.m_data[PIdx::uy] = uyp[i];
		    p.m_data[PIdx::uz] = uzp[i];
                });
            BL_PROFILE_VAR_STOP(blp_copy);
	}
    }
}

Real
LaserParticleContainer::ComputeWeight (Real Sx, Real Sy) const
{
    constexpr Real eps = 0.1;
    constexpr Real fac = 1.0/(2.0*3.1415926535897932*PhysConst::mu0*PhysConst::c*PhysConst::c*eps);
    return fac * wavelength * Sx * Sy / std::min(Sx,Sy) * e_max;
}
