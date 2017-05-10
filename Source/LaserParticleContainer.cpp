
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <ParticleContainer.H>

using namespace amrex;

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

	// Parse the type of laser profile and set the corresponding flag `profile`
	std::string laser_type_s;
	pp.get("profile", laser_type_s);
	std::transform(laser_type_s.begin(), laser_type_s.end(), laser_type_s.begin(), ::tolower);
	if (laser_type_s == "gaussian") {
	    profile = laser_t::Gaussian;
	} else {
	    amrex::Abort("Unknown laser type");
	}

	// Parse the properties of the antenna
	pp.getarr("position", position);
	pp.getarr("direction", nvec);
	pp.getarr("polarization", p_X);
	pp.query("pusher_algo", pusher_algo);
	pp.get("e_max", e_max);
	pp.get("wavelength", wavelength);

	if ( profile == laser_t::Gaussian ) {
	    // Parse the properties of the Gaussian profile
	   pp.get("profile_waist", profile_waist);
	   pp.get("profile_duration", profile_duration);
	   pp.get("profile_t_peak", profile_t_peak);
	   pp.get("profile_focal_distance", profile_focal_distance);
	}

	// Plane normal
	Real s = 1.0/std::sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
	nvec = { nvec[0]*s, nvec[1]*s, nvec[2]*s };

	// The first polarization vector
	s = 1.0/std::sqrt(p_X[0]*p_X[0] + p_X[1]*p_X[1] + p_X[2]*p_X[2]);
	p_X = { p_X[0]*s, p_X[1]*s, p_X[2]*s };

	Real dp = std::inner_product(nvec.begin(), nvec.end(), p_X.begin(), 0.0);
	if (std::abs(dp) > 1.e-14) {
	    amrex::Abort("Laser plane vector is not perpendicular to the main polarization vector");
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
LaserParticleContainer::InitData ()
{
    for (int lev = 0; lev < maxLevel(); ++lev) {
        InitData(lev);
    }
}

void
LaserParticleContainer::InitData (int lev)
{
    const Geometry& geom = Geom(lev);
    const RealBox& prob_domain = geom.ProbDomain();

    // spacing of laser particles in the laser plane.
    // has to be done after geometry is set up.
    Real S_X, S_Y;
    ComputeSpacing(lev, S_X, S_Y);
    ComputeWeightMobility(S_X, S_Y);

    auto Transform = [&](int i, int j) -> Array<Real>
    {
#if (BL_SPACEDIM == 3)
	return { position[0] + (S_X*(i+0.5))*u_X[0] + (S_Y*(j+0.5))*u_Y[0],
		 position[1] + (S_X*(i+0.5))*u_X[1] + (S_Y*(j+0.5))*u_Y[1],
		 position[2] + (S_X*(i+0.5))*u_X[2] + (S_Y*(j+0.5))*u_Y[2] };
#else
	return { position[0] + (S_X*(i+0.5))*u_X[0],
		 0.0,
		 position[2] + (S_X*(i+0.5))*u_X[2] };
#endif
    };

    // Given the "lab" frame coordinates, return the real coordinates in the laser plane coordinates
    auto InverseTransform = [&](const Array<Real>& pos) -> Array<Real>
    {
#if (BL_SPACEDIM == 3)
	return {u_X[0]*(pos[0]-position[0])+u_X[1]*(pos[1]-position[1])+u_X[2]*(pos[2]-position[2]),
		u_Y[0]*(pos[0]-position[0])+u_Y[1]*(pos[1]-position[1])+u_Y[2]*(pos[2]-position[2])};
#else
	return {u_X[0]*(pos[0]-position[0])+u_X[2]*(pos[2]-position[2]), 0.0};
#endif
    };

    Array<int> plane_lo(2, std::numeric_limits<int>::max());
    Array<int> plane_hi(2, std::numeric_limits<int>::min());
    {
	auto compute_min_max = [&](Real x, Real y, Real z)
	{
	    const Array<Real>& pos_plane = InverseTransform({x, y, z});
	    int i = pos_plane[0]/S_X;
	    int j = pos_plane[1]/S_Y;
	    plane_lo[0] = std::min(plane_lo[0], i);
	    plane_lo[1] = std::min(plane_lo[1], j);
	    plane_hi[0] = std::max(plane_hi[0], i);
	    plane_hi[1] = std::max(plane_hi[1], j);
	};

	const Real* prob_lo = prob_domain.lo();
	const Real* prob_hi = prob_domain.hi();
#if (BL_SPACEDIM == 3)
	compute_min_max(prob_lo[0], prob_lo[1], prob_lo[2]);
	compute_min_max(prob_hi[0], prob_lo[1], prob_lo[2]);
	compute_min_max(prob_lo[0], prob_hi[1], prob_lo[2]);
	compute_min_max(prob_hi[0], prob_hi[1], prob_lo[2]);
	compute_min_max(prob_lo[0], prob_lo[1], prob_hi[2]);
	compute_min_max(prob_hi[0], prob_lo[1], prob_hi[2]);
	compute_min_max(prob_lo[0], prob_hi[1], prob_hi[2]);
	compute_min_max(prob_hi[0], prob_hi[1], prob_hi[2]);
#else
	compute_min_max(prob_lo[0], 0.0, prob_lo[1]);
	compute_min_max(prob_hi[0], 0.0, prob_lo[1]);
	compute_min_max(prob_lo[0], 0.0, prob_hi[1]);
	compute_min_max(prob_hi[0], 0.0, prob_hi[1]);
#endif
    }

    const int nprocs = ParallelDescriptor::NProcs();
    const int myproc = ParallelDescriptor::MyProc();

#if (BL_SPACEDIM == 3)
    const Box plane_box {IntVect(plane_lo[0],plane_lo[1],0),
                         IntVect(plane_hi[0],plane_hi[1],0)};
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
#else
    BoxArray plane_ba { Box {IntVect(plane_lo[0],0), IntVect(plane_hi[0],0)} };
#endif

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
#if (BL_SPACEDIM == 3)
		const Real* x = pos.data();
#else
		const Real x[2] = {pos[0], pos[2]};
#endif
		if (prob_domain.contains(x))
		{
		    for (int k = 0; k<2; ++k) {
			particle_x.push_back(pos[0]);
			particle_y.push_back(pos[1]);
			particle_z.push_back(pos[2]);
		    }
		    particle_w.push_back( weight);
		    particle_w.push_back(-weight);
		}
	    }
	}
    }
    const int np = particle_z.size();
    Array<Real> particle_ux(np, 0.0);
    Array<Real> particle_uy(np, 0.0);
    Array<Real> particle_uz(np, 0.0);

    if (Verbose()) amrex::Print() << "Adding laser particles\n";
    AddNParticles(lev,
                  np, particle_x.data(), particle_y.data(), particle_z.data(),
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

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    // WarpX assumes the same number of guard cells for Jx, Jy, Jz
    long ngJ  = jx.nGrow();

    BL_ASSERT(OnSameGrids(lev,jx));

    {
	Array<Real> xp, yp, zp, giv, plane_Xp, plane_Yp, amplitude_E;

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
	    const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto&  wp = attribs[PIdx::w ];
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];

	    const long np  = pti.numParticles();

	    // Data on the grid
	    FArrayBox& jxfab = jx[pti];
	    FArrayBox& jyfab = jy[pti];
	    FArrayBox& jzfab = jz[pti];

	            giv.resize(np);
               plane_Xp.resize(np);
               plane_Yp.resize(np);
            amplitude_E.resize(np);

	    //
	    // copy data from particle container to temp arrays
	    //
	    BL_PROFILE_VAR_START(blp_copy);
            pti.GetPosition(xp, yp, zp);
	    BL_PROFILE_VAR_STOP(blp_copy);

	    for (int i = 0; i < np; ++i)
            {
                // Find the coordinates of the particles in the emission plane
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
            }

            const std::array<Real,3>& xyzmin = WarpX::LowerCorner(box, lev);

	    //
	    // Particle Push
	    //
	    BL_PROFILE_VAR_START(blp_pxr_pp);
	    // Calculate the laser amplitude to be emitted,
	    // at the position of the emission plane
	    if (profile == laser_t::Gaussian) {
		warpx_gaussian_laser( &np, plane_Xp.data(), plane_Yp.data(),
				      &t, &wavelength, &e_max, &profile_waist, &profile_duration,
				      &profile_t_peak, &profile_focal_distance, amplitude_E.data() );
	    }

	    // Calculate the corresponding momentum and position for the particles
            for (int i = 0; i < np; ++i)
            {
                // Calculate the velocity according to the amplitude of E
                Real sign_charge = std::copysign( 1.0, wp[i] );
                Real v_over_c = sign_charge * mobility * amplitude_E[i];
                BL_ASSERT( v_over_c < 1 );
                giv[i] = std::sqrt( 1 - v_over_c * v_over_c );
                Real gamma = 1./giv[i];
                // The velocity is along the laser polarization p_X
                Real vx = PhysConst::c * v_over_c * p_X[0];
                Real vy = PhysConst::c * v_over_c * p_X[1];
                Real vz = PhysConst::c * v_over_c * p_X[2];
                // Get the corresponding momenta
                uxp[i] = gamma * vx;
                uyp[i] = gamma * vy;
                uzp[i] = gamma * vz;
                // Push the the particle positions
                xp[i] += vx * dt;
#if (BL_SPACEDIM == 3)
                yp[i] += vy * dt;
#endif
                zp[i] += vz * dt;
            }

	    BL_PROFILE_VAR_STOP(blp_pxr_pp);

	    //
	    // Current Deposition
	    // xxxxx this part needs to be thread safe if we have OpenMP over tiles
	    //
	    long lvect = 8;
	    BL_PROFILE_VAR_START(blp_pxr_cd);
	    warpx_current_deposition(
        jxfab.dataPtr(), &ngJ, jxfab.length(),
        jyfab.dataPtr(), &ngJ, jyfab.length(),
        jzfab.dataPtr(), &ngJ, jzfab.length(),
	    &np, xp.data(), yp.data(), zp.data(),
	    uxp.data(), uyp.data(), uzp.data(),
	    giv.data(), wp.data(), &this->charge,
	    &xyzmin[0], &xyzmin[1], &xyzmin[2],
	    &dt, &dx[0], &dx[1], &dx[2],
	    &WarpX::nox,&WarpX::noy,&WarpX::noz,
	    &lvect,&WarpX::current_deposition_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_cd);

	    //
	    // copy particle data back
	    //
	    BL_PROFILE_VAR_START(blp_copy);
            pti.SetPosition(xp, yp, zp);
            BL_PROFILE_VAR_STOP(blp_copy);
	}
    }

}

void
LaserParticleContainer::PostRestart ()
{
    Real Sx, Sy;
    const int lev = 0; // xxxxx
    ComputeSpacing(lev, Sx, Sy);
    ComputeWeightMobility(Sx, Sy);
}

void
LaserParticleContainer::ComputeSpacing (int lev, Real& Sx, Real& Sy) const
{
    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    const Real eps = dx[0]*1.e-50;
#if (BL_SPACEDIM == 3)
    Sx = std::min(std::min(dx[0]/(std::abs(u_X[0])+eps),
			   dx[1]/(std::abs(u_X[1])+eps)),
		           dx[2]/(std::abs(u_X[2])+eps));
    Sy = std::min(std::min(dx[0]/(std::abs(u_Y[0])+eps),
			   dx[1]/(std::abs(u_Y[1])+eps)),
		           dx[2]/(std::abs(u_Y[2])+eps));
#else
    Sx = std::min(dx[0]/(std::abs(u_X[0])+eps),
		  dx[2]/(std::abs(u_X[2])+eps));
    Sy = 1.0;
#endif
}

void
LaserParticleContainer::ComputeWeightMobility (Real Sx, Real Sy)
{
    constexpr Real eps = 0.1;
    constexpr Real fac = 1.0/(2.0*3.1415926535897932*PhysConst::mu0*PhysConst::c*PhysConst::c*eps);
    weight = fac * wavelength * Sx * Sy / std::min(Sx,Sy) * e_max;

    // The mobility is the constant of proportionality between the field to
    // be emitted, and the corresponding velocity that the particles need to have.
    mobility = (Sx * Sy)/(weight * PhysConst::mu0 * PhysConst::c * PhysConst::c);
}
