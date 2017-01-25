#include <REAL.H>
#include <BLFort.H>

extern "C"
{
	// Charge deposition
	void warpx_charge_deposition(Real* rho,
                        const long* np,	const Real* xp, const Real* yp, const Real* zp,	const Real* w,
			const Real* q, const Real* xmin, const Real* ymin, const Real* zmin,
			const Real* dx, const Real* dy, const Real* dz,
			const long* nx, const long* ny, const long* nz,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nox, const long* noy,const long* noz,
			const long* lvect, const long* charge_depo_algo);

	// Current deposition

	void warpx_current_deposition(Real* jx, Real* jy, Real* jz, const long* np,
			const Real* xp, const Real* yp, const Real* zp,
			const Real* uxp, const Real* uyp,const Real* uzp,
			const Real* gip, const Real* w, const Real* q,
			const Real* xmin, const Real* ymin, const Real* zmin,
			const Real* dt,
			const Real* dx, const Real* dy, const Real* dz,
			const long* nx, const long* ny, const long* nz,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nox, const long* noy,const long* noz,
			const long* lvect, const long* current_depo_algo);

	// Field gathering

	void warpx_geteb_energy_conserving(const long* np,
			const Real* xp, const Real* yp, const Real* zp,
			Real* exp, Real* eyp, Real* ezp,
		        Real* bxp, Real* byp, Real* bzp,
			const Real* xmin, const Real* ymin, const Real* zmin,
			const Real* dx, const Real* dy, const Real* dz,
			const long* nx, const long* ny, const long* nz,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nox, const long* noy, const long* noz,
			const Real* exg, const Real* eyg, const Real* ezg,
			const Real* bxg, const Real* byg, const Real* bzg,
			const int* ll4symtry, const int* l_lower_order_in_v,
			const long* lvect,
			const long* field_gathe_algo);

	// Particle pusher

	void warpx_particle_pusher(const long* np,
			Real* xp, Real* yp, Real* zp,
			Real* uxp, Real* uyp, Real* uzp, Real* gaminv,
			const Real* exp, const Real* eyp,const Real* ezp,
			const Real* bxp, const Real* byp,const Real* bzp,
			const Real* charge, const Real* mass, const Real* dt,
			const long* particle_pusher_algo);

  // Laser pusher

	void warpx_gaussian_laser( const long* np,
				Real* Xp, Real* Yp, Real* t, Real* wavelength, Real* waist,
				Real* duration, Real* t_peak, Real* f, Real* amplitude );

	// Maxwell solver

	void warpx_push_evec( Real* ex, Real* ey, Real* ez,
			const Real* bx, const Real* by, const Real* bz,
			const Real* jx, const Real* jy, const Real* jz,
			const Real* mudt,
			const Real* dtsdx, const Real* dtsdy, const Real* dtsdz,
			const long* nxlocal, const long* nylocal, const long* nzlocal,
			const long* nxorder, const long* nyorder, const long* nzorder,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nxs, const long* nys, const long* nzs,
			const int* l_nodal );

	void warpx_push_bvec( const Real* ex, const Real* ey, const Real* ez,
			Real* bx, Real* by, Real* bz,
			const Real* dtsdx, const Real* dtsdy, const Real* dtsdz,
			const long* nxlocal, const long* nylocal, const long* nzlocal,
			const long* nxorder, const long* nyorder, const long* nzorder,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nxs, const long* nys, const long* nzs,
			const int* l_nodal );

#ifdef USE_OPENBC_POISSON
    void warpx_openbc_potential (Real* rho, Real* phi, const Real* dx);
    void warpx_openbc_decompose (const int*, const int*, int*, int*);
#endif

    void warpx_compute_E (const int* lo, const int* hi,
			  const BL_FORT_FAB_ARG_3D(phi),
			  BL_FORT_FAB_ARG_3D(Ex),
			  BL_FORT_FAB_ARG_3D(Ey),
			  BL_FORT_FAB_ARG_3D(Ez),
			  const Real* dx);

}
