#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>

extern "C" 
{
	// Charge deposition
	void warpx_charge_deposition(amrex::Real* rho, 
                        const long* np,	const amrex::Real* xp, const amrex::Real* yp, const amrex::Real* zp,	const amrex::Real* w,
			const amrex::Real* q, const amrex::Real* xmin, const amrex::Real* ymin, const amrex::Real* zmin,
			const amrex::Real* dx, const amrex::Real* dy, const amrex::Real* dz,
			const long* nx, const long* ny, const long* nz,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nox, const long* noy,const long* noz,
			const long* lvect, const long* charge_depo_algo);

	// Current deposition

	void warpx_current_deposition(amrex::Real* jx, amrex::Real* jy, amrex::Real* jz, const long* np,
			const amrex::Real* xp, const amrex::Real* yp, const amrex::Real* zp,
			const amrex::Real* uxp, const amrex::Real* uyp,const amrex::Real* uzp,
			const amrex::Real* gip, const amrex::Real* w, const amrex::Real* q,
			const amrex::Real* xmin, const amrex::Real* ymin, const amrex::Real* zmin,
			const amrex::Real* dt, 
			const amrex::Real* dx, const amrex::Real* dy, const amrex::Real* dz,
			const long* nx, const long* ny, const long* nz,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nox, const long* noy,const long* noz,
			const long* lvect, const long* current_depo_algo);

	// Field gathering

	void warpx_geteb_energy_conserving(const long* np, 
			const amrex::Real* xp, const amrex::Real* yp, const amrex::Real* zp,
			const amrex::Real* exp, const amrex::Real* eyp,const amrex::Real* ezp,
			const amrex::Real* bxp, const amrex::Real* byp,const amrex::Real* bzp,
			const amrex::Real* xmin, const amrex::Real* ymin, const amrex::Real* zmin,
			const amrex::Real* dx, const amrex::Real* dy, const amrex::Real* dz,
			const long* nx, const long* ny, const long* nz,
			const long* nxguard, const long* nyguard, const long* nzguard,
			const long* nox, const long* noy, const long* noz,
			const amrex::Real* exg, const amrex::Real* eyg, const amrex::Real* ezg,
			const amrex::Real* bxg, const amrex::Real* byg, const amrex::Real* bzg,
			const int* ll4symtry, const int* l_lower_order_in_v,
			const long* lvect,  
			const long* field_gathe_algo);

	// Particle pusher

	void warpx_particle_pusher(const long* np,
			const amrex::Real* xp, const amrex::Real* yp, const amrex::Real* zp,
			const amrex::Real* uxp, const amrex::Real* uyp, const amrex::Real* uzp,
			const amrex::Real* gaminv,
			const amrex::Real* exp, const amrex::Real* eyp,const amrex::Real* ezp,
			const amrex::Real* bxp, const amrex::Real* byp,const amrex::Real* bzp,
			const amrex::Real* charge, const amrex::Real* mass, amrex::Real* dt,
			const long* particle_pusher_algo);

	// Maxwell solver

	void warpx_push_evec( amrex::Real* ex, amrex::Real* ey, amrex::Real* ez,
			amrex::Real* bx, amrex::Real* by, amrex::Real* bz,
			const amrex::Real* jx, const amrex::Real* jy, const amrex::Real* jz,
			const amrex::Real* mudt,
			const amrex::Real* dtsdx, const amrex::Real* dtsdy, const amrex::Real* dtsdz, 
			const long* nxlocal, const long* nylocal, const long* nzlocal,
			const long* nxorder, const long* nyorder, const long* nzorder,
			const long* nxguard, const long* nyguard, const long* nzguard, 
			const long* nxs, const long* nys, const long* nzs,
			const int* l_nodal );

	void warpx_push_bvec( amrex::Real* ex, amrex::Real* ey, amrex::Real* ez,
			amrex::Real* bx, amrex::Real* by, amrex::Real* bz,
			const amrex::Real* dtsdx, const amrex::Real* dtsdy, const amrex::Real* dtsdz, 
			const long* nxlocal, const long* nylocal, const long* nzlocal,
			const long* nxorder, const long* nyorder, const long* nzorder,
			const long* nxguard, const long* nyguard, const long* nzguard, 
			const long* nxs, const long* nys, const long* nzs,
			const int* l_nodal );

#ifdef USE_OPENBC_POISSON
    void warpx_openbc_potential (amrex::Real* rho, amrex::Real* phi, const amrex::Real* dx);
    void warpx_openbc_decompose (const int*, const int*, int*, int*);
#endif

    void warpx_compute_E (const int* lo, const int* hi,
			  const BL_FORT_FAB_ARG_3D(phi),
			  BL_FORT_FAB_ARG_3D(Ex),
			  BL_FORT_FAB_ARG_3D(Ey),
			  BL_FORT_FAB_ARG_3D(Ez),
			  const amrex::Real* dx);
			  			  
}
