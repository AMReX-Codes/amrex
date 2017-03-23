
#include <AMReX_BLFort.H>

#if (BL_SPACEDIM == 3)

#define WRPX_PXR_PUSH_BVEC               warpx_pxr_push_em3d_bvec
#define WRPX_PXR_PUSH_EVEC               warpx_pxr_push_em3d_evec

#elif (BL_SPACEDIM == 2)

#define WRPX_PXR_PUSH_BVEC               warpx_pxr_push_em2d_bvec
#define WRPX_PXR_PUSH_EVEC               warpx_pxr_push_em2d_evec

#endif

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
			amrex::Real* exp, amrex::Real* eyp, amrex::Real* ezp,
		        amrex::Real* bxp, amrex::Real* byp, amrex::Real* bzp,
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

	// Particle pusher (velocity and position)

	void warpx_particle_pusher(const long* np,
			amrex::Real* xp, amrex::Real* yp, amrex::Real* zp,
			amrex::Real* uxp, amrex::Real* uyp, amrex::Real* uzp, amrex::Real* gaminv,
			const amrex::Real* exp, const amrex::Real* eyp,const amrex::Real* ezp,
			const amrex::Real* bxp, const amrex::Real* byp,const amrex::Real* bzp,
			const amrex::Real* charge, const amrex::Real* mass, const amrex::Real* dt,
			const long* particle_pusher_algo);

	// Particle pusher (position)

	void warpx_particle_pusher_positions(const long* np,
                                             amrex::Real* xp, amrex::Real* yp, amrex::Real* zp,
                                             amrex::Real* uxp, amrex::Real* uyp, amrex::Real* uzp, amrex::Real* gaminv,
                                             const amrex::Real* dt);

    // Laser pusher

	void warpx_gaussian_laser( const long* np,
				amrex::Real* Xp, amrex::Real* Yp, amrex::Real* t, amrex::Real* wavelength, amrex::Real* e_max, amrex::Real* waist,
				amrex::Real* duration, amrex::Real* t_peak, amrex::Real* f, amrex::Real* amplitude );

	// Maxwell solver

        void WRPX_PXR_PUSH_EVEC(const int* lo, const int* hi,
				BL_FORT_FAB_ARG_3D(ex),
				BL_FORT_FAB_ARG_3D(ey),
				BL_FORT_FAB_ARG_3D(ez),
				const BL_FORT_FAB_ARG_3D(bx),
				const BL_FORT_FAB_ARG_3D(by),
				const BL_FORT_FAB_ARG_3D(bz),
				const BL_FORT_FAB_ARG_3D(jx),
				const BL_FORT_FAB_ARG_3D(jy),
				const BL_FORT_FAB_ARG_3D(jz),
				const amrex::Real* mudt,
				const amrex::Real* dtsdx, const amrex::Real* dtsdy, const amrex::Real* dtsdz,
				const int* norder);

        void WRPX_PXR_PUSH_BVEC(const int* lo, const int* hi,
				const BL_FORT_FAB_ARG_3D(ex),
				const BL_FORT_FAB_ARG_3D(ey),
				const BL_FORT_FAB_ARG_3D(ez),
				BL_FORT_FAB_ARG_3D(bx),
				BL_FORT_FAB_ARG_3D(by),
				BL_FORT_FAB_ARG_3D(bz),
				const amrex::Real* dtsdx, const amrex::Real* dtsdy, const amrex::Real* dtsdz,
				const int* norder);

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
