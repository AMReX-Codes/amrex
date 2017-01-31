#include <REAL.H>
#include <BLFort.H>

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
				Real* Xp, Real* Yp, Real* t, Real* wavelength, Real* e_max, Real* waist,
				Real* duration, Real* t_peak, Real* f, Real* amplitude );

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
				const Real* mudt,
				const Real* dtsdx, const Real* dtsdy, const Real* dtsdz,
				const int* norder);

        void WRPX_PXR_PUSH_BVEC(const int* lo, const int* hi,
				const BL_FORT_FAB_ARG_3D(ex),
				const BL_FORT_FAB_ARG_3D(ey),
				const BL_FORT_FAB_ARG_3D(ez),
				BL_FORT_FAB_ARG_3D(bx),
				BL_FORT_FAB_ARG_3D(by),
				BL_FORT_FAB_ARG_3D(bz),
				const Real* dtsdx, const Real* dtsdy, const Real* dtsdz,
				const int* norder);

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
