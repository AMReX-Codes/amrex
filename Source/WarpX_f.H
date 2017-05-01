
#include <AMReX_BLFort.H>

#if (BL_SPACEDIM == 3)

#define WRPX_PXR_PUSH_BVEC               warpx_pxr_push_em3d_bvec
#define WRPX_PXR_PUSH_EVEC               warpx_pxr_push_em3d_evec
#define WRPX_PUSH_PML_BVEC               warpx_push_pml_bvec_3d
#define WRPX_PUSH_PML_EVEC               warpx_push_pml_evec_3d

#elif (BL_SPACEDIM == 2)

#define WRPX_PXR_PUSH_BVEC               warpx_pxr_push_em2d_bvec
#define WRPX_PXR_PUSH_EVEC               warpx_pxr_push_em2d_evec
#define WRPX_PUSH_PML_BVEC               warpx_push_pml_bvec_2d
#define WRPX_PUSH_PML_EVEC               warpx_push_pml_evec_2d

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
	void warpx_current_deposition(
			amrex::Real* jx, const long* jx_ng, const int* jx_ntot,
			amrex::Real* jy, const long* jy_ng, const int* jy_ntot,
			amrex::Real* jz, const long* jz_ng, const int* jz_ntot,
			const long* np,
			const amrex::Real* xp, const amrex::Real* yp, const amrex::Real* zp,
			const amrex::Real* uxp, const amrex::Real* uyp,const amrex::Real* uzp,
			const amrex::Real* gip, const amrex::Real* w, const amrex::Real* q,
			const amrex::Real* xmin, const amrex::Real* ymin, const amrex::Real* zmin,
			const amrex::Real* dt,
			const amrex::Real* dx, const amrex::Real* dy, const amrex::Real* dz,
			const long* nox, const long* noy,const long* noz,
			const long* lvect, const long* current_depo_algo);

	// Field gathering

	void warpx_geteb_energy_conserving(const long* np,
			const amrex::Real* xp, const amrex::Real* yp, const amrex::Real* zp,
			amrex::Real* exp, amrex::Real* eyp, amrex::Real* ezp,
		    amrex::Real* bxp, amrex::Real* byp, amrex::Real* bzp,
			const amrex::Real* xmin, const amrex::Real* ymin, const amrex::Real* zmin,
			const amrex::Real* dx, const amrex::Real* dy, const amrex::Real* dz,
			const long* nox, const long* noy, const long* noz,
			const amrex::Real* exg, const long* exg_ng, const int* exg_ntot,
			const amrex::Real* eyg, const long* eyg_ng, const int* eyg_ntot,
			const amrex::Real* ezg, const long* ezg_ng, const int* ezg_ntot,
			const amrex::Real* bxg, const long* bxg_ng, const int* bxg_ntot,
			const amrex::Real* byg, const long* byg_ng, const int* byg_ntot,
			const amrex::Real* bzg, const long* bzg_ng, const int* bzg_ntot,
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

        void WRPX_PXR_PUSH_EVEC(
				const int* xlo, const int* xhi,
				const int* ylo, const int* yhi,
				const int* zlo, const int* zhi,
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

        void WRPX_PXR_PUSH_BVEC(const int* xlo, const int* xhi,
				const int* ylo, const int* yhi,
				const int* zlo, const int* zhi,
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

//  These functions are used in electrostatic mode

    void warpx_compute_E_nodal (const int* lo, const int* hi,
                                const amrex::Real* phi,
                                amrex::Real* Ex, amrex::Real* Ey, 
                                amrex::Real* Ez, amrex::Real* dx);
    
    void warpx_deposit_cic(const amrex::Real* particles, int ns, int np,
                           const amrex::Real* weights, 
                           const amrex::Real* charge,
                           amrex::Real* rho, const int* lo, const int* hi,
                           const amrex::Real* plo, const amrex::Real* dx);
    
    void warpx_interpolate_cic(const amrex::Real* particles, int ns, int np,
                               amrex::Real* Ex_p, amrex::Real* Ey_p, amrex::Real* Ez_p,
                               const amrex::Real* Ex, const amrex::Real* Ey, const amrex::Real* Ez,
                               const int* lo, const int* hi,
                               const amrex::Real* plo, const amrex::Real* dx);

    void warpx_push_leapfrog(amrex::Real* particles, int ns, int np, 
                             amrex::Real* vx_p, amrex::Real* vy_p, amrex::Real* vz_p,
                             const amrex::Real* Ex_p, const amrex::Real* Ey_p, const amrex::Real* Ez_p,
                             const amrex::Real* charge, const amrex::Real* mass, const amrex::Real* dt);

    void warpx_push_leapfrog_positions(amrex::Real* particles, int ns, int np, 
                                       const amrex::Real* vx_p, const amrex::Real* vy_p, 
                                       const amrex::Real* vz_p, const amrex::Real* dt);

//  These functions are used to evolve E and B in the PML

    void WRPX_PUSH_PML_BVEC(const int* xlo, const int* xhi,
                            const int* ylo, const int* yhi,
                            const int* zlo, const int* zhi,
                            const BL_FORT_FAB_ARG_3D(ex),
                            const BL_FORT_FAB_ARG_3D(ey),
                            const BL_FORT_FAB_ARG_3D(ez),
                            BL_FORT_FAB_ARG_3D(bx),
                            BL_FORT_FAB_ARG_3D(by),
                            BL_FORT_FAB_ARG_3D(bz),
                            const amrex::Real* sigx1, int sigx1_lo, int sigx1_hi, 
                            const amrex::Real* sigx2, int sigx2_lo, int sigx2_hi,
#if (BL_SPACEDIM == 3)
                            const amrex::Real* sigy1, int sigy1_lo, int sigy1_hi, 
                            const amrex::Real* sigy2, int sigy2_lo, int sigy2_hi, 
#endif
                            const amrex::Real* sigz1, int sigz1_lo, int sigz1_hi, 
                            const amrex::Real* sigz2, int sigz2_lo, int sigz2_hi);

    void WRPX_PUSH_PML_EVEC(const int* xlo, const int* xhi,
                            const int* ylo, const int* yhi,
                            const int* zlo, const int* zhi,
                            BL_FORT_FAB_ARG_3D(ex),
                            BL_FORT_FAB_ARG_3D(ey),
                            BL_FORT_FAB_ARG_3D(ez),
                            const BL_FORT_FAB_ARG_3D(bx),
                            const BL_FORT_FAB_ARG_3D(by),
                            const BL_FORT_FAB_ARG_3D(bz),
                            const amrex::Real* sigx1, int sigx1_lo, int sigx1_hi, 
                            const amrex::Real* sigx2, int sigx2_lo, int sigx2_hi, 
#if (BL_SPACEDIM == 3)
                            const amrex::Real* sigy1, int sigy1_lo, int sigy1_hi, 
                            const amrex::Real* sigy2, int sigy2_lo, int sigy2_hi, 
#endif
                            const amrex::Real* sigz1, int sigz1_lo, int sigz1_hi, 
                            const amrex::Real* sigz2, int sigz2_lo, int sigz2_hi);
}
