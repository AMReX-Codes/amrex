#include <AMReX_BLFort.H>

#ifdef __cplusplus

#if AMREX_SPACEDIM==2
#define WARPX_ARLIM_ANYD(x) std::array<int,3>{(x)[0], 0, (x)[1]}.data()
#else
#define WARPX_ARLIM_ANYD(x) x
#endif

#define WARPX_TO_FORTRAN_BOX(x) WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())
#define WARPX_TO_FORTRAN_ANYD(x) (x).dataPtr(), WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())
#define WARPX_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())

#endif

#if (AMREX_SPACEDIM == 3)

#define WRPX_SYNC_CURRENT                warpx_sync_current_3d
#define WRPX_SYNC_RHO                    warpx_sync_rho_3d

#define WRPX_PUSH_PML_BVEC               warpx_push_pml_bvec_3d
#define WRPX_PUSH_PML_EVEC               warpx_push_pml_evec_3d
#define WRPX_PUSH_PML_EVEC_F             warpx_push_pml_evec_f_3d
#define WRPX_PUSH_PML_F                  warpx_push_pml_f_3d
#define WRPX_DAMP_PML                    warpx_damp_pml_3d
#define WRPX_DAMP_PML_F                  warpx_damp_pml_f_3d
#define WRPX_SUM_FINE_TO_CRSE_NODAL      warpx_sum_fine_to_crse_nodal_3d
#define WRPX_ZERO_OUT_BNDRY              warpx_zero_out_bndry_3d
#define WRPX_BUILD_MASK                  warpx_build_mask_3d
#define WRPX_COMPUTE_E_NODAL             warpx_compute_E_nodal_3d
#define WRPX_DEPOSIT_CIC                 warpx_deposit_cic_3d
#define WRPX_INTERPOLATE_CIC             warpx_interpolate_cic_3d
#define WRPX_INTERPOLATE_CIC_TWO_LEVELS  warpx_interpolate_cic_two_levels_3d
#define WRPX_PUSH_LEAPFROG               warpx_push_leapfrog_3d
#define WRPX_PUSH_LEAPFROG_POSITIONS     warpx_push_leapfrog_positions_3d

#elif (AMREX_SPACEDIM == 2)

#define WRPX_SYNC_CURRENT                warpx_sync_current_2d
#define WRPX_SYNC_RHO                    warpx_sync_rho_2d

#define WRPX_PUSH_PML_BVEC               warpx_push_pml_bvec_2d
#define WRPX_PUSH_PML_EVEC               warpx_push_pml_evec_2d
#define WRPX_PUSH_PML_EVEC_F             warpx_push_pml_evec_f_2d
#define WRPX_PUSH_PML_F                  warpx_push_pml_f_2d
#define WRPX_DAMP_PML                    warpx_damp_pml_2d
#define WRPX_DAMP_PML_F                  warpx_damp_pml_f_2d

#define WRPX_SUM_FINE_TO_CRSE_NODAL      warpx_sum_fine_to_crse_nodal_2d
#define WRPX_ZERO_OUT_BNDRY              warpx_zero_out_bndry_2d
#define WRPX_BUILD_MASK                  warpx_build_mask_2d
#define WRPX_COMPUTE_E_NODAL             warpx_compute_E_nodal_2d
#define WRPX_DEPOSIT_CIC                 warpx_deposit_cic_2d
#define WRPX_INTERPOLATE_CIC             warpx_interpolate_cic_2d
#define WRPX_INTERPOLATE_CIC_TWO_LEVELS  warpx_interpolate_cic_two_levels_2d
#define WRPX_PUSH_LEAPFROG               warpx_push_leapfrog_2d
#define WRPX_PUSH_LEAPFROG_POSITIONS     warpx_push_leapfrog_positions_2d

#endif

#ifdef __cplusplus
extern "C"
{
#endif

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

///
/// These functions are used in electrostatic mode.
///

    void WRPX_SUM_FINE_TO_CRSE_NODAL(const int* lo, const int* hi, const int* lrat,
                                     amrex::Real* crse, const int* clo, const int* chi,
                                     const amrex::Real* fine, const int* flo, const int* fhi);

    void WRPX_ZERO_OUT_BNDRY(const int* lo, const int* hi,
                             amrex::Real* input_data, amrex::Real* bndry_data,
                             const int* mask);

    void WRPX_BUILD_MASK(const int* lo, const int* hi,
                         const int* tmp_mask, int* mask, const int* ncells);

    void WRPX_COMPUTE_E_NODAL(const int* lo, const int* hi,
                              const amrex::Real* phi,
                              amrex::Real* Ex, amrex::Real* Ey,
#if (AMREX_SPACEDIM == 3)
                          amrex::Real* Ez,
#endif
                          const amrex::Real* dx);

    void WRPX_DEPOSIT_CIC(const amrex::Real* particles, int ns, int np,
                          const amrex::Real* weights,
                          const amrex::Real* charge,
                          amrex::Real* rho, const int* lo, const int* hi,
                          const amrex::Real* plo, const amrex::Real* dx,
                          const int* ng);

    void WRPX_INTERPOLATE_CIC_TWO_LEVELS(const amrex::Real* particles, int ns, int np,
                                         amrex::Real* Ex_p, amrex::Real* Ey_p,
#if (AMREX_SPACEDIM == 3)
                                         amrex::Real* Ez_p,
#endif
                                         const amrex::Real* Ex, const amrex::Real* Ey,
#if (AMREX_SPACEDIM == 3)
                                         const amrex::Real* Ez,
#endif
                                         const int* lo, const int* hi, const amrex::Real* dx,
                                         const amrex::Real* cEx, const amrex::Real* cEy,
#if (AMREX_SPACEDIM == 3)
                                         const amrex::Real* cEz,
#endif
                                         const int* mask,
                                         const int* clo, const int* chi, const amrex::Real* cdx,
                                         const amrex::Real* plo, const int* ng, const int* lev);

    void WRPX_INTERPOLATE_CIC(const amrex::Real* particles, int ns, int np,
                              amrex::Real* Ex_p, amrex::Real* Ey_p,
#if (AMREX_SPACEDIM == 3)
                              amrex::Real* Ez_p,
#endif
                              const amrex::Real* Ex, const amrex::Real* Ey,
#if (AMREX_SPACEDIM == 3)
                              const amrex::Real* Ez,
#endif
                              const int* lo, const int* hi,
                              const amrex::Real* plo, const amrex::Real* dx,
                              const int* ng);

    void WRPX_PUSH_LEAPFROG(amrex::Real* particles, int ns, int np,
                            amrex::Real* vx_p, amrex::Real* vy_p,
#if (AMREX_SPACEDIM == 3)
                            amrex::Real* vz_p,
#endif
                            const amrex::Real* Ex_p, const amrex::Real* Ey_p,
#if (AMREX_SPACEDIM == 3)
                            const amrex::Real* Ez_p,
#endif
                            const amrex::Real* charge, const amrex::Real* mass, const amrex::Real* dt,
                            const amrex::Real* prob_lo, const amrex::Real* prob_hi);

    void WRPX_PUSH_LEAPFROG_POSITIONS(amrex::Real* particles, int ns, int np,
                                      amrex::Real* vx_p, amrex::Real* vy_p,
#if (AMREX_SPACEDIM == 3)
                                      amrex::Real* vz_p,
#endif
                                      const amrex::Real* dt, const amrex::Real* prob_lo,
                                      const amrex::Real* prob_hi);

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
                            const amrex::Real* dtsdx,
			    const amrex::Real* dtsdy,
			    const amrex::Real* dtsdz,
			    const int* maxwell_fdtd_solver_id);


    void WRPX_PUSH_PML_EVEC(const int* xlo, const int* xhi,
                            const int* ylo, const int* yhi,
                            const int* zlo, const int* zhi,
                            BL_FORT_FAB_ARG_3D(ex),
                            BL_FORT_FAB_ARG_3D(ey),
                            BL_FORT_FAB_ARG_3D(ez),
                            const BL_FORT_FAB_ARG_3D(bx),
                            const BL_FORT_FAB_ARG_3D(by),
                            const BL_FORT_FAB_ARG_3D(bz),
                            const amrex::Real* dtsdx,
                            const amrex::Real* dtsdy,
                            const amrex::Real* dtsdz);

    void WRPX_PUSH_PML_EVEC_F(const int* xlo, const int* xhi,
                              const int* ylo, const int* yhi,
                              const int* zlo, const int* zhi,
                              BL_FORT_FAB_ARG_3D(ex),
                              BL_FORT_FAB_ARG_3D(ey),
                              BL_FORT_FAB_ARG_3D(ez),
                              const BL_FORT_FAB_ARG_3D(f),
                              const amrex::Real* dtsdx,
                              const amrex::Real* dtsdy,
                              const amrex::Real* dtsdz,
                              const int* maxwell_fdtd_solver_id);

    void WRPX_PUSH_PML_F(const int* lo, const int* hi,
                         BL_FORT_FAB_ARG_3D(f),
                         const BL_FORT_FAB_ARG_3D(ex),
                         const BL_FORT_FAB_ARG_3D(ey),
                         const BL_FORT_FAB_ARG_3D(ez),
                         const amrex::Real* dtsdx,
                         const amrex::Real* dtsdy,
                         const amrex::Real* dtsdz);

    void WRPX_DAMP_PML (const int* texlo, const int* texhi,
                        const int* teylo, const int* teyhi,
                        const int* tezlo, const int* tezhi,
                        const int* tbxlo, const int* tbxhi,
                        const int* tbylo, const int* tbyhi,
                        const int* tbzlo, const int* tbzhi,
                        amrex::Real* ex, const int* exlo, const int* exhi,
                        amrex::Real* ey, const int* eylo, const int* eyhi,
                        amrex::Real* ez, const int* ezlo, const int* ezhi,
                        amrex::Real* bx, const int* bxlo, const int* bxhi,
                        amrex::Real* by, const int* bylo, const int* byhi,
                        amrex::Real* bz, const int* bzlo, const int* bzhi,
                        const amrex::Real* sigex, int sigex_lo, int sigex_hi,
#if (AMREX_SPACEDIM == 3)
                        const amrex::Real* sigey, int sigey_lo, int sigey_hi,
#endif
                        const amrex::Real* sigez, int sigez_lo, int sigez_hi,
                        const amrex::Real* sigbx, int sigbx_lo, int sigbx_hi,
#if (AMREX_SPACEDIM == 3)
                        const amrex::Real* sigby, int sigby_lo, int sigby_hi,
#endif
                        const amrex::Real* sigbz, int sigbz_lo, int sigbz_hi);

    void WRPX_DAMP_PML_F (const int* tndlo, const int* tndhi,
                          amrex::Real* F,  const int* flo, const int* fhi,
                          const amrex::Real* sigex, int sigex_lo, int sigex_hi,
#if (AMREX_SPACEDIM == 3)
                          const amrex::Real* sigey, int sigey_lo, int sigey_hi,
#endif
                          const amrex::Real* sigez, int sigez_lo, int sigez_hi,
                          const amrex::Real* sigbx, int sigbx_lo, int sigbx_hi,
#if (AMREX_SPACEDIM == 3)
                          const amrex::Real* sigby, int sigby_lo, int sigby_hi,
#endif
                          const amrex::Real* sigbz, int sigbz_lo, int sigbz_hi);


    void WRPX_SYNC_CURRENT (const int* lo, const int* hi,
                             BL_FORT_FAB_ARG_ANYD(crse),
                             const BL_FORT_FAB_ARG_ANYD(fine),
                             const int* dir);

    void WRPX_SYNC_RHO (const int* lo, const int* hi,
                        BL_FORT_FAB_ARG_ANYD(crse),
                        const BL_FORT_FAB_ARG_ANYD(fine),
                        const int* ncomp);

#ifdef WARPX_USE_PSATD
    void warpx_fft_mpi_init (int fcomm);
    void warpx_fft_domain_decomp (int* warpx_local_nz, int* warpx_local_z0,
                                  const int* global_lo, const int* global_hi);
    void warpx_fft_dataplan_init (const int* nox, const int* noy, const int* noz,
                                  void* fft_data, const int* ndata,
                                  const amrex_real* dx_w, const amrex_real* dt_w,
                                  const int* fftw_plan_measure, const int* do_nodal );
    void warpx_fft_nullify ();
    void warpx_fft_push_eb (amrex_real* ex_w, const int* exlo, const int* exhi,
                            amrex_real* ey_w, const int* eylo, const int* eyhi,
                            amrex_real* ez_w, const int* ezlo, const int* ezhi,
                            amrex_real* bx_w, const int* bxlo, const int* bxhi,
                            amrex_real* by_w, const int* bylo, const int* byhi,
                            amrex_real* bz_w, const int* bzlo, const int* bzhi,
                            amrex_real* jx_w, const int* jxlo, const int* jxhi,
                            amrex_real* jy_w, const int* jylo, const int* jyhi,
                            amrex_real* jz_w, const int* jzlo, const int* jzhi,
                            amrex_real* rhoold_w, const int* r1lo, const int* r1hi,
                            amrex_real* rho_w, const int* r2lo, const int* r2hi);

#endif

    void warpx_build_buffer_masks (const int* lo, const int* hi,
                                   int* msk, const int* mlo, const int* mhi,
                                   const int* gmsk, const int* glo, const int* ghi, const int* ng);

#ifdef __cplusplus
}
#endif
