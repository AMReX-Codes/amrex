#ifndef CNS_F_H_
#define CNS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void cns_initdata (const int* level, const amrex_real* time,
                       const int* lo, const int* hi,
                       amrex_real* u, const int* ulo, const int* uhi,
                       const amrex_real* dx, const amrex_real* problo);

    void cns_init_fort (const int* phys_bc_lo, const int* phys_bc_hi,
                        int, int, int, int, int, int, int,
                        const amrex_real* problo, const amrex_real* probhi);

    void cns_hypfill (BL_FORT_FAB_ARG_3D(state),
                      const int* dlo, const int* dhi,
                      const amrex_real* dx, const amrex_real* glo, 
                      const amrex_real* time, const int* bc);
    
    void cns_denfill (BL_FORT_FAB_ARG_3D(state),
                      const int* dlo, const int* dhi,
                      const amrex_real* dx, const amrex_real* glo, 
                      const amrex_real* time, const int* bc);

    void cns_nullfill (BL_FORT_FAB_ARG_3D(state),
                       const int* dlo, const int* dhi,
                       const amrex_real* dx, const amrex_real* glo, 
                       const amrex_real* time, const int* bc);

    void cns_derpres (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                      const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const amrex_real* delta, const amrex_real* xlo,
                      const amrex_real* time, const amrex_real* dt, const int* bcrec, 
                      const int* level, const int* grid_no);

    void cns_dervel (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const amrex_real* delta, const amrex_real* xlo,
                     const amrex_real* time, const amrex_real* dt, const int* bcrec, 
                     const int* level, const int* grid_no);

    void cns_compute_temperature (const int* lo, const int* hi,
                                  amrex_real* u, const int* ulo, const int* uhi);

    void cns_estdt (const int* lo, const int* hi,
                    const amrex_real* u, const int* ulo, const int* uhi,
                    const amrex_real* dx, amrex_real* dt);
                    

    void cns_compute_dudt (const int* lo, const int* hi,
                           amrex_real* dudt, const int* utlo, const int* uthi,
                           const amrex_real* u, const int* ulo, const int* uhi,
                           amrex_real* fx, const int* fxlo, const int* fxhi,
                           amrex_real* fy, const int* fylo, const int* fyhi,
                           amrex_real* fz, const int* fzlo, const int* fzhi,
                           const amrex_real* dx, const amrex_real* dt, const int* level);

    void cns_eb_compute_dudt (const int* lo, const int* hi,
                              amrex_real* dudt, const int* utlo, const int* uthi,
                              const amrex_real* u, const int* ulo, const int* uhi,
                              amrex_real* fx, const int* fxlo, const int* fxhi,
                              amrex_real* fy, const int* fylo, const int* fyhi,
                              amrex_real* fz, const int* fzlo, const int* fzhi,
                              const void* flag, const int* fglo, const int* fghi,
                              const amrex_real* volfrac, const int* vlo, const int* vhi,
                              const amrex_real* bcent, const int* blo, const int* bhi,
                              const amrex_real* ax, const int* axlo, const int* axhi,
                              const amrex_real* ay, const int* aylo, const int* ayhi,
                              const amrex_real* az, const int* azlo, const int* azhi,
                              const amrex_real* cx, const int* cxlo, const int* cxhi,
                              const amrex_real* cy, const int* cylo, const int* cyhi,
                              const amrex_real* cz, const int* czlo, const int* czhi,
                              const int* as_crse,
                              amrex_real* rr_drho_crse, const int* rdclo, const int* rdchi,
                              const int* rr_flag_crse, const int* rfclo, const int* rfchi,
                              const int* as_fine,
                              amrex_real* dm_as_fine, const int* dflo, const int* dfhi,
                              const int* levmsk, const int* lmlo, const int* lmhi,
                              const amrex_real* dx, const amrex_real* dt, const int* level);

    void cns_tag_denerror (const int* lo, const int* hi,
                           char* tag, const int* tlo, const int* thi,
                           amrex_real* rho, const int* rlo, const int* rhi,
                           const void* flag, const int* flo, const int* fhi,
                           const amrex_real* dengrad, const char* tagval, const char* clearval);


#ifdef __cplusplus
}
#endif

#endif
