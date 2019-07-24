#ifndef AMREX_EB2_F_H_
#define AMREX_EB2_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_eb2_coarsen_from_fine (const int* lo, const int* hi,
                                      const int* xlo, const int* xhi,
                                      const int* ylo, const int* yhi,
#if (AMREX_SPACEDIM == 3)
                                      const int* zlo, const int* zhi,
#endif
                                      amrex_real* cvol, const int* cvlo, const int* cvhi,
                                      const amrex_real* fvol, const int* fvlo, const int* fvhi,
                                      amrex_real* ccent, const int* cclo, const int* cchi,
                                      const amrex_real* fcent, const int* fclo, const int* fchi,
                                      amrex_real* cba, const int* cbalo, const int* cbahi,
                                      const amrex_real* fba, const int* fbalo, const int* fbahi,
                                      amrex_real* cbc, const int* cbclo, const int* cbchi,
                                      const amrex_real* fbc, const int* fbclo, const int* fbchi,
                                      amrex_real* cbn, const int* cbnlo, const int* cbnhi,
                                      const amrex_real* fbn, const int* fbnlo, const int* fbnhi,
                                      amrex_real* capx, const int* caxlo, const int* caxhi,
                                      const amrex_real* fapx, const int* faxlo, const int* faxhi,
                                      amrex_real* capy, const int* caylo, const int* cayhi,
                                      const amrex_real* fapy, const int* faylo, const int* fayhi,
#if (AMREX_SPACEDIM == 3)
                                      amrex_real* capz, const int* cazlo, const int* cazhi,
                                      const amrex_real* fapz, const int* fazlo, const int* fazhi,
#endif
                                      amrex_real* cfcx, const int* cfxlo, const int* cfxhi,
                                      const amrex_real* ffcx, const int* ffxlo, const int* ffxhi,
                                      amrex_real* cfcy, const int* cfylo, const int* cfyhi,
                                      const amrex_real* ffcy, const int* ffylo, const int* ffyhi,
#if (AMREX_SPACEDIM == 3)
                                      amrex_real* cfcz, const int* cfzlo, const int* cfzhi,
                                      const amrex_real* ffcz, const int* ffzlo, const int* ffzhi,
#endif
                                      const void* cflag, const int* cflo, const int* cfhi,
                                      const void* fflag, const int* fflo, const int* ffhi,
                                      int* tile_error);


    void amrex_eb2_build_cellflag_from_ap (const int* lo, const int*hi,
                                           const amrex_real* apx, const int* axlo, const int* axhi,
                                           const amrex_real* apy, const int* aylo, const int* ayhi,
#if (AMREX_SPACEDIM == 3)
                                           const amrex_real* apz, const int* azlo, const int* azhi,
#endif
                                           void* flag, const int* flo, const int* fhi);
    

    void amrex_eb2_check_mvmc (const int* cclo, const int* cchi,
                               const int* ndlo, const int* ndhi,
                               amrex_real* cls, const int* clo, const int* chi,
                               const amrex_real* fls, const int* flo, const int* fhi,
#if (AMREX_SPACEDIM == 3)
                               int* ncuts, const int* tlo, const int* thi,
#endif
                               int* error);

#ifdef __cplusplus
}
#endif

#endif
