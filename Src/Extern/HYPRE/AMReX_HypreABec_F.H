#ifndef AMREX_HYPREABEC_F_H_
#define AMREX_HYPREABEC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_hpacoef (const int* lo, const int* hi, amrex_real* mat,
                        const amrex_real* acoefs, const int* alo, const int* ahi,
                        const amrex_real* sa);

    void amrex_hpbcoef (const int* lo, const int* hi, amrex_real* mat,
                        const amrex_real* bcoefs, const int* blo, const int* bhi,
                        const amrex_real* sb, const amrex_real* dx, const int* idim);

    void amrex_hpmat (const int* lo, const int* hi, amrex_real* mat,
                      const amrex_real* bcoefs, const int* blo, const int* bhi,
                      const int* mask, const int* mlo, const int* mhi,
                      const amrex_real* sb, const amrex_real* dx, const int* cdir,
                      const int* bct, const amrex_real* bcl);

    void amrex_hmac_ij(BL_FORT_FAB_ARG(acoefs),
                ARLIM_P(reglo), ARLIM_P(reghi),
                const amrex_real& alpha, void *A, BL_FORT_IFAB_ARG(Index));
    
    void amrex_hmbc_ij(BL_FORT_FAB_ARG(bcoefs),
                    ARLIM_P(reglo), ARLIM_P(reghi),
                    const amrex_real& beta, const amrex_real* dx, const int& n, 
                    void* A, BL_FORT_IFAB_ARG(Index));

    void amrex_hmmat_ij(ARLIM_P(reglo), ARLIM_P(reghi),
                     const int& cdir, const int& bct,
                     const int& bho, const amrex_real& bcl,
                     const BL_FORT_IFAB_ARG(mask),
                     BL_FORT_FAB_ARG(bcoefs),
                     const amrex_real& beta, const amrex_real* dx,
                     void* A, BL_FORT_IFAB_ARG(Index));

    void amrex_hmmat3_ij(ARLIM_P(reglo), ARLIM_P(reghi),
                      const int& cdir, const int& bctype,
                      const int& bho, const amrex_real& bcl,
                      const BL_FORT_IFAB_ARG(mask),
                      BL_FORT_FAB_ARG(bcoefs),
                      const amrex_real& beta, const amrex_real* dx,
                      void* A, BL_FORT_IFAB_ARG(Index));

    
    void amrex_BuildGlobalIndex(BL_FORT_IFAB_ARG(Index),
                                ARLIM_P(reglo), ARLIM_P(reghi),
                                const int& numCellsProc);
    
    void amrex_MatConvertGlobal(amrex_real* mat, int* rows, int* cols,
                                amrex_real* values, int* numCols,
                                const int& nCSR, const int& nRows,
                                ARLIM_P(reglo), ARLIM_P(reghi),
                                BL_FORT_IFAB_ARG(Index), const int& nDim);
    
    void amrex_conv_Vec_Local_Global(void* x, amrex_real* vec,
                                   const int& nRows,
                                   ARLIM_P(reglo), ARLIM_P(reghi),
                                   BL_FORT_IFAB_ARG(Index) );

    void amrex_conv_Vec_Global_Local(BL_FORT_FAB_ARG(vec), amrex_real* VecGB,
                                     const int& nRows,
                                     ARLIM_P(reglo), ARLIM_P(reghi));

    
#ifdef __cplusplus
};
#endif

#endif /* _HYPREABEC_F_H_ */
