#ifndef _HYPREABEC_F_H_
#define _HYPREABEC_F_H_

#include <AMReX_REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_HPACOEF          hpacoef_2d
#define FORT_HPBCOEF          hpbcoef_2d
#define FORT_HPBMAT3          hpbmat3_2d
#define FORT_HPBVEC3          hpbvec3_2d
#elif (BL_SPACEDIM == 3)
#define FORT_HPACOEF          hpacoef_3d
#define FORT_HPBCOEF          hpbcoef_3d
#define FORT_HPBMAT3          hpbmat3_3d
#define FORT_HPBVEC3          hpbvec3_3d
#endif

#else

#if (BL_SPACEDIM == 2)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_HPACOEF          HPACOEF_2d 
#define FORT_HPBCOEF          HPBCOEF_2d 
#define FORT_HPBMAT3          HPBMAT3_2d
#define FORT_HPBVEC3          HPBVEC3_2d
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_HPACOEF          hpacoef_2d 
#define FORT_HPBCOEF          hpbcoef_2d 
#define FORT_HPBMAT3          hpbmat3_2d
#define FORT_HPBVEC3          hpbvec3_2d
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_HPACOEF          hpacoef_2d_
#define FORT_HPBCOEF          hpbcoef_2d_
#define FORT_HPBMAT3          hpbmat3_2d_
#define FORT_HPBVEC3          hpbvec3_2d_
#endif
#endif

#if (BL_SPACEDIM == 3)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_HPACOEF          HPACOEF_3d 
#define FORT_HPBCOEF          HPBCOEF_3d 
#define FORT_HPBMAT3          HPBMAT3_3d
#define FORT_HPBVEC3          HPBVEC3_3d
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_HPACOEF          hpacoef_3d 
#define FORT_HPBCOEF          hpbcoef_3d 
#define FORT_HPBMAT3          hpbmat3_3d
#define FORT_HPBVEC3          hpbvec3_3d
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_HPACOEF          hpacoef_3d_
#define FORT_HPBCOEF          hpbcoef_3d_ 
#define FORT_HPBMAT3          hpbmat3_3d_
#define FORT_HPBVEC3          hpbvec3_3d_
#endif
#endif

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>
#include "HYPRE_parcsr_ls.h" 
#include "_hypre_parcsr_mv.h" 

#ifdef __cplusplus
extern "C" {
#endif
  void FORT_HPACOEF(amrex_real* mat, amrex_real* acoefs, ARLIM_P(alo), ARLIM_P(ahi),
		    const int* lo, const int* hi, const amrex_real& sa);
  void FORT_HPBCOEF(amrex_real* mat, amrex_real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const amrex_real& sb,
		    const amrex_real* dx, const int& idim);
  void FORT_HPBMAT3(amrex_real* mat, amrex_real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const amrex_real& sb, const amrex_real* dx, 
		    const int& idim, const int& bct, const amrex_real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi));
  void FORT_HPBVEC3(amrex_real* vec, amrex_real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const amrex_real& sb, const amrex_real* dx, 
		    const int& idim, const int& bct, const amrex_real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi),
		    const amrex_real* bcv, ARLIM_P(bvlo), ARLIM_P(bvhi));

    void amrex_hmac(amrex::Real* mat, 
                    BL_FORT_FAB_ARG(acoefs),
                    ARLIM_P(reglo), ARLIM_P(reghi),
                    const amrex::Real& alpha);
    
    void amrex_hmac_ij(BL_FORT_FAB_ARG(acoefs),
                ARLIM_P(reglo), ARLIM_P(reghi),
                const amrex::Real& alpha, HYPRE_IJMatrix &A, BL_FORT_IFAB_ARG(Index));
    
    void amrex_hmbc(amrex::Real* mat,
                    BL_FORT_FAB_ARG(bcoefs),
                    ARLIM_P(reglo), ARLIM_P(reghi),
                    const amrex::Real& beta, const amrex::Real* dx, const int& n);

    void amrex_hmbc_ij(BL_FORT_FAB_ARG(bcoefs),
                    ARLIM_P(reglo), ARLIM_P(reghi),
                    const amrex::Real& beta, const amrex::Real* dx, const int& n, 
                    HYPRE_IJMatrix &A, BL_FORT_IFAB_ARG(Index));
    
    void amrex_hmmat(amrex::Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
                     const int& cdir, const int& bct,
                     const int& bho, const amrex::Real& bcl,
                     const BL_FORT_IFAB_ARG(mask),
                     BL_FORT_FAB_ARG(bcoefs),
                     const amrex::Real& beta, const amrex::Real* dx);

    void amrex_hmmat_ij(ARLIM_P(reglo), ARLIM_P(reghi),
                     const int& cdir, const int& bct,
                     const int& bho, const amrex::Real& bcl,
                     const BL_FORT_IFAB_ARG(mask),
                     BL_FORT_FAB_ARG(bcoefs),
                     const amrex::Real& beta, const amrex::Real* dx,
                     HYPRE_IJMatrix &A, BL_FORT_IFAB_ARG(Index));

    void amrex_hmmat3(amrex::Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
                      const int& cdir, const int& bctype,
                      const int& bho, const amrex::Real& bcl,
                      const BL_FORT_IFAB_ARG(mask),
                      BL_FORT_FAB_ARG(bcoefs),
                      const amrex::Real& beta, const amrex::Real* dx);
    
    void amrex_hmmat3_ij(ARLIM_P(reglo), ARLIM_P(reghi),
                      const int& cdir, const int& bctype,
                      const int& bho, const amrex::Real& bcl,
                      const BL_FORT_IFAB_ARG(mask),
                      BL_FORT_FAB_ARG(bcoefs),
                      const amrex::Real& beta, const amrex::Real* dx,
                      HYPRE_IJMatrix &A, BL_FORT_IFAB_ARG(Index));

    
    void amrex_hbvec (amrex::Real* vec, ARLIM_P(reglo), ARLIM_P(reghi),
                const int& cdir, const int& bctype,
                const int& bho, const amrex::Real& bcl,
                const amrex::Real* fs, ARLIM_P(fslo), ARLIM_P(fshi),
                const BL_FORT_IFAB_ARG(mask),
                const BL_FORT_FAB_ARG(bcoefs),
                const amrex::Real& beta, const amrex::Real* dx);

    void amrex_hbvec3(amrex::Real* vec, ARLIM_P(reglo), ARLIM_P(reghi),
                const int& cdir, const int& bctype,
                const int& bho, const amrex::Real& bcl,
                const amrex::Real* fs, ARLIM_P(fslo), ARLIM_P(fshi),
                const BL_FORT_IFAB_ARG(mask),
                const BL_FORT_FAB_ARG(bcoefs),
                const amrex::Real& beta, const amrex::Real* dx);

    void amrex_BuildGlobalIndex(BL_FORT_IFAB_ARG(Index),ARLIM_P(reglo), ARLIM_P(reghi),
                                const int& numCellsProc);
    
    void amrex_MatConvertGlobal(amrex::Real* mat, int* rows, int* cols, amrex::Real* values, int* numCols,
                                const int& nCSR, const int& nRows,
                                ARLIM_P(reglo), ARLIM_P(reghi),
                                BL_FORT_IFAB_ARG(Index), const int& nDim);
    
    void amrex_conv_Vec_Local_Global(amrex::Real* vec,amrex::Real* VecGB,int* VecIndices,const int& nRows,
                                     ARLIM_P(reglo), ARLIM_P(reghi),
                                     BL_FORT_IFAB_ARG(Index) );

    void amrex_conv_Vec_Global_Local(BL_FORT_FAB_ARG(vec), amrex::Real* VecGB, const int& nRows,
                                     ARLIM_P(reglo), ARLIM_P(reghi));

    
#ifdef __cplusplus
};
#endif

#endif

#endif /* _HYPREABEC_F_H_ */
