#ifndef SPECIALIZE_H
#define SPECIALIZE_H

#include <REAL.H>
#include <ArrayLim.H>

#if defined(BL_LANG_FORT)

#define FORT_FASTCOPY        fastcopy
#define FORT_FASTCOPYTOMEM   fastcopytomem
#define FORT_FASTCOPYFROMMEM fastcopyfrommem
#define FORT_FASTSETVAL      fastsetval
#define FORT_FASTZERONORM    fastzeronorm
#define FORT_FASTONENORM     fastonenorm
#define FORT_FASTSUM         fastsum
#define FORT_FASTPLUS        fastplus
#define FORT_FASTMULT        fastmult
#define FORT_FASTSAXPY       fastsaxpy
#define FORT_FASTMINUS       fastminus
#define FORT_FASTDIVIDE      fastdivide
#define FORT_FASTPROTDIVIDE  fastprotdivide
#define FORT_FASTINVERT      fastinvert
#else

#if defined(BL_FORT_USE_UPPERCASE)
#define FORT_FASTCOPY        FASTCOPY
#define FORT_FASTCOPYTOMEM   FASTCOPYTOMEM
#define FORT_FASTCOPYFROMMEM FASTCOPYFROMMEM
#define FORT_FASTSETVAL      FASTSETVAL
#define FORT_FASTZERONORM    FASTZERONORM
#define FORT_FASTONENORM     FASTONENORM
#define FORT_FASTSUM         FASTSUM
#define FORT_FASTPLUS        FASTPLUS
#define FORT_FASTMULT        FASTMULT
#define FORT_FASTSAXPY       FASTSAXPY
#define FORT_FASTMINUS       FASTMINUS
#define FORT_FASTDIVIDE      FASTDIVIDE
#define FORT_FASTPROTDIVIDE  FASTPROTDIVIDE
#define FORT_FASTINVERT      FASTINVERT
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_FASTSAXPY       fastsaxpy
#define FORT_FASTCOPY        fastcopy
#define FORT_FASTCOPYTOMEM   fastcopytomem
#define FORT_FASTCOPYFROMMEM fastcopyfrommem
#define FORT_FASTSETVAL      fastsetval
#define FORT_FASTZERONORM    fastzeronorm
#define FORT_FASTONENORM     fastonenorm
#define FORT_FASTSUM         fastsum
#define FORT_FASTPLUS        fastplus
#define FORT_FASTMULT        fastmult
#define FORT_FASTMINUS       fastminus
#define FORT_FASTDIVIDE      fastdivide
#define FORT_FASTPROTDIVIDE  fastprotdivide
#define FORT_FASTINVERT      fastinvert
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_FASTCOPY        fastcopy_
#define FORT_FASTCOPYTOMEM   fastcopytomem_
#define FORT_FASTCOPYFROMMEM fastcopyfrommem_
#define FORT_FASTSETVAL      fastsetval_
#define FORT_FASTZERONORM    fastzeronorm_
#define FORT_FASTONENORM     fastonenorm_
#define FORT_FASTSUM         fastsum_
#define FORT_FASTPLUS        fastplus_
#define FORT_FASTMULT        fastmult_
#define FORT_FASTSAXPY       fastsaxpy_
#define FORT_FASTMINUS       fastminus_
#define FORT_FASTDIVIDE      fastdivide_
#define FORT_FASTPROTDIVIDE  fastprotdivide_
#define FORT_FASTINVERT      fastinvert_
#endif

extern "C"
{
    void FORT_FASTCOPY (Real* dest,
                        ARLIM_P(dlo), ARLIM_P(dhi),
                        D_DECL(const int&,const int&,const int&),
                        D_DECL(const int&,const int&,const int&),
                        const Real* src,
                        ARLIM_P(slo), ARLIM_P(shi),
                        D_DECL(const int&,const int&,const int&),
                        const int* ncomp);

    void FORT_FASTCOPYTOMEM (const int*  lo,
                             const int*  hi,
                             const Real* data,
                             ARLIM_P(dlo), ARLIM_P(dhi),
                             const int*  ncomp,
                             Real* dst);

    void FORT_FASTCOPYFROMMEM (const int*  lo,
                               const int*  hi,
                               Real* data,
                               ARLIM_P(dlo), ARLIM_P(dhi),
                               const int*  ncomp,
                               const Real* dst);

    void FORT_FASTSETVAL (const Real* val,
                          const int*  lo,
                          const int*  hi,
                          Real*       dest,
                          ARLIM_P(dlo), ARLIM_P(dhi),
                          const int*  ncomp);


    void FORT_FASTZERONORM (const Real* dest,
                            ARLIM_P(dlo), ARLIM_P(dhi),
                            const int*  lo,
                            const int*  hi,
                            const int*  ncomp,
                            const Real* nrm);

    void FORT_FASTONENORM (const Real* dest,
                           ARLIM_P(dlo), ARLIM_P(dhi),
                           const int*  lo,
                           const int*  hi,
                           const int*  ncomp,
                           const Real* nrm);

    void FORT_FASTSUM (const Real* dest,
                       ARLIM_P(dlo), ARLIM_P(dhi),
                       const int*  lo,
                       const int*  hi,
                       const int*  ncomp,
                       const Real* nrm);

    void FORT_FASTPLUS (Real* dest,
                        ARLIM_P(dlo), ARLIM_P(dhi),
                        D_DECL(const int&,const int&,const int&),
                        D_DECL(const int&,const int&,const int&),
                        const Real* src,
                        ARLIM_P(slo), ARLIM_P(shi),
                        D_DECL(const int&,const int&,const int&),
                        const int* ncomp);

    void FORT_FASTMULT (Real* dest,
                        ARLIM_P(dlo), ARLIM_P(dhi),
                        D_DECL(const int&,const int&,const int&),
                        D_DECL(const int&,const int&,const int&),
                        const Real* src,
                        ARLIM_P(slo), ARLIM_P(shi),
                        D_DECL(const int&,const int&,const int&),
                        const int* ncomp);

    void FORT_FASTSAXPY (Real* dest,
                         ARLIM_P(dlo), ARLIM_P(dhi),
                         D_DECL(const int&,const int&,const int&),
                         D_DECL(const int&,const int&,const int&),
                         const Real* a,
                         const Real* src,
                         ARLIM_P(slo), ARLIM_P(shi),
                         D_DECL(const int&,const int&,const int&),
                         const int* ncomp);

    void FORT_FASTMINUS (Real* dest,
                         ARLIM_P(dlo), ARLIM_P(dhi),
                         D_DECL(const int&,const int&,const int&),
                         D_DECL(const int&,const int&,const int&),
                         const Real* src,
                         ARLIM_P(slo), ARLIM_P(shi),
                         D_DECL(const int&,const int&,const int&),
                         const int* ncomp);

    void FORT_FASTDIVIDE (Real* dest,
                          ARLIM_P(dlo), ARLIM_P(dhi),
                          D_DECL(const int&,const int&,const int&),
                          D_DECL(const int&,const int&,const int&),
                          const Real* src,
                          ARLIM_P(slo), ARLIM_P(shi),
                          D_DECL(const int&,const int&,const int&),
                          const int* ncomp);

    void FORT_FASTPROTDIVIDE (Real* dest,
                              ARLIM_P(dlo), ARLIM_P(dhi),
                              D_DECL(const int&,const int&,const int&),
                              D_DECL(const int&,const int&,const int&),
                              const Real* src,
                              ARLIM_P(slo), ARLIM_P(shi),
                              D_DECL(const int&,const int&,const int&),
                              const int* ncomp);

    void FORT_FASTINVERT (const Real* dest,
                          ARLIM_P(dlo), ARLIM_P(dhi),
                          const int*  lo,
                          const int*  hi,
                          const Real* val,
                          const int*  ncomp);
}

#endif

#endif /*SPECIALIZE_H*/
