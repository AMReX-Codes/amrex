#ifdef BL_LANG_FORT
      integer dbunit
      logical dbnorm, dbmat, dbon
      common /hgdbg/ dbunit, dbnorm, dbmat, dbon
      save hgdbg
#else
#ifndef _HG_DEBUG_F_H_
#define _HG_DEBUG_F_H_
      
#if defined( BL_FORT_USE_UNDERSCORE )
#define FORT_HGDEBUGINIT  hgdebuginit_
#define FORT_HGDEBUGCTL   hgdebugctl_
#define FORT_HGDEBUGON    hgdebugon_
#define FORT_HGDEBUGOFF   hgdebugoff_
#elif defined( BL_FORT_USE_UPPERCASE )
#define FORT_HGDEBUGINIT  HGDEBUGINIT
#define FORT_HGDEBUGCTL   HGDEBUGCTL
#define FORT_HGDEBUGON    HGDEBUGON
#define FORT_HGDEBUGOFF   HGDEBUGOFF
#elif defined( BL_FORT_USE_LOWERCASE )
#define FORT_HGDEBUGINIT  hgdebuginit
#define FORT_HGDEBUGCTL   hgdebugctl
#define FORT_HGDEBUGON    hgdebugon
#define FORT_HGDEBUGOFF   hgdebugoff
#else
#error "none of BL_FORT_USE_{UNDERSCORE,UPPERCASE,LOWERCASE} defined"
#endif

extern "C" 
{
    void FORT_HGDEBUGINIT (const int*, const int*);
    void FORT_HGDEBUGCTL (const int*, const int*, const int*);
    void FORT_HGDEBUGON ();
    void FORT_HGDEBUGOFF ();
}

#endif
#endif      
