/*
** $Id: SPECIALIZE_F.H,v 1.2 1998-10-07 21:18:28 vince Exp $
*/

#ifndef SPECIALIZE_H
#define SPECIALIZE_H

#include <REAL.H>
#include <ArrayLim.H>

#ifdef BL_LANG_FORT

#define FORT_FASTCOPY   fastcopy
#define FORT_FASTSETVAL fastsetval

#endif /*BL_LANG_FORT*/


#ifdef BL_LANG_CC

#ifdef BL_FORT_USE_UPPERCASE
#define FORT_FASTCOPY       FASTCOPY
#define FORT_FASTSETVAL     FASTSETVAL
#else
#  ifdef BL_AIX
#define FORT_FASTCOPY       fastcopy
#define FORT_FASTSETVAL     fastsetval
#else
#define FORT_FASTCOPY       fastcopy_
#define FORT_FASTSETVAL     fastsetval_
#endif
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
                        D_DECL(const int&,const int&,const int&),
                        const int& ncomp);

    void FORT_FASTSETVAL (const Real* val,
                          const int*  lo,
                          const int*  hi,
                          Real*       dest,
                          ARLIM_P(dlo), ARLIM_P(dhi),
                          const int&  ncomp);
};
#endif /*BL_LANG_CC*/

#endif /*SPECIALIZE_H*/
