#ifndef _MYNORM_H_
#define _MYNORM_H_
Real MFNorm (const MultiFab& mfab,
             const BoxArray& normBoxes,
             const int       exponent,
             const int       srcComp,
             const int       numComp);
#endif
