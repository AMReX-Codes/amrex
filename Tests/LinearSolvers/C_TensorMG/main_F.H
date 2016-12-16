#ifndef _main_F_H_
#define _main_F_H_ 

#include <AMReX_ArrayLim.H>

#if (BL_SPACEDIM == 2)
#define FORT_FILLCRSE  fillcrse2d
#define FORT_FILLFINE  fillfine2d
#define FORT_MAKEMU    makemu2d
#define FORT_FILLRHS    fillrhs2d
#endif

#if (BL_SPACEDIM == 3)
#define FORT_FILLCRSE  fillcrse3d
#define FORT_FILLFINE  fillfine3d
#define FORT_MAKEMU    makemu3d
#define FORT_FILLRHS    fillrhs3d
#endif

#if !defined(BL_LANG_FORT)
extern "C"
{
    void FORT_FILLCRSE (
	Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
	const Real* h, const int* nc
	);
    void FORT_FILLFINE (
	Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
	const Real* h, const int* nc
	);
    void FORT_FILLRHS (
	Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
	const Real* h, const int* nc
	);
    void FORT_MAKEMU(
		     Real*mu, ARLIM_P(mulo),ARLIM_P(muhi),
		     const Real* h, int&idir );
};
#endif

#endif
