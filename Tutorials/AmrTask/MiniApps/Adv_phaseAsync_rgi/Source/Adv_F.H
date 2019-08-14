#ifndef _Adv_F_H_
#define _Adv_F_H_
#include <AMReX_BLFort.H>

extern "C" 
{
  void get_tagging_params(const int* name, const int* namlen);

  void nullfill(BL_FORT_FAB_ARG_3D(state),
		const int* dlo, const int* dhi,
		const amrex::Real* dx, const amrex::Real* glo, 
		const amrex::Real* time, const int* bc);
  
  void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		   BL_FORT_FAB_ARG_3D(state),
		   const int* tagval, const int* clearval,
		   const int* lo, const int* hi,
		   const amrex::Real* dx, const amrex::Real* problo,
		   const amrex::Real* time, const int* level);
  
  void initdata(const int& level, const amrex::Real& time, 
		const int* lo, const int* hi,
		BL_FORT_FAB_ARG_3D(state),
		const amrex::Real* dx, const amrex::Real* problo);
  
  void get_face_velocity(const int& level, const amrex::Real& time, 
			 D_DECL(BL_FORT_FAB_ARG(xvel),
				BL_FORT_FAB_ARG(yvel),
				BL_FORT_FAB_ARG(zvel)),
			 const amrex::Real* dx, const amrex::Real* problo);
  
  void advect(const amrex::Real& time, const int* lo, const int*hi,
	      const BL_FORT_FAB_ARG_3D(statein),
	      BL_FORT_FAB_ARG_3D(stateout),
	      D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
		     const BL_FORT_FAB_ARG_3D(yvel),
		     const BL_FORT_FAB_ARG_3D(zvel)),
	      D_DECL(BL_FORT_FAB_ARG_3D(fx),
		     BL_FORT_FAB_ARG_3D(fy),
		     BL_FORT_FAB_ARG_3D(fz)),
	      const amrex::Real* dx, const amrex::Real& dt);
}

#endif
