#ifndef AMREX_INTERP_FACE_REG_C_H_
#define AMREX_INTERP_FACE_REG_C_H_

#include <AMReX_Array4.H>

#if (AMREX_SPACEDIM == 2)
#include <AMReX_InterpFaceReg_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_InterpFaceReg_3D_C.H>
#endif

#endif
