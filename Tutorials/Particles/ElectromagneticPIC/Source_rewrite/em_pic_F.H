#ifndef _EM_PIC_F_H_
#define _EM_PIC_F_H_

#include <AMReX_BLFort.H>

// A macro that add either _acc or _omp to function names
#if defined (AMREX_USE_ACC) && ! defined (AMREX_OMP_OFFLOAD)
#define FUNC_SUFFIX(x) x##_acc
#elif !defined (AMREX_USE_ACC) && defined (AMREX_OMP_OFFLOAD)
#define FUNC_SUFFIX(x) x##_omp
#else
#define FUNC_SUFFIX(x) x
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void FUNC_SUFFIX(push_momentum_boris)
        (const int np,
         amrex::Real* uxp,
         amrex::Real* uyp,
         amrex::Real* uzp,
         amrex::Real* gaminv,
         const amrex::Real* exp,
         const amrex::Real* eyp,
         const amrex::Real* ezp,
         const amrex::Real* bxp,
         const amrex::Real* byp,
         const amrex::Real* bzp,
         const amrex::Real q,
         const amrex::Real m,
         const amrex::Real dt);

    void FUNC_SUFFIX(push_position_boris)
        (const int np, void* structs,
         const amrex::Real* uxp,
         const amrex::Real* uyp,
         const amrex::Real* uzp,
         const amrex::Real* gaminv,
         const amrex::Real dt);

    void FUNC_SUFFIX(set_gamma)(const int np,
                         const amrex::Real* uxp,
                         const amrex::Real* uyp,
                         const amrex::Real* uzp,
                         amrex::Real* gaminv);

    void FUNC_SUFFIX(deposit_current)(BL_FORT_FAB_ARG_3D(jx),
                               BL_FORT_FAB_ARG_3D(jy),
                               BL_FORT_FAB_ARG_3D(jz),
                               const int np, const void* structs,
                               const amrex::Real* uxp,
                               const amrex::Real* uyp,
                               const amrex::Real* uzp,
                               const amrex::Real* gaminv,
                               const amrex::Real* w,
                               const amrex::Real q,
                               const amrex::Real* plo,
                               const amrex::Real dt,
                               const amrex::Real* dx);

    void FUNC_SUFFIX(gather_magnetic_field)(const int np, const void* structs,
                                     amrex::Real* bx,
                                     amrex::Real* by,
                                     amrex::Real* bz,
                                     const BL_FORT_FAB_ARG_3D(bxg),
                                     const BL_FORT_FAB_ARG_3D(byg),
                                     const BL_FORT_FAB_ARG_3D(bzg),
                                     const amrex::Real* plo, const amrex::Real* dx);


    void FUNC_SUFFIX(gather_electric_field)(const int np, const void* structs,
                                     amrex::Real* ex,
                                     amrex::Real* ey,
                                     amrex::Real* ez,
                                     const BL_FORT_FAB_ARG_3D(exg),
                                     const BL_FORT_FAB_ARG_3D(eyg),
                                     const BL_FORT_FAB_ARG_3D(ezg),
                                     const amrex::Real* plo, const amrex::Real* dx);

    void FUNC_SUFFIX(push_electric_field_x)(const int* xlo, const int* xhi,
                                     BL_FORT_FAB_ARG_3D(ex),
                                     const BL_FORT_FAB_ARG_3D(by),
                                     const BL_FORT_FAB_ARG_3D(bz),
                                     const BL_FORT_FAB_ARG_3D(jx),
                                     const amrex::Real mudt,
                                     const amrex::Real dtsdy,
                                     const amrex::Real dtsdz);

    void FUNC_SUFFIX(push_electric_field_y)(const int* ylo, const int* yhi,
                                     BL_FORT_FAB_ARG_3D(ey),
                                     const BL_FORT_FAB_ARG_3D(bx),
                                     const BL_FORT_FAB_ARG_3D(bz),
                                     const BL_FORT_FAB_ARG_3D(jy),
                                     const amrex::Real mudt,
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdz);

    void FUNC_SUFFIX(push_electric_field_z)(const int* zlo, const int* zhi,
                                     BL_FORT_FAB_ARG_3D(ez),
                                     const BL_FORT_FAB_ARG_3D(bx),
                                     const BL_FORT_FAB_ARG_3D(by),
                                     const BL_FORT_FAB_ARG_3D(jz),
                                     const amrex::Real mudt,
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdy);

    void FUNC_SUFFIX(push_magnetic_field_x)(const int* xlo, const int* xhi,
                                     BL_FORT_FAB_ARG_3D(bx),
                                     const BL_FORT_FAB_ARG_3D(ey),
                                     const BL_FORT_FAB_ARG_3D(ez),
                                     const amrex::Real dtsdy,
                                     const amrex::Real dtsdz);

    void FUNC_SUFFIX(push_magnetic_field_y)(const int* ylo, const int* yhi,
                                     BL_FORT_FAB_ARG_3D(by),
                                     const BL_FORT_FAB_ARG_3D(ex),
                                     const BL_FORT_FAB_ARG_3D(ez),
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdz);

    void FUNC_SUFFIX(push_magnetic_field_z)(const int* zlo, const int* zhi,
                                     BL_FORT_FAB_ARG_3D(bz),
                                     const BL_FORT_FAB_ARG_3D(ex),
                                     const BL_FORT_FAB_ARG_3D(ey),
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdy);

    void check_langmuir_solution(const int* boxlo,  const int* boxhi,
                                 const int* testlo, const int* testhi,
                                 const BL_FORT_FAB_ARG_3D(jx), const amrex::Real time,
                                 amrex::Real* max_error);


#ifdef __cplusplus
};
#endif

#endif /*_EM_PIC_F_H_*/
