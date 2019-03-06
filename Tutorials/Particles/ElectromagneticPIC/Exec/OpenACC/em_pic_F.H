#ifndef _EM_PIC_F_H_
#define _EM_PIC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void push_momentum_boris
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

    void push_position_boris
        (const int np, void* structs,
         const amrex::Real* uxp,
         const amrex::Real* uyp,
         const amrex::Real* uzp,
         const amrex::Real* gaminv,
         const amrex::Real dt);

    void set_gamma(const int np,
                         const amrex::Real* uxp,
                         const amrex::Real* uyp,
                         const amrex::Real* uzp,
                         amrex::Real* gaminv);

    void deposit_current(BL_FORT_FAB_ARG_3D(jx),
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

    void gather_magnetic_field(const int np, const void* structs,
                                     amrex::Real* bx,
                                     amrex::Real* by,
                                     amrex::Real* bz,
                                     const BL_FORT_FAB_ARG_3D(bxg),
                                     const BL_FORT_FAB_ARG_3D(byg),
                                     const BL_FORT_FAB_ARG_3D(bzg),
                                     const amrex::Real* plo, const amrex::Real* dx);


    void gather_electric_field(const int np, const void* structs,
                                     amrex::Real* ex,
                                     amrex::Real* ey,
                                     amrex::Real* ez,
                                     const BL_FORT_FAB_ARG_3D(exg),
                                     const BL_FORT_FAB_ARG_3D(eyg),
                                     const BL_FORT_FAB_ARG_3D(ezg),
                                     const amrex::Real* plo, const amrex::Real* dx);

    void push_electric_field_x(const int* xlo, const int* xhi,
                                     BL_FORT_FAB_ARG_3D(ex),
                                     const BL_FORT_FAB_ARG_3D(by),
                                     const BL_FORT_FAB_ARG_3D(bz),
                                     const BL_FORT_FAB_ARG_3D(jx),
                                     const amrex::Real mudt,
                                     const amrex::Real dtsdy,
                                     const amrex::Real dtsdz);

    void push_electric_field_y(const int* ylo, const int* yhi,
                                     BL_FORT_FAB_ARG_3D(ey),
                                     const BL_FORT_FAB_ARG_3D(bx),
                                     const BL_FORT_FAB_ARG_3D(bz),
                                     const BL_FORT_FAB_ARG_3D(jy),
                                     const amrex::Real mudt,
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdz);

    void push_electric_field_z(const int* zlo, const int* zhi,
                                     BL_FORT_FAB_ARG_3D(ez),
                                     const BL_FORT_FAB_ARG_3D(bx),
                                     const BL_FORT_FAB_ARG_3D(by),
                                     const BL_FORT_FAB_ARG_3D(jz),
                                     const amrex::Real mudt,
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdy);

    void push_magnetic_field_x(const int* xlo, const int* xhi,
                                     BL_FORT_FAB_ARG_3D(bx),
                                     const BL_FORT_FAB_ARG_3D(ey),
                                     const BL_FORT_FAB_ARG_3D(ez),
                                     const amrex::Real dtsdy,
                                     const amrex::Real dtsdz);

    void push_magnetic_field_y(const int* ylo, const int* yhi,
                                     BL_FORT_FAB_ARG_3D(by),
                                     const BL_FORT_FAB_ARG_3D(ex),
                                     const BL_FORT_FAB_ARG_3D(ez),
                                     const amrex::Real dtsdx,
                                     const amrex::Real dtsdz);

    void push_magnetic_field_z(const int* zlo, const int* zhi,
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
