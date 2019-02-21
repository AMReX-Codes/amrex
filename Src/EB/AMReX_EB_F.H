#ifndef AMREX_EB_F_H_
#define AMREX_EB_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_eb_compute_normals(
            const int* lo , const int* hi,
            const void* flag,   const int* fglo, const int* fghi,
            amrex_real* normal, const int* nlo,  const int* nhi,
            const amrex_real* ax, const int* axlo, const int* axhi,
            const amrex_real* ay, const int* aylo, const int* ayhi,
            const amrex_real* az, const int* azlo, const int* azhi
        );

    void amrex_eb_init_levelset(
            const int * lo,   const int * hi,
            amrex_real * phi, const int * phlo, const int * phhi
        );

    void amrex_eb_fill_levelset(
            const int * lo,             const int * hi,
            const amrex_real * eb_list, const int * l_eb,
            int * valid,                const int * vlo,  const int * vhi,
            amrex_real * phi,           const int * phlo, const int * phhi,
            const amrex_real * dx,      const amrex_real * dx_eb
        );

    void amrex_eb_fill_levelset_loc(
            const int * lo,              const int * hi,
            const amrex_real * eb_list,  const int * l_eb,
            int * valid,                 const int * v_lo,      const int * v_hi,
            amrex_real * phi,            const int * p_lo,      const int * p_hi,
            const amrex_real * ls_guess, const int * lsg_lo,    const int * lsg_hi,
            const amrex_real * ls_thres, const amrex_real * dx, const amrex_real * dx_eb
        );

    void amrex_eb_fill_levelset_bcs(
            amrex_real * phi, const int * philo, const int * phihi,
            int * valid,      const int * vlo,   const int * vhi,
            const int * peri, const int * domlo, const int * domhi,
            const amrex_real * eb_list, const int * l_eb,
            const amrex_real * dx, const amrex_real * dx_eb
        );

    void amrex_eb_update_levelset_intersection(
            const int * lo,           const int * hi,
            const int * valid_in,     const int * vilo, const int * vihi,
            const amrex_real * ls_in, const int * lslo, const int * lshi,
            int * valid,              const int * vlo,  const int * vhi,
            amrex_real * phi,         const int * phlo, const int * phhi
        );

    void amrex_eb_update_levelset_intersection_bcs(
            const int * valid_in,     const int * vilo,  const int * vihi,
            const amrex_real * ls_in, const int * lslo,  const int * lshi,
            int * valid,              const int * vlo,   const int * vhi,
            amrex_real * phi,         const int * phlo,  const int * phhi,
            const int * periodic,     const int * domlo, const int * domhi
        );

    void amrex_eb_update_levelset_union(
            const int * lo,           const int * hi,
            const int * valid_in,     const int * vilo, const int * vihi,
            const amrex_real * ls_in, const int * lslo, const int * lshi,
            int * valid,              const int * vlo,  const int * vhi,
            amrex_real * phi,         const int * phlo, const int * phhi
        );

    void amrex_eb_update_levelset_union_bcs(
            const int * valid_in,     const int * vilo,  const int * vihi,
            const amrex_real * ls_in, const int * lslo,  const int * lshi,
            int * valid,              const int * vlo,   const int * vhi,
            amrex_real * phi,         const int * phlo,  const int * phhi,
            const int * periodic,     const int * domlo, const int * domhi
        );

    void amrex_eb_fill_valid(
            const int * lo,         const int * hi,
            int * valid,            const int * vlo,  const int * vhi,
            const amrex_real * phi, const int * phlo, const int * phhi,
            const int * n_pad
        );

    void amrex_eb_fill_valid_bcs(
            int * valid,          const int * vlo,   const int * vhi,
            const int * periodic, const int * domlo, const int * domhi
        );

    void amrex_eb_count_facets(
            const int * lo,    const int * hi,
            const void * flag, const int * flo, const int * fhi,
            int * n_facets
        );

    void amrex_eb_as_list(
            const int * lo,           const int * hi,  int * c_facets,
            const void * flag,        const int * flo, const int * fhi,
            const amrex_real * norm,  const int * nlo, const int * nhi,
            const amrex_real * bcent, const int * blo, const int * bhi,
            amrex_real * list_out,    const int * lsize,
            const amrex_real * dx
        );

    void amrex_eb_threshold_levelset(
            const int * lo,   const int * hi,     const amrex_real * threshold,
            amrex_real * phi, const int * phi_lo, const int * phi_hi
        );

    void amrex_eb_validate_levelset(
            const int * lo,          const int * hi,   const int * n_pad,
            const amrex_real * impf, const int * imlo, const int * imhi,
            const int * valid,       const int * vlo,  const int * vhi,
            amrex_real * phi,        const int * phlo, const int * phhi
        );

    void amrex_eb_validate_levelset_bcs(
            amrex_real * phi,        const int * phlo,  const int * phhi,
            const int * valid,       const int * vlo,   const int * vhi,
            const int * periodic,    const int * domlo, const int * domhi,
            const amrex_real * impf, const int * imlo,  const int * imhi
        );

#ifdef __cplusplus
}
#endif

#endif
