#ifndef _SMC_F_H_
#define _SMC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void derivative_stencil_init();
    void chemistry_init();
    void variables_init();

    int get_num_ghost_cells();
    int get_num_species();
    int get_num_cons();
    int get_num_prim();

    void chemterm_3d(const int* lo, const int* hi, 
		     const BL_FORT_FAB_ARG_3D(q), BL_FORT_FAB_ARG_3D(uprime)); 
    void comp_courno_3d(const int* lo, const int* hi, const amrex_real* dx, 
			const BL_FORT_FAB_ARG_3D(q), amrex_real& courno);
    void ctoprim_3d(const int* tlo, const int* thi, const int* lo, const int* hi,
		    amrex_real* u, amrex_real* q, const int& ngu, const int& ngq);
    void get_trans_prop_3d(const int* tlo, const int*thi, const int* lo, const int* hi,
			   const amrex_real* q, amrex_real* mu, amrex_real* xi, amrex_real* lam, amrex_real* Ddiag, int& ng);
    void init_data_3d(const int* tlo, const int* thi, const int* lo, const int* hi, 
		      amrex_real* u, const int& ng, const amrex_real dx[], const amrex_real plo[], const amrex_real phi[]);
    void make_plot_3d(const int* lo, const int* hi,
		      BL_FORT_FAB_ARG_3D(pv), const BL_FORT_FAB_ARG_3D(u));
    void reset_rho_3d(const int* tlo, const int* thi, const int* lo, const int* hi, 
		      amrex_real* u, const int& ng);
    void hypterm_3d(const int* tlo, const int*thi, const amrex_real* dx,
		    const BL_FORT_FAB_ARG_3D(u),
		    const BL_FORT_FAB_ARG_3D(q),
		    BL_FORT_FAB_ARG_3D(uprime));
    void narrow_diffterm_3d(const int* tlo, const int*thi, const amrex_real* dx,
			    const BL_FORT_FAB_ARG_3D(q),
			    BL_FORT_FAB_ARG_3D(uprime),
			    const amrex_real* mu, const amrex_real* xi, const amrex_real* lam, const amrex_real* Ddiag);
#ifdef __cplusplus
}
#endif

#endif

