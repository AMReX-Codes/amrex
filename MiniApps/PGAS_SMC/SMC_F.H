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
    void comp_courno_3d(const int* lo, const int* hi, const amrex::Real* dx, 
			const BL_FORT_FAB_ARG_3D(q), amrex::Real& courno);
    void ctoprim_3d(const int* tlo, const int* thi, const int* lo, const int* hi,
		    amrex::Real* u, amrex::Real* q, const int& ngu, const int& ngq);
    void get_trans_prop_3d(const int* tlo, const int*thi, const int* lo, const int* hi,
			   const amrex::Real* q, amrex::Real* mu, amrex::Real* xi, amrex::Real* lam, amrex::Real* Ddiag, int& ng);
    void init_data_3d(const int* tlo, const int* thi, const int* lo, const int* hi, 
		      amrex::Real* u, const int& ng, const amrex::Real dx[], const amrex::Real plo[], const amrex::Real phi[]);
    void make_plot_3d(const int* lo, const int* hi,
		      BL_FORT_FAB_ARG_3D(pv), const BL_FORT_FAB_ARG_3D(u));
    void reset_rho_3d(const int* tlo, const int* thi, const int* lo, const int* hi, 
		      amrex::Real* u, const int& ng);
    void hypterm_3d(const int* tlo, const int*thi, const amrex::Real* dx,
		    const BL_FORT_FAB_ARG_3D(u),
		    const BL_FORT_FAB_ARG_3D(q),
		    BL_FORT_FAB_ARG_3D(uprime));
    void narrow_diffterm_3d(const int* tlo, const int*thi, const amrex::Real* dx,
			    const BL_FORT_FAB_ARG_3D(q),
			    BL_FORT_FAB_ARG_3D(uprime),
			    const amrex::Real* mu, const amrex::Real* xi, const amrex::Real* lam, const amrex::Real* Ddiag);
#ifdef __cplusplus
}
#endif

#endif

