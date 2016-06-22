#ifndef _SMC_F_H_
#define _SMC_F_H_

#include <BLFort.H>

extern "C"
{
    void derivative_stencil_init();
    void chemistry_init();
    void variables_init();

    int get_num_ghost_cells();
    int get_num_species();
    int get_num_cons();
    int get_num_prim();

    void chemterm_3d(const int* lo, const int* hi, 
		     const BL_FORT_FAB_ARG_3D(q), BL_FORT_FAB_ARG_3D(uprime)); 
    void comp_courno_3d(const int* lo, const int* hi, const Real* dx, 
			const BL_FORT_FAB_ARG_3D(q), Real& courno);
    void ctoprim_3d(const int* tlo, const int* thi, const int* lo, const int* hi,
		    Real* u, Real* q, const int& ngu, const int& ngq);
    void get_trans_prop_3d(const int* tlo, const int*thi, const int* lo, const int* hi,
			   const Real* q, Real* mu, Real* xi, Real* lam, Real* Ddiag, int& ng);
    void init_data_3d(const int* tlo, const int* thi, const int* lo, const int* hi, 
		      Real* u, const int& ng, const Real dx[], const Real plo[], const Real phi[]);
    void make_plot_3d(const int* lo, const int* hi,
		      BL_FORT_FAB_ARG_3D(pv), const BL_FORT_FAB_ARG_3D(u));
    void reset_rho_3d(const int* tlo, const int* thi, const int* lo, const int* hi, 
		      Real* u, const int& ng);
    void hypterm_3d(const int* tlo, const int*thi, const Real* dx,
		    const BL_FORT_FAB_ARG_3D(u),
		    const BL_FORT_FAB_ARG_3D(q),
		    BL_FORT_FAB_ARG_3D(uprime));
    void narrow_diffterm_3d(const int* tlo, const int*thi, const Real* dx,
			    const BL_FORT_FAB_ARG_3D(q),
			    BL_FORT_FAB_ARG_3D(uprime),
			    const Real* mu, const Real* xi, const Real* lam, const Real* Ddiag);
}

#endif

