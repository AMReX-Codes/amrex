#ifndef _SMC_F_H_
#define _SMC_F_H_

extern "C"
{
    void derivative_stencil_init();
    void chemistry_init();
    void variables_init();

    int get_num_ghost_cells();
    int get_num_species();
    int get_num_cons();
    int get_num_prim();

    void init_data_3d(const int* lo, const int* hi, const int* clo, const int* chi, const int* ng, 
		      const Real dx[], Real* u, const Real plo[], const Real phi[]);
    void reset_rho_3d(const int* lo, const int* hi, const int* ng, Real* u);
}

#endif

