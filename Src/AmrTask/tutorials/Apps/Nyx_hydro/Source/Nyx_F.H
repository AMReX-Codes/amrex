#ifndef _Nyx_F_H_
#define _Nyx_F_H_

#include <AMReX_BLFort.H>
#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void fort_alloc_simd_vec();
  void fort_dealloc_simd_vec();

  void fort_integrate_comoving_a
    (amrex::Real* old_a, amrex::Real* new_a, amrex::Real* dt);

  void fort_integrate_comoving_a_to_z 
    (amrex::Real* old_a, amrex::Real* z_value, amrex::Real* dt);

  void fort_integrate_comoving_a_to_a 
    (amrex::Real* old_a, amrex::Real* a_value, amrex::Real* dt);

  void set_simd(const int *simd_width);

  //  void fort_get_omm    (amrex::Real* omm);
  //  void fort_get_omb    (amrex::Real* frac);
  //  void fort_get_hubble (amrex::Real* hubble);

  void fort_set_omm    (const amrex::Real& omm);
  void fort_set_omb    (const amrex::Real& frac);
  void fort_set_hubble (const amrex::Real& hubble);

  void fort_estdt
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const amrex::Real dx[], amrex::Real* dt, amrex::Real* comoving_a);

  void fort_estdt_comoving_a
    (amrex::Real* old_a, amrex::Real* new_dummy_a, amrex::Real* dt, amrex::Real* change_allowed,
     amrex::Real* fixed_da_interval, amrex::Real* final_a, int* dt_modified);

  void fort_get_method_params(int* HYP_GROW);

  void fort_set_method_params
    (const int& dm, const int& NumAdv, const int& Ndiag, const int& do_hydro,
     const int& ppm_type, const int& ppm_ref,
     const int& ppm_flatten_before_integrals,
     const int& use_colglaz, const int& use_flattening,
     const int& corner_coupling,
     const int& version_2, const int& use_const_species,
     const amrex::Real& gamma_in, const int& normalize_species,
     const int& heat_cool_type, const int& inhomo_reion);

  void fort_tabulate_rates();

  void filcc
    (const amrex::Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
     const int * domlo, const int * domhi,
     const amrex::Real * dx_crse, const amrex::Real * xlo,
     const int * bc);

  void fort_network_init();

  void fort_set_xhydrogen(amrex::Real& xhydrogen_in);

  void fort_get_num_spec(int* nspec);

  void fort_get_num_aux(int* naux);

  void fort_get_spec_names(int* spec_names, int* ispec, int* len);

  void fort_get_aux_names(int* aux_names, int* iaux, int* len);

  void fort_set_eos_params(const amrex::Real& h_species_in, const amrex::Real& he_species_in);

  void fort_set_small_values
    (const amrex::Real* average_dens, const amrex::Real* average_temp,
     const amrex::Real* comoving_a,
     const amrex::Real* small_dens,
     const amrex::Real* small_temp, const amrex::Real* small_pres);

  void fort_set_problem_params
    (const int& dm, const int* physbc_lo, const int* physbc_hi,
     const int& Outflow_value, const int& Symmetry_value,
     const int& coord_type);

  void fort_initdata
    (const int& level, const amrex::Real& time, const int* lo, const int* hi,
     const int& num_state, BL_FORT_FAB_ARG(state),
     const int& num_diag,  BL_FORT_FAB_ARG(diag_eos), const amrex::Real dx[],
     const amrex::Real xlo[], const amrex::Real xhi[]);

  void fort_init_zhi
    (const int* lo, const int* hi, 
     const int& num_diag, BL_FORT_FAB_ARG(diag_eos),
     const int& ratio, BL_FORT_FAB_ARG(zhi));

  void fort_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

  void fort_init_e_from_t
     (BL_FORT_FAB_ARG(state), const int* num_state,
      BL_FORT_FAB_ARG( diag), const int* num_diag,
      const int* lo, const int* hi, amrex::Real* comoving_a);

  void fort_init_e_from_rhoe
     (BL_FORT_FAB_ARG(state), const int* num_state,
      const int* lo, const int* hi, amrex::Real* comoving_a);

  void fort_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

  void fort_enforce_nonnegative_species
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[],
     const int* print_fortran_warnings);

  void reset_internal_e
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(D_new),
     BL_FORT_FAB_ARG(reset_e_src),
     const int* print_fortran_warnings,
     amrex::Real* comoving_a);

  void hypfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void denfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void generic_fill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void fort_make_hydro_sources
    (const amrex::Real* time,
     const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(ugdnvx),
     const BL_FORT_FAB_ARG(ugdnvy),
     const BL_FORT_FAB_ARG(ugdnvz),
     const BL_FORT_FAB_ARG(src),
           BL_FORT_FAB_ARG(hydro_src),
           BL_FORT_FAB_ARG(divu_cc),
     const BL_FORT_FAB_ARG(grav),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(const BL_FORT_FAB_ARG(xflux),
            const BL_FORT_FAB_ARG(yflux),
            const BL_FORT_FAB_ARG(zflux)),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const int* print_fortran_warnings);

  void fort_update_state
    ( const int lo[], const int hi[],
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(src),
      const BL_FORT_FAB_ARG(hydro_src),
      const BL_FORT_FAB_ARG(divu_cc),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new,
      const int* print_fortran_warnings);

  void fort_add_grav_source
    ( const int lo[], const int hi[],
      const BL_FORT_FAB_ARG(u_in),
            BL_FORT_FAB_ARG(u_out),
      const BL_FORT_FAB_ARG(grav),
      const amrex::Real* dt,
      const amrex::Real* a_old,
      const amrex::Real* a_new);

  void time_center_sources
    (const int lo[], const int hi[], BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const amrex::Real* dt, const int* print_fortran_warnings);

  void adjust_heat_cool
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     BL_FORT_FAB_ARG(ext_src_old), BL_FORT_FAB_ARG(ext_src_new),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const amrex::Real* dt);

  void fort_correct_gsrc
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(grav_src_old),
     const BL_FORT_FAB_ARG(grad_phi_cc),
     const BL_FORT_FAB_ARG(S_old), BL_FORT_FAB_ARG(S_new),
     const amrex::Real* a_old, const amrex::Real* a_new,
     const amrex::Real* dt);

  void fort_syncgsrc
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(gPhi),
     const BL_FORT_FAB_ARG(gdPhi),
     const BL_FORT_FAB_ARG(S),
     const BL_FORT_FAB_ARG(dS),
     BL_FORT_FAB_ARG(src),
     const amrex::Real* a_new,
     const amrex::Real& dt);

  void sum_over_level
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const amrex::Real dx[], amrex::Real* sum);

  void sum_product
    (BL_FORT_FAB_ARG(fab1), BL_FORT_FAB_ARG(fab2),
     const int lo[], const int hi[], const amrex::Real dx[], amrex::Real* sum);

  void sum_prod_prod
    (BL_FORT_FAB_ARG(fab1), BL_FORT_FAB_ARG(fab2),
     BL_FORT_FAB_ARG(fab3),
     const int lo[], const int hi[], const amrex::Real dx[], amrex::Real* sum);

  void fort_avgdown
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(fine_fab),
     const int ovlo[], const int ovhi[], const int rat[]);

#ifdef FORCING
  void fort_alloc_spect
    (const int* length);

  void fort_set_wavevector
    (const int k[], const int* indx);

  void fort_set_modes
    (const amrex::Real even[], const amrex::Real odd[],
     const int* length, const int* comp);
#endif

#ifdef AGN
  void fort_ext_src
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(old_diag),
     BL_FORT_FAB_ARG(new_diag),
     BL_FORT_FAB_ARG(ext_src),
     const amrex::Real* particle_locs_and_mass,
     const amrex::Real* particle_data,
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* z,
     const amrex::Real* dt);
#else
  void fort_ext_src
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(old_diag),
     BL_FORT_FAB_ARG(new_diag),
     BL_FORT_FAB_ARG(ext_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* z,
     const amrex::Real* dt);
#endif

  void integrate_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(diag_eos),
     const amrex::Real* a, const amrex::Real* delta_time,
     const int* min_iter, const int* max_iter);

  void integrate_state_with_source
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     BL_FORT_FAB_ARG(diag_eos),
     BL_FORT_FAB_ARG(hydro_src),
     BL_FORT_FAB_ARG(reset_src),
     BL_FORT_FAB_ARG(IR),
     const amrex::Real* a, const amrex::Real* delta_time,
     const int* min_iter, const int* max_iter);

  void integrate_state_fcvode_with_source
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     BL_FORT_FAB_ARG(diag_eos),
     BL_FORT_FAB_ARG(hydro_src),
     BL_FORT_FAB_ARG(reset_src),
     BL_FORT_FAB_ARG(IR),
     const amrex::Real* a, const amrex::Real* delta_time,
     const int* min_iter, const int* max_iter);

  void fort_compute_temp
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* comoving_a,
     const int* print_fortran_warnings);

  void fort_compute_temp_vec
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* comoving_a,
     const int* print_fortran_warnings);

  void fort_interp_to_this_z
    (const amrex::Real* z);

  void fort_setup_eos_params
    (amrex::Real* eos_nr_eps,
     amrex::Real* vode_rtol,
     amrex::Real* vode_atol_scaled);

  void fort_compute_max_temp_loc
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     const amrex::Real* max_temp, const amrex::Real* den_maxt,
     const int* imax, const int* jmax, const int* kmax);

  void fort_compute_rho_temp
    (const int lo[], const int hi[], const amrex::Real dx[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* rho_ave, amrex::Real* rho_T_sum, amrex::Real* T_sum,
     amrex::Real* Tinv_sum, amrex::Real* T_meanrho_sum, amrex::Real* rho_sum,
     amrex::Real* vol_sum, amrex::Real* vol_mn_sum);

  void fort_compute_gas_frac
    (const int lo[], const int hi[], const amrex::Real dx[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(diag_eos),
     amrex::Real* rho_ave, amrex::Real* T_cut, amrex::Real* rho_cut,
     amrex::Real* whim_mass, amrex::Real* whim_vol,
     amrex::Real* hh_mass, amrex::Real* hh_vol,
     amrex::Real* igm_mass, amrex::Real* igm_vol,
     amrex::Real* mass_sum, amrex::Real* vol_sum);

#ifdef AUX_UPDATE
  void auxupdate
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const amrex::Real * dt);
#endif

  void get_rhoe
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(rhoe),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(ye),
   const BL_FORT_FAB_ARG(state));
#ifdef __cplusplus
}
#endif

#endif
