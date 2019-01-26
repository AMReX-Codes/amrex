#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void update_agn_particles(const int* np, void* particles,
                              const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
                              const amrex::Real* prob_lo,
                              const amrex::Real* dx, const amrex::Real& dt,
                              const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);

    void nyx_compute_overlap(const int* np, void* particles, 
                             const int* ng, void*    ghosts, const amrex::Real* dx);

    void agn_merge_particles(const int* np, void* particles, 
                             const int* ng, void*    ghosts, const amrex::Real* dx);

    void agn_particle_velocity(const int* np, void* particles,
                               const amrex::Real* state_old, const int* sold_lo, const int* sold_hi,
                               const amrex::Real* state_new, const int* snew_lo, const int* snew_hi,
                               const amrex::Real* dx, const int* add_energy);

  void agn_accrete_mass(const int* np, void* particles,
                        const amrex::Real* state, const amrex::Real* density_lost,
                        const int* s_lo, const int* s_hi,
                        const amrex::Real* dt, const amrex::Real* dx);

  void agn_release_energy(const int* np, void* particles,
                          const amrex::Real* state,
                          const int* s_lo, const int* s_hi,
                          const amrex::Real* diag_eos,
                          const int* d_lo, const int* d_hi,
                          const amrex::Real* a,
                          const amrex::Real* dx);

  void init_uniform01_rng();

#ifdef __cplusplus
}
#endif
