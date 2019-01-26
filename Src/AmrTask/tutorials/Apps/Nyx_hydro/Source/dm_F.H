#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void update_dm_particles(const int* np, void* particles,
                             const amrex::Real* accel, const int* accel_lo, const int* accel_hi,
                             const amrex::Real* prob_lo, 
                             const amrex::Real* dx, const amrex::Real& dt, 
                             const amrex::Real& a_prev, const amrex::Real& a_cur, const int* do_move);
#ifdef __cplusplus
}
#endif
