#ifndef AMREX_RUNGE_KUTTA_H_
#define AMREX_RUNGE_KUTTA_H_
#include <AMReX_Config.H>

#include <AMReX_FabArray.H>

namespace amrex::RungeKutta {

/**
 * \brief Functions for Runge-Kutta methods
 *
 * This namespace RungeKutta has functions for a number RK methods, RK2, RK3
 * and RK4.  Here, RK2 refers to the explicit trapezoid rule, RK3 refers to
 * the SSPRK3
 * (https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Third-order_Strong_Stability_Preserving_Runge-Kutta_(SSPRK3)),
 * and RK4 is the classical fourth-order method
 * (https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Classic_fourth-order_method).
 * The function templates take the old data in FabArray/MultiFab as input,
 * and evolve the system for one time step.  The result is stored in another
 * FabArray/MultiFab.  These two FabArrays must have ghost cells if they are
 * needed for evaluating the right-hand side.  The functions take three
 * callable objects for computing the right-hand side, filling ghost cells,
 * and optionally post-processing RK stage results.  For RK3 and RK4, they
 * also need a callable object for storing the data needed for filling
 * coarse/fine boundaries in AMR simulations.
 *
 * The callable object for right-hand side has the signature of `void(int
 * stage, MF& dudt, MF const& u, Real t, Real dt)`, where `stage` is the RK
 * stage number starting from 1, `dudt` is the output, `u` is the input, `t`
 * is the first-order approximate time of the stage, and `dt` is the
 * sub-time step, which can be used for reflux operations in AMR
 * simulations.
 *
 * The callable object for filling ghost cells has the signature of
 * `void(int stage, MF& u, Real t)`, where `stage` is the RK stage number
 * starting from 1, `u` is a FabArray/MultiFab whose ghost cells need to be
 * filled, and `t` is the first-order approximate time of the data at that
 * stage.  The FillPatcher class can be useful for implementing such a
 * callable.  See AmrLevel::RK for an example.
 *
 * The callable object for post-processing stage results is optional.  It's
 * no-op by default.  Its function signature is `void(int stage, MF& u)`,
 * where `stage` is the RK stage number and `u` is the result of that stage.
 *
 * For RK3 and RK4, one must also provide a callable object with the
 * signature of `void(Array<MF,order> const& rkk)`, where `order` is the RK
 * order and `rkk` contains the right-hand side at all the RK stages.  The
 * FillPatcher class can be useful for implementing such a callable.  See
 * AmrLevel::RK for an example.
 */

struct PostStageNoOp {
    template <typename MF>
    std::enable_if_t<IsFabArray<MF>::value> operator() (int, MF&) const {}
};

namespace detail {
//! Unew = Uold + dUdt * dt
template <typename MF>
void rk_update (MF& Unew, MF const& Uold, MF const& dUdt, Real dt)
{
    auto const& snew = Unew.arrays();
    auto const& sold = Uold.const_arrays();
    auto const& sdot = dUdt.const_arrays();
    amrex::ParallelFor(Unew, IntVect(0), Unew.nComp(), [=] AMREX_GPU_DEVICE
                       (int bi, int i, int j, int k, int n) noexcept
    {
            snew[bi](i,j,k,n) = sold[bi](i,j,k,n) + dt*sdot[bi](i,j,k,n);
    });
    Gpu::streamSynchronize();
}

//! Unew = Uold + (dUdt1 + dUdt2) * dt
template <typename MF>
void rk_update (MF& Unew, MF const& Uold, MF const& dUdt1, MF const& dUdt2, Real dt)
{
    auto const& snew = Unew.arrays();
    auto const& sold = Uold.const_arrays();
    auto const& sdot1 = dUdt1.const_arrays();
    auto const& sdot2 = dUdt2.const_arrays();
    amrex::ParallelFor(Unew, IntVect(0), Unew.nComp(), [=] AMREX_GPU_DEVICE
                       (int bi, int i, int j, int k, int n) noexcept
    {
            snew[bi](i,j,k,n) = sold[bi](i,j,k,n) + dt*(sdot1[bi](i,j,k,n) +
                                                        sdot2[bi](i,j,k,n));
    });
    Gpu::streamSynchronize();
}

//! Unew = (Uold+Unew)/2 + dUdt * dt/2
template <typename MF>
void rk2_update_2 (MF& Unew, MF const& Uold, MF const& dUdt, Real dt)
{
    auto const& snew = Unew.arrays();
    auto const& sold = Uold.const_arrays();
    auto const& sdot = dUdt.const_arrays();
    amrex::ParallelFor(Unew, IntVect(0), Unew.nComp(), [=] AMREX_GPU_DEVICE
                       (int bi, int i, int j, int k, int n) noexcept
    {
        snew[bi](i,j,k,n) = Real(0.5)*(snew[bi](i,j,k,n) +
                                       sold[bi](i,j,k,n) +
                                       sdot[bi](i,j,k,n) * dt);
    });
    Gpu::streamSynchronize();
}

//! Unew = Uold + (k1 + k2 + 4*k3) * dt6, where dt6 = dt/6
template <typename MF>
void rk3_update_3 (MF& Unew, MF const& Uold, Array<MF,3> const& rkk, Real dt6)
{
    auto const& snew = Unew.arrays();
    auto const& sold = Uold.const_arrays();
    auto const& k1 = rkk[0].const_arrays();
    auto const& k2 = rkk[1].const_arrays();
    auto const& k3 = rkk[2].const_arrays();
    amrex::ParallelFor(Unew, IntVect(0), Unew.nComp(), [=] AMREX_GPU_DEVICE
                       (int bi, int i, int j, int k, int n) noexcept
    {
        snew[bi](i,j,k,n) = sold[bi](i,j,k,n)
            + dt6 * (k1[bi](i,j,k,n) + k2[bi](i,j,k,n)
                     +      Real(4.) * k3[bi](i,j,k,n));
    });
    Gpu::streamSynchronize();
}

//! Unew = Uold + (k1+k4+2*(k2+k3))*dt6, where dt6 = dt/6
template <typename MF>
void rk4_update_4 (MF& Unew, MF const& Uold, Array<MF,4> const& rkk, Real dt6)
{
    auto const& snew = Unew.arrays();
    auto const& sold = Uold.const_arrays();
    auto const& k1 = rkk[0].const_arrays();
    auto const& k2 = rkk[1].const_arrays();
    auto const& k3 = rkk[2].const_arrays();
    auto const& k4 = rkk[3].const_arrays();
    amrex::ParallelFor(Unew, IntVect(0), Unew.nComp(), [=] AMREX_GPU_DEVICE
                       (int bi, int i, int j, int k, int n) noexcept
    {
        snew[bi](i,j,k,n) = sold[bi](i,j,k,n)
            + dt6 * (            k1[bi](i,j,k,n) + k4[bi](i,j,k,n)
                     + Real(2.)*(k2[bi](i,j,k,n) + k3[bi](i,j,k,n)));
    });
    Gpu::streamSynchronize();
}
}

/**
 * \brief Time stepping with RK2
 *
 * \param Uold       input FabArray/MultiFab data at time
 * \param Unew       output FabArray/MultiFab data at time+dt
 * \param time       time at the beginning of the step
 * \param dt         time step
 * \param frhs       computing the right-hand side
 * \param fillbndry  filling ghost cells
 * \param post_stage post-processing stage results
 */
template <typename MF, typename F, typename FB, typename P = PostStageNoOp>
void RK2 (MF& Uold, MF& Unew, Real time, Real dt, F&& frhs, FB&& fillbndry,
          P&& post_stage = PostStageNoOp())
{
    BL_PROFILE("RungeKutta2");

    MF dUdt(Unew.boxArray(), Unew.DistributionMap(), Unew.nComp(), 0,
            MFInfo(), Unew.Factory());

    // RK2 stage 1
    fillbndry(1, Uold, time);
    frhs(1, dUdt, Uold, time, Real(0.5)*dt);
    // Unew = Uold + dt * dUdt
    detail::rk_update(Unew, Uold, dUdt, dt);
    post_stage(1, Unew);

    // RK2 stage 2
    fillbndry(2, Unew, time+dt);
    frhs(2, dUdt, Unew, time, Real(0.5)*dt);
    // Unew = (Uold+Unew)/2 + dUdt_2 * dt/2,
    // which is Unew = Uold + dt/2 * (dUdt_1 + dUdt_2)
    detail::rk2_update_2(Unew, Uold, dUdt, dt);
    post_stage(2, Unew);
}

/**
 * \brief Time stepping with RK3
 *
 * \param Uold            input FabArray/MultiFab data at time
 * \param Unew            output FabArray/MultiFab data at time+dt
 * \param time            time at the beginning of the step
 * \param dt              time step
 * \param frhs            computing the right-hand side
 * \param fillbndry       filling ghost cells
 * \param store_crse_data storing right-hand side data for AMR
 * \param post_stage      post-processing stage results
 */
template <typename MF, typename F, typename FB, typename R,
          typename P = PostStageNoOp>
void RK3 (MF& Uold, MF& Unew, Real time, Real dt, F&& frhs, FB&& fillbndry,
          R&& store_crse_data, P&& post_stage = PostStageNoOp())
{
    BL_PROFILE("RungeKutta3");

    Array<MF,3> rkk;
    for (auto& mf : rkk) {
        mf.define(Unew.boxArray(), Unew.DistributionMap(), Unew.nComp(), 0,
                  MFInfo(), Unew.Factory());
    }

    // RK3 stage 1
    fillbndry(1, Uold, time);
    frhs(1, rkk[0], Uold, time, dt/Real(6.));
    // Unew = Uold + k1 * dt
    detail::rk_update(Unew, Uold, rkk[0], dt);
    post_stage(1, Unew);

    // RK3 stage 2
    fillbndry(2, Unew, time+dt);
    frhs(2, rkk[1], Unew, time+dt, dt/Real(6.));
    // Unew = Uold + (k1+k2) * dt/4
    detail::rk_update(Unew, Uold, rkk[0], rkk[1], Real(0.25)*dt);
    post_stage(2, Unew);

    // RK3 stage 3
    Real t_half = time + Real(0.5)*dt;
    fillbndry(3, Unew, t_half);
    frhs(3, rkk[2], Unew, t_half, dt*Real(2./3.));
    // Unew = Uold + (k1/6 + k2/6 + k3*(2/3)) * dt
    detail::rk3_update_3(Unew, Uold, rkk, Real(1./6.)*dt);
    post_stage(3, Unew);

    store_crse_data(rkk);
}

/**
 * \brief Time stepping with RK4
 *
 * \param Uold            input FabArray/MultiFab data at time
 * \param Unew            output FabArray/MultiFab data at time+dt
 * \param time            time at the beginning of the step
 * \param dt              time step
 * \param frhs            computing the right-hand side
 * \param fillbndry       filling ghost cells
 * \param store_crse_data storing right-hand side data for AMR
 * \param post_stage      post-processing stage results
 */
template <typename MF, typename F, typename FB, typename R,
          typename P = PostStageNoOp>
void RK4 (MF& Uold, MF& Unew, Real time, Real dt, F&& frhs, FB&& fillbndry,
          R&& store_crse_data, P&& post_stage = PostStageNoOp())
{
    BL_PROFILE("RungeKutta4");

    Array<MF,4> rkk;
    for (auto& mf : rkk) {
        mf.define(Unew.boxArray(), Unew.DistributionMap(), Unew.nComp(), 0,
                  MFInfo(), Unew.Factory());
    }

    // RK4 stage 1
    fillbndry(1, Uold, time);
    frhs(1, rkk[0], Uold, time, dt/Real(6.));
    // Unew = Uold + k1 * dt/2
    detail::rk_update(Unew, Uold, rkk[0], Real(0.5)*dt);
    post_stage(1, Unew);

    // RK4 stage 2
    Real t_half = time + Real(0.5)*dt;
    fillbndry(2, Unew, t_half);
    frhs(2, rkk[1], Unew, t_half, dt/Real(3.));
    // Unew = Uold + k2 * dt/2
    detail::rk_update(Unew, Uold, rkk[1], Real(0.5)*dt);
    post_stage(2, Unew);

    // RK4 stage 3
    fillbndry(3, Unew, t_half);
    frhs(3, rkk[2], Unew, t_half, dt/Real(3.));
    // Unew = Uold + k3 * dt;
    detail::rk_update(Unew, Uold, rkk[2], dt);
    post_stage(3, Unew);

    // RK4 stage 4
    fillbndry(4, Unew, time+dt);
    frhs(4, rkk[3], Unew, time+dt, dt/Real(6.));
    // Unew = Uold + (k1/6 + k2/3 + k3/3 + k4/6) * dt
    detail::rk4_update_4(Unew, Uold, rkk, Real(1./6.)*dt);
    post_stage(4, Unew);

    store_crse_data(rkk);
}

}

#endif
