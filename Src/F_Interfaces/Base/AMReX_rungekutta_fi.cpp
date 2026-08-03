#include <AMReX_MultiFab.H>
#include <AMReX_RungeKutta.H>

#include <array>

using namespace amrex;

extern "C"
{
    using RKFillFunc = void (*) (void*, int, MultiFab*, Real);
    using RKRhsFunc = void (*) (void*, int, MultiFab*, MultiFab const*, Real, Real);
    using RKStore3Func = void (*) (void*, MultiFab const*, MultiFab const* const*);
    using RKStore4Func = void (*) (void*, MultiFab const*, MultiFab const* const*);
    using RKPostStageFunc = void (*) (void*, int, MultiFab*);

    void amrex_fi_rungekutta_rk2 (MultiFab* old_state, MultiFab* new_state,
                                  Real time, Real dt, void* ctx,
                                  RKRhsFunc rhs, RKFillFunc fill,
                                  RKPostStageFunc post_stage)
    {
        RungeKutta::RK2(*old_state, *new_state, time, dt,
            [=] (int stage, MultiFab& dudt, MultiFab const& state, Real t, Real dtsub)
            {
                rhs(ctx, stage, &dudt, &state, t, dtsub);
            },
            [=] (int stage, MultiFab& state, Real t)
            {
                fill(ctx, stage, &state, t);
            },
            [=] (int stage, MultiFab& state)
            {
                if (post_stage) {
                    post_stage(ctx, stage, &state);
                }
            });
    }

    void amrex_fi_rungekutta_rk3 (MultiFab* old_state, MultiFab* new_state,
                                  Real time, Real dt, void* ctx,
                                  RKRhsFunc rhs, RKFillFunc fill,
                                  RKStore3Func store, RKPostStageFunc post_stage)
    {
        RungeKutta::RK3(*old_state, *new_state, time, dt,
            [=] (int stage, MultiFab& dudt, MultiFab const& state, Real t, Real dtsub)
            {
                rhs(ctx, stage, &dudt, &state, t, dtsub);
            },
            [=] (int stage, MultiFab& state, Real t)
            {
                fill(ctx, stage, &state, t);
            },
            [=] (Array<MultiFab,3> const& rkk)
            {
                if (store) {
                    std::array<MultiFab const*,3> rkp{{rkk.data(), &rkk[1], &rkk[2]}};
                    store(ctx, old_state, rkp.data());
                }
            },
            [=] (int stage, MultiFab& state)
            {
                if (post_stage) {
                    post_stage(ctx, stage, &state);
                }
            });
    }

    void amrex_fi_rungekutta_rk4 (MultiFab* old_state, MultiFab* new_state,
                                  Real time, Real dt, void* ctx,
                                  RKRhsFunc rhs, RKFillFunc fill,
                                  RKStore4Func store, RKPostStageFunc post_stage)
    {
        RungeKutta::RK4(*old_state, *new_state, time, dt,
            [=] (int stage, MultiFab& dudt, MultiFab const& state, Real t, Real dtsub)
            {
                rhs(ctx, stage, &dudt, &state, t, dtsub);
            },
            [=] (int stage, MultiFab& state, Real t)
            {
                fill(ctx, stage, &state, t);
            },
            [=] (Array<MultiFab,4> const& rkk)
            {
                if (store) {
                    std::array<MultiFab const*,4> rkp{{rkk.data(), &rkk[1], &rkk[2], &rkk[3]}};
                    store(ctx, old_state, rkp.data());
                }
            },
            [=] (int stage, MultiFab& state)
            {
                if (post_stage) {
                    post_stage(ctx, stage, &state);
                }
            });
    }
}
