#ifndef AMREX_INTEGRATOR_BASE_H
#define AMREX_INTEGRATOR_BASE_H
#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#if defined(AMREX_PARTICLES)
#include <AMReX_Particles.H>
#endif

#include <functional>
#include <type_traits>

namespace amrex {

template<class T, typename Tv = void> struct IntegratorOps;

#if defined(AMREX_PARTICLES)
template<class T>
struct IntegratorOps<T, std::enable_if_t<std::is_base_of_v<amrex::ParticleContainerBase, T> > >
{

    static void CreateLike (amrex::Vector<std::unique_ptr<T> >& V, const T& Other)
    {
        // Emplace a new T in V with the same size as Other and get a reference
        V.emplace_back(std::make_unique<T>(Other.Geom(0), Other.ParticleDistributionMap(0), Other.ParticleBoxArray(0)));
        T& pc = *V[V.size()-1];

        // We want the particles to have all the same position, cpu, etc.
        // as in Other, so do a copy from Other to our new particle container.
        Copy(pc, Other);
    }

    static void Copy (T& Y, const T& Other)
    {
        // Copy the contents of Other into Y
        const bool local = true;
        Y.copyParticles(Other, local);
    }

    static void Saxpy (T& Y, const amrex::Real a, T& X)
    {
        // Calculate Y += a * X using a particle-level saxpy function supplied by the particle container T
        using TParIter = amrex::ParIter<T::NStructReal, T::NStructInt, T::NArrayReal, T::NArrayInt>;
        using ParticleType = amrex::Particle<T::NStructReal, T::NStructInt>;

        int lev = 0;
        TParIter pty(Y, lev);
        TParIter ptx(X, lev);

        auto checkValid = [&]() -> bool {
            bool pty_v = pty.isValid();
            bool ptx_v = ptx.isValid();
            AMREX_ASSERT(pty_v == ptx_v);
            return pty_v && ptx_v;
        };

        auto ptIncrement = [&](){ ++pty; ++ptx; };

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (; checkValid(); ptIncrement())
        {
            const int npy  = pty.numParticles();
            const int npx  = ptx.numParticles();
            AMREX_ALWAYS_ASSERT(npy == npx);

            ParticleType* psy = &(pty.GetArrayOfStructs()[0]);
            ParticleType* psx = &(ptx.GetArrayOfStructs()[0]);

            auto particle_apply_rhs = T::particle_apply_rhs;

            amrex::ParallelFor ( npy, [=] AMREX_GPU_DEVICE (int i) {
                ParticleType& py = psy[i];
                const ParticleType& px = psx[i];
                particle_apply_rhs(py, a, px);
            });
        }
    }

};
#endif

template<class T>
struct IntegratorOps<T, std::enable_if_t<std::is_same_v<amrex::Vector<amrex::MultiFab>, T> > >
{

    static void CreateLike (amrex::Vector<std::unique_ptr<T> >& V, const T& Other, bool Grow = false)
    {
        // Emplace a new T in V with the same size as Other
        V.emplace_back(std::make_unique<T>());
        for (auto const& other_mf : Other) {
            IntVect nGrow = Grow ? other_mf.nGrowVect() : IntVect(0);
            V.back()->push_back(amrex::MultiFab(other_mf.boxArray(), other_mf.DistributionMap(), other_mf.nComp(), nGrow));
        }
    }

    static void Copy (T& Y, const T& Other, const Vector<int>& scomp={}, const Vector<int>& ncomp={}, bool Grow = true)
    {
        // Copy the contents of Other into Y
        const int size = Y.size();
        bool specify_components = !scomp.empty() && ncomp.size() == scomp.size();
        for (int i = 0; i < size; ++i) {
            IntVect nGrow = Grow ? Other[i].nGrowVect() : IntVect(0);
            const int iscomp = specify_components ? scomp[i] : 0;
            const int incomp = specify_components ? ncomp[i] : Other[i].nComp();
            if (incomp > 0) {
                amrex::MultiFab::Copy(Y[i], Other[i], iscomp, iscomp, incomp, nGrow);
            }
        }
    }

    static void Saxpy (T& Y, const amrex::Real a, const T& X, const Vector<int>& scomp={}, const Vector<int>& ncomp={}, bool Grow = false)
    {
        // Calculate Y += a * X
        const int size = Y.size();
        bool specify_components = !scomp.empty() && ncomp.size() == scomp.size();
        for (int i = 0; i < size; ++i) {
            IntVect nGrow = Grow ? X[i].nGrowVect() : IntVect(0);
            const int iscomp = specify_components ? scomp[i] : 0;
            const int incomp = specify_components ? ncomp[i] : X[i].nComp();
            if (incomp > 0) {
                amrex::MultiFab::Saxpy(Y[i], a, X[i], iscomp, iscomp, incomp, nGrow);
            }
        }
    }

};

template<class T>
struct IntegratorOps<T, std::enable_if_t<std::is_same_v<amrex::MultiFab, T> > >
{

    static void CreateLike (amrex::Vector<std::unique_ptr<T> >& V, const T& Other, bool Grow = false)
    {
        // Emplace a new T in V with the same size as Other
        IntVect nGrow = Grow ? Other.nGrowVect() : IntVect(0);
        V.emplace_back(std::make_unique<T>(Other.boxArray(), Other.DistributionMap(), Other.nComp(), nGrow));
    }

    static void Copy (T& Y, const T& Other, const int scomp=0, const int ncomp=-1, bool Grow = true)
    {
        // Copy the contents of Other into Y
        IntVect nGrow = Grow ? Other.nGrowVect() : IntVect(0);
        const int mf_ncomp = ncomp > 0 ? ncomp : Other.nComp();
        amrex::MultiFab::Copy(Y, Other, scomp, scomp, mf_ncomp, nGrow);
    }

    static void Saxpy (T& Y, const amrex::Real a, const T& X, const int scomp=0, const int ncomp=-1, bool Grow = false)
    {
        // Calculate Y += a * X
        IntVect nGrow = Grow ? X.nGrowVect() : IntVect(0);
        const int mf_ncomp = ncomp > 0 ? ncomp : X.nComp();
        amrex::MultiFab::Saxpy(Y, a, X, scomp, scomp, mf_ncomp, nGrow);
    }

};

template<class T>
class IntegratorBase
{
protected:
    /**
     * \brief Rhs is the right-hand-side function the integrator will use.
     */
    std::function<void(T& rhs, const T& state, const amrex::Real time)> Rhs;

    /**
     * \brief RhsIm is the implicit right-hand-side function an ImEx integrator
     * will use.
     */
    std::function<void(T& rhs, const T& state, const amrex::Real time)> RhsIm;

    /**
     * \brief RhsEx is the explicit right-hand-side function an ImEx integrator
     * will use.
     */
    std::function<void(T& rhs, const T& state, const amrex::Real time)> RhsEx;

    /**
     * \brief RhsFast is the fast timescale right-hand-side function a multirate
     * integrator will use.
     */
    std::function<void(T& rhs, const T& state, const amrex::Real time)> RhsFast;

    /**
     * \brief The pre_rhs_action function is called by the integrator on state
     * data before using it to evaluate a right-hand side.
     */
    std::function<void (T&, amrex::Real)> pre_rhs_action;

    /**
     * \brief The post_stage_action function is called by the integrator on
     * the computed stage just after it is computed
     */
    std::function<void (T&, amrex::Real)> post_stage_action;

    /**
     * \brief The post_step_action function is called by the integrator on
     * the computed state just after it is computed
     */
    std::function<void (T&, amrex::Real)> post_step_action;

    /**
     * \brief The post_stage_action function is called by the integrator on
     * the computed stage just after it is computed
     */
    std::function<void (T&, amrex::Real)> post_fast_stage_action;

    /**
     * \brief The post_step_action function is called by the integrator on
     * the computed state just after it is computed
     */
    std::function<void (T&, amrex::Real)> post_fast_step_action;

    /**
     * \brief Flag to enable/disable adaptive time stepping in single rate
     * methods or at the slow time scale in multirate methods (bool)
     */
    bool use_adaptive_time_step = false;

    /**
     * \brief Current integrator time step size (Real)
     */
    amrex::Real time_step;

    /**
     * \brief Step size of the last completed step (Real)
     */
    amrex::Real previous_time_step;

    /**
     * \brief Flag to enable/disable adaptive time stepping at the fast time
     * scale in multirate methods (bool)
     */
    bool use_adaptive_fast_time_step = false;

    /**
     * \brief Current integrator fast time scale time step size with multirate
     * methods (Real)
     */
    amrex::Real fast_time_step;

    /**
     * \brief Number of integrator time steps (Long)
     */
    amrex::Long num_steps = 0;

    /**
     * \brief Max number of internal steps before an error is returned (Long)
     */
    int max_steps = 500;

    /**
     * \brief Relative tolerance for adaptive time stepping (Real)
     */
    amrex::Real rel_tol = 1.0e-4;

    /**
     * \brief Absolute tolerance for adaptive time stepping (Real)
     */
    amrex::Real abs_tol = 1.0e-9;

    /**
     * \brief Relative tolerance for adaptive time stepping at the fast time
     * scale (Real)
     */
    amrex::Real fast_rel_tol = 1.0e-4;

    /**
     * \brief Absolute tolerance for adaptive time stepping at the fast time
     * scale (Real)
     */
    amrex::Real fast_abs_tol = 1.0e-9;


public:
    IntegratorBase () = default;

    IntegratorBase (const T& /* S_data */) {}

    virtual ~IntegratorBase () = default;

    void set_rhs (std::function<void(T&, const T&, const amrex::Real)> F)
    {
        Rhs = F;
    }

    void set_imex_rhs (std::function<void(T&, const T&, const amrex::Real)> Fi,
                       std::function<void(T&, const T&, const amrex::Real)> Fe)
    {
        RhsIm = Fi;
        RhsEx = Fe;
    }

    void set_fast_rhs (std::function<void(T&, const T&, const amrex::Real)> F)
    {
        RhsFast = F;
    }

    void set_pre_rhs_action (std::function<void (T&, amrex::Real)> A)
    {
        pre_rhs_action = A;
    }

    void set_post_stage_action (std::function<void (T&, amrex::Real)> A)
    {
        post_stage_action = A;
    }

    void set_post_step_action (std::function<void (T&, amrex::Real)> A)
    {
        post_step_action = A;
    }

    void set_post_fast_stage_action (std::function<void (T&, amrex::Real)> A)
    {
        post_fast_stage_action = A;
    }

    void set_post_fast_step_action (std::function<void (T&, amrex::Real)> A)
    {
        post_fast_step_action = A;
    }

    void set_post_update (std::function<void (T&, amrex::Real)> A)
    {
        set_post_stage_action(A);
        set_post_step_action(A);
    }

    amrex::Real get_time_step ()
    {
        return time_step;
    }

    void set_time_step (amrex::Real dt)
    {
        time_step = dt;
        use_adaptive_time_step = false;
    }

    void set_adaptive_step ()
    {
        use_adaptive_time_step = true;
    }

    void set_fast_time_step (amrex::Real dt)
    {
        fast_time_step = dt;
        use_adaptive_fast_time_step = false;
    }

    void set_adaptive_fast_step ()
    {
        use_adaptive_fast_time_step = true;
    }

    void set_max_steps (int steps)
    {
        max_steps = steps;
    }

    void set_tolerances (amrex::Real rtol, amrex::Real atol)
    {
        rel_tol = rtol;
        abs_tol = atol;
    }

    void set_fast_tolerances (amrex::Real rtol, amrex::Real atol)
    {
        fast_rel_tol = rtol;
        fast_abs_tol = atol;
    }

    /**
     * \brief Take a single time step from (time, S_old) to (time + dt, S_new)
     * with the given step size.
     */
    virtual amrex::Real advance (T& S_old, T& S_new, amrex::Real time,
                                 amrex::Real dt) = 0;

    /**
     * \brief Evolve the current (internal) integrator state to time_out
     */
    virtual void evolve (T& S_out, const amrex::Real time_out) = 0;

    virtual void time_interpolate (const T& S_new, const T& S_old,
                                   amrex::Real timestep_fraction, T& data) = 0;

    virtual void map_data (std::function<void(T&)> Map) = 0;
};

}

#endif
