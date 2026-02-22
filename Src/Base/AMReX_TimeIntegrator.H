#ifndef AMREX_TIME_INTEGRATOR_H
#define AMREX_TIME_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <AMReX_FEIntegrator.H>
#include <AMReX_RKIntegrator.H>

#ifdef AMREX_USE_SUNDIALS
#include <AMReX_SundialsIntegrator.H>
#endif

#include <functional>

/**
 * \file AMReX_TimeIntegrator.H
 * \brief Runtime-configurable wrapper that selects among supported integrators.
 */

namespace amrex {

/**
 * \enum IntegratorTypes
 * \brief Enumerates the concrete integrator backends available at runtime.
 */
enum struct IntegratorTypes {
    ForwardEuler = 0,
    ExplicitRungeKutta,
    Sundials
};

/**
 * \class amrex::TimeIntegrator
 * \brief Runtime-selected wrapper that forwards to a concrete integrator backend.
 *
 * The template parameter \c T represents the state container operated on by the
 * integrator.  Callers configure RHS callbacks, tolerances, and optional hooks,
 * then drive the integration via advance/evolve/integrate helpers.
 */
template<class T>
class TimeIntegrator
{
private:
    std::unique_ptr<IntegratorBase<T> > integrator_ptr;

    IntegratorTypes read_parameters ()
    {
        amrex::ParmParse pp("integration");

        int integrator_type;
        std::string integrator_str;
        pp.get("type", integrator_str);

        if (integrator_str == "ForwardEuler") {
            integrator_type = static_cast<int>(IntegratorTypes::ForwardEuler);
        } else if (integrator_str == "RungeKutta") {
            integrator_type = static_cast<int>(IntegratorTypes::ExplicitRungeKutta);
        } else if (integrator_str == "SUNDIALS") {
            integrator_type = static_cast<int>(IntegratorTypes::Sundials);
        } else {
            try {
                integrator_type = std::stoi(integrator_str, nullptr);
            } catch (const std::invalid_argument& ia) {
                Print() << "Invalid integration.type: " << ia.what() << '\n';
                Error("Failed to initialize AMReX TimeIntegrator class.");
            }

            AMREX_ALWAYS_ASSERT(integrator_type >= static_cast<int>(IntegratorTypes::ForwardEuler) &&
                                integrator_type <= static_cast<int>(IntegratorTypes::Sundials));
        }

#ifndef AMREX_USE_SUNDIALS
        if (integrator_type == static_cast<int>(IntegratorTypes::Sundials)) {
            Error("AMReX has not been compiled with SUNDIALS. Recompile with USE_SUNDIALS=TRUE.");
        }
#endif

        return static_cast<IntegratorTypes>(integrator_type);
    }

    void set_default_functions ()
    {
        // Abort if the RHS function is called but was not set
        set_rhs([](T& /* S_rhs */, T& /* S_data */, const amrex::Real /* time */) { amrex::Abort("RHS function not set!"); });
        set_imex_rhs([](T& /* S_rhs */, T& /* S_data */, const amrex::Real /* time */) { amrex::Abort("Implicit RHS function not set!"); },
                     [](T& /* S_rhs */, T& /* S_data */, const amrex::Real /* time */) { amrex::Abort("Explicit RHS function not set!"); } );
        set_fast_rhs([](T& /* S_rhs */, T& /* S_data */, const amrex::Real /* time */) { amrex::Abort("Fast RHS function not set!"); } );

        // In general, the following functions can be used to fill BCs. Which
        // function to set will depend on the method type and intended use case

        // By default, do nothing after a stage or step
        set_post_stage_action([](T& /* S_data */, const amrex::Real /* time */){});
        set_post_step_action([](T& /* S_data */, const amrex::Real /* time */){});

        // By default, do nothing after a stage or step
        set_post_fast_stage_action([](T& /* S_data */, const amrex::Real /* time */){});
        set_post_fast_step_action([](T& /* S_data */, const amrex::Real /* time */){});
    }

public:

    //! Construct an empty TimeIntegrator and install no-op callbacks.
    TimeIntegrator () {
        // initialize functions to do nothing
        set_default_functions();
    }

    /**
     * \brief Construct using a specific integrator type and initial state.
     *
     * \param integrator_type Which backend to instantiate.
     * \param S_data          Initial state copied into the backend.
     * \param time            Initial physical time associated with \p S_data.
     */
    TimeIntegrator (IntegratorTypes integrator_type, const T& S_data, const amrex::Real time = 0.0)
    {
        // initialize the integrator class corresponding to the desired type
        initialize_integrator(integrator_type, S_data, time);

        // initialize functions to do nothing
        set_default_functions();
    }

    /**
     * \brief Construct using the runtime parameter \c integration.type.
     *
     * \param S_data initial state copied into the backend.
     * \param time   initial physical time associated with \p S_data.
     */
    TimeIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        // initialize the integrator class corresponding to the input parameter selection
        IntegratorTypes integrator_type = read_parameters();
        initialize_integrator(integrator_type, S_data, time);

        // initialize functions to do nothing
        set_default_functions();
    }

    virtual ~TimeIntegrator () {}

    /**
     * \brief Instantiate the backend integrator corresponding to \p integrator_type.
     * \param integrator_type enum value returned by read_parameters().
     * \param S_data          state used to size backend storage.
     * \param time            initial physical time.
     */
    void initialize_integrator (IntegratorTypes integrator_type, const T& S_data, const amrex::Real time = 0.0)
    {
        switch (integrator_type)
        {
            case IntegratorTypes::ForwardEuler:
                integrator_ptr = std::make_unique<FEIntegrator<T> >(S_data, time);
                break;
            case IntegratorTypes::ExplicitRungeKutta:
                integrator_ptr = std::make_unique<RKIntegrator<T> >(S_data, time);
                break;
#ifdef AMREX_USE_SUNDIALS
            case IntegratorTypes::Sundials:
                integrator_ptr = std::make_unique<SundialsIntegrator<T> >(S_data, time);
                break;
#endif
            default:
                amrex::Error("integrator type did not match a valid integrator type.");
                break;
        }
    }

    /**
     * \brief Set the RHS callback used by fully explicit methods.
     *
     * The functor receives `(S_rhs, S_data, time)` where `S_rhs` is the
     * destination for the derivative, `S_data` is the current state, and
     * `time` is the physical time.
     */
    void set_rhs (std::function<void(T&, T&, const amrex::Real)> F)
    {
        integrator_ptr->set_rhs(F);
    }

    /**
     * \brief Set the implicit and explicit RHS callbacks used by IMEX methods.
     *
     * Each functor receives `(S_rhs, S_data, time)` and fills `S_rhs` with the
     * corresponding implicit or explicit contribution.
     */
    void set_imex_rhs (std::function<void(T&, T&, const amrex::Real)> Fi,
                       std::function<void(T&, T&, const amrex::Real)> Fe)
    {
        integrator_ptr->set_imex_rhs(Fi, Fe);
    }

    /**
     * \brief Set the RHS callback for the fast subsystem in multirate methods.
     *
     * The functor receives `(S_rhs, S_data, time)` just like `set_rhs`.
     */
    void set_fast_rhs (std::function<void(T&, T&, const amrex::Real)> F)
    {
        integrator_ptr->set_fast_rhs(F);
    }

    /**
     * \brief Install a functor executed after each stage (e.g., to fill ghost cells).
     *
     * The functor receives `(S_data, time_stage)` where `S_data` contains the
     * stage result.
     */
    void set_post_stage_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_stage_action(A);
    }

    /**
     * \brief Install a functor executed after each completed timestep.
     *
     * The functor receives `(S_data, time_new)` once the step is accepted.
     */
    void set_post_step_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_step_action(A);
    }

    /**
     * \brief Install a functor executed after each fast stage in multirate schemes.
     *
     * Invoked with `(S_data, time_stage)` where `S_data` contains the fast-stage state.
     */
    void set_post_fast_stage_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_fast_stage_action(A);
    }

    /**
     * \brief Install a functor executed after each fast subsystem step.
     *
     * Invoked with `(S_data, time_step_end)` once the fast step completes.
     */
    void set_post_fast_step_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_fast_step_action(A);
    }

    //! Return the current timestep selected by the backend.
    amrex::Real get_time_step ()
    {
        return integrator_ptr->get_time_step();
    }

    //! Set a fixed timestep to be used by the backend.
    void set_time_step (amrex::Real dt)
    {
        integrator_ptr->set_time_step(dt);
    }

    //! Switch the backend into adaptive-step mode.
    void set_adaptive_step ()
    {
        integrator_ptr->set_adaptive_step();
    }

    //! Set the fast-timescale timestep for multirate methods.
    void set_fast_time_step (amrex::Real dt)
    {
        integrator_ptr->set_fast_time_step(dt);
    }

    //! Switch the fast subsystem into adaptive-step mode.
    void set_adaptive_fast_step ()
    {
        integrator_ptr->set_adaptive_fast_step();
    }

    //! Limit the number of steps taken by adaptive integrators.
    void set_max_steps (int steps)
    {
        integrator_ptr->set_max_steps(steps);
    }

    //! Set relative and absolute tolerances for error-controlled timestepping.
    void set_tolerances (amrex::Real rtol, amrex::Real atol)
    {
        integrator_ptr->set_tolerances(rtol, atol);
    }

    //! Set tolerances for the fast subsystem in multirate methods.
    void set_fast_tolerances (amrex::Real rtol, amrex::Real atol)
    {
        integrator_ptr->set_fast_tolerances(rtol, atol);
    }

    /**
     * \brief Advance one step from \p time to `time+dt`.
     *
     * \param S_old input state at the start of the step.
     * \param S_new output state after the advance.
     * \param time  current physical time.
     * \param dt    timestep size.
     */
    void advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real dt)
    {
        integrator_ptr->advance(S_old, S_new, time, dt);
    }

    /**
     * \brief Repeatedly advance until the backend reaches \p time_out.
     *
     * \param S_out state updated in-place to the requested time.
     * \param time_out target physical time.
     */
    void evolve (T& S_out, const amrex::Real time_out)
    {
        integrator_ptr->evolve(S_out, time_out);
    }

    /**
     * \brief Convenience driver that advances repeatedly up to \p end_time.
     *
     * The method swaps \p S_old and \p S_new between steps so repeated calls can
     * reuse the same buffers.
     */
    void integrate (T& S_old, T& S_new, amrex::Real start_time, const amrex::Real start_timestep,
                    const amrex::Real end_time, const int start_step, const int max_steps)
    {
        amrex::Real m_time = start_time;
        amrex::Real m_timestep = start_timestep;
        bool stop = false;
        for (int m_step_number = start_step; m_step_number < max_steps && !stop; ++m_step_number)
        {
            if (end_time - m_time < m_timestep) {
                m_timestep = end_time - m_time;
                stop = true;
            }

            if (m_step_number > 0) {
                std::swap(S_old, S_new);
            }

            // Call the time integrator advance
            integrator_ptr->advance(S_old, S_new, m_time, m_timestep);

            // Update our time variable
            m_time += m_timestep;
        }
    }

    /**
     * \brief Linearly interpolate between \p S_old and \p S_new.
     *
     * \param S_new             State at the end of the step.
     * \param S_old             State at the beginning of the step.
     * \param timestep_fraction Fraction in [0,1] measuring the interpolation point.
     * \param data              Destination object that receives the interpolated state.
     */
    void time_interpolate (const T& S_new, const T& S_old, amrex::Real timestep_fraction, T& data)
    {
        integrator_ptr->time_interpolate(S_new, S_old, timestep_fraction, data);
    }

    /**
     * \brief Apply a mapping functor to all state held internally by the backend.
     * \param Map Callable invoked as `void(T&)` on each owned container.
     */
    void map_data (std::function<void(T&)> Map)
    {
        integrator_ptr->map_data(Map);
    }
};

}

#endif
