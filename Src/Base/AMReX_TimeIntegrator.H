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

namespace amrex {

enum struct IntegratorTypes {
    ForwardEuler = 0,
    ExplicitRungeKutta,
    Sundials
};

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
        // By default, do nothing in the RHS
        set_rhs([](T& /* S_rhs */, const T& /* S_data */, const amrex::Real /* time */){});
        set_imex_rhs([](T& /* S_rhs */, const T& /* S_data */, const amrex::Real /* time */){},
                     [](T& /* S_rhs */, const T& /* S_data */, const amrex::Real /* time */){});
        set_fast_rhs([](T& /* S_rhs */, const T& /* S_data */, const amrex::Real /* time */){});

        // In general, the following functions can be used to fill BCs. Which
        // function to set will depend on the method type and intended use case

        // By default, do nothing before calling the RHS
        set_pre_rhs_action([](T& /* S_data */, amrex::Real /* time */){});

        // By default, do nothing after a stage or step
        set_post_stage_action([](T& /* S_data */, const amrex::Real /* time */){});
        set_post_step_action([](T& /* S_data */, const amrex::Real /* time */){});

        // By default, do nothing after a stage or step
        set_post_fast_stage_action([](T& /* S_data */, const amrex::Real /* time */){});
        set_post_fast_step_action([](T& /* S_data */, const amrex::Real /* time */){});
    }

public:

    TimeIntegrator () {
        // initialize functions to do nothing
        set_default_functions();
    }

    TimeIntegrator (IntegratorTypes integrator_type, const T& S_data, const amrex::Real time = 0.0)
    {
        // initialize the integrator class corresponding to the desired type
        initialize_integrator(integrator_type, S_data, time);

        // initialize functions to do nothing
        set_default_functions();
    }

    TimeIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        // initialize the integrator class corresponding to the input parameter selection
        IntegratorTypes integrator_type = read_parameters();
        initialize_integrator(integrator_type, S_data, time);

        // initialize functions to do nothing
        set_default_functions();
    }

    virtual ~TimeIntegrator () {}

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

    void set_rhs (std::function<void(T&, const T&, const amrex::Real)> F)
    {
        integrator_ptr->set_rhs(F);
    }

    void set_imex_rhs (std::function<void(T&, const T&, const amrex::Real)> Fi,
                       std::function<void(T&, const T&, const amrex::Real)> Fe)
    {
        integrator_ptr->set_imex_rhs(Fi, Fe);
    }

    void set_fast_rhs (std::function<void(T&, const T&, const amrex::Real)> F)
    {
        integrator_ptr->set_fast_rhs(F);
    }

    void set_pre_rhs_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_pre_rhs_action(A);
    }

    void set_post_stage_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_stage_action(A);
    }

    void set_post_step_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_step_action(A);
    }

    void set_post_fast_stage_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_fast_stage_action(A);
    }

    void set_post_fast_step_action (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_fast_step_action(A);
    }

    void set_post_update (std::function<void (T&, amrex::Real)> A)
    {
        integrator_ptr->set_post_update(A);
    }

    amrex::Real get_time_step ()
    {
        return integrator_ptr->get_time_step();
    }

    void set_time_step (amrex::Real dt)
    {
        integrator_ptr->set_time_step(dt);
    }

    void set_adaptive_step ()
    {
        integrator_ptr->set_adaptive_step();
    }

    void set_fast_time_step (amrex::Real dt)
    {
        integrator_ptr->set_fast_time_step(dt);
    }

    void set_adaptive_fast_step ()
    {
        integrator_ptr->set_adaptive_fast_step();
    }

    void set_max_steps (int steps)
    {
        integrator_ptr->set_max_steps(steps);
    }

    void set_tolerances (amrex::Real rtol, amrex::Real atol)
    {
        integrator_ptr->set_tolerances(rtol, atol);
    }

    void set_fast_tolerances (amrex::Real rtol, amrex::Real atol)
    {
        integrator_ptr->set_fast_tolerances(rtol, atol);
    }

    void advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real dt)
    {
        integrator_ptr->advance(S_old, S_new, time, dt);
    }

    void evolve (T& S_out, const amrex::Real time_out)
    {
        integrator_ptr->evolve(S_out, time_out);
    }

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

    void time_interpolate (const T& S_new, const T& S_old, amrex::Real timestep_fraction, T& data)
    {
        integrator_ptr->time_interpolate(S_new, S_old, timestep_fraction, data);
    }

    void map_data (std::function<void(T&)> Map)
    {
        integrator_ptr->map_data(Map);
    }
};

}

#endif
