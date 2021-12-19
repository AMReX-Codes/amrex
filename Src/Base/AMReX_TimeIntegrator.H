#ifndef AMREX_TIME_INTEGRATOR_H
#define AMREX_TIME_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <AMReX_FEIntegrator.H>
#include <AMReX_RKIntegrator.H>
#include <functional>

namespace amrex {

enum struct IntegratorTypes {
    ForwardEuler = 0,
    ExplicitRungeKutta
};

template<class T>
class TimeIntegrator
{
private:
    std::unique_ptr<IntegratorBase<T> > integrator_ptr;
    std::function<void ()> post_timestep;

    IntegratorTypes read_parameters ()
    {
        amrex::ParmParse pp("integration");

        int integrator_type = static_cast<int>(IntegratorTypes::ForwardEuler);
        pp.get("type", integrator_type);

        return static_cast<IntegratorTypes>(integrator_type);
    }

    void set_default_functions ()
    {
        // By default, do nothing post-timestep
        set_post_timestep([](){});

        // By default, do nothing after updating the state
        // In general, this is where BCs should be filled
        set_post_update([](T& /* S_data */, amrex::Real /* S_time */){});

        // By default, do nothing
        set_rhs([](T& /* S_rhs */, const T& /* S_data */, const amrex::Real /* time */){});
    }

public:

    TimeIntegrator () {
        // initialize functions to do nothing
        set_default_functions();
    }

    TimeIntegrator (IntegratorTypes integrator_type, const T& S_data)
    {
        // initialize the integrator class corresponding to the desired type
        initialize_integrator(integrator_type, S_data);

        // initialize functions to do nothing
        set_default_functions();
    }

    TimeIntegrator (const T& S_data)
    {
        // initialize the integrator class corresponding to the input parameter selection
        IntegratorTypes integrator_type = read_parameters();
        initialize_integrator(integrator_type, S_data);

        // initialize functions to do nothing
        set_default_functions();
    }

    virtual ~TimeIntegrator () {}

    void initialize_integrator (IntegratorTypes integrator_type, const T& S_data)
    {
        switch (integrator_type)
        {
            case IntegratorTypes::ForwardEuler:
                integrator_ptr = std::make_unique<FEIntegrator<T> >(S_data);
                break;
            case IntegratorTypes::ExplicitRungeKutta:
                integrator_ptr = std::make_unique<RKIntegrator<T> >(S_data);
                break;
            default:
                amrex::Error("integrator type did not match a valid integrator type.");
                break;
        }
    }

    void set_post_timestep (std::function<void ()> F)
    {
        post_timestep = F;
    }

    void set_post_update (std::function<void (T&, amrex::Real)> F)
    {
        integrator_ptr->set_post_update(F);
    }

    void set_rhs (std::function<void(T&, const T&, const amrex::Real)> F)
    {
        integrator_ptr->set_rhs(F);
    }

    std::function<void ()> get_post_timestep ()
    {
        return post_timestep;
    }

    std::function<void (T&, amrex::Real)> get_post_update ()
    {
        return integrator_ptr->get_post_update();
    }

    std::function<void(T&, const T&, const amrex::Real)> get_rhs ()
    {
        return integrator_ptr->get_rhs();
    }

    void advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real timestep)
    {
        integrator_ptr->advance(S_old, S_new, time, timestep);
    }

    void integrate (T& S_old, T& S_new, amrex::Real start_time, const amrex::Real start_timestep,
                    const amrex::Real end_time, const int nsteps)
    {
        amrex::Real time = start_time;
        amrex::Real timestep = start_timestep;
        bool stop_advance = false;
        for (int step_number = 0; step_number < nsteps && !stop_advance; ++step_number)
        {
            if (end_time - time < timestep) {
                timestep = end_time - time;
                stop_advance = true;
            }

            if (step_number > 0) {
                std::swap(S_old, S_new);
            }

            // Call the time integrator advance
            integrator_ptr->advance(S_old, S_new, time, timestep);

            // Update our time variable
            time += timestep;

            // Call the post-timestep hook
            post_timestep();
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
