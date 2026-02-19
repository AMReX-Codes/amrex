#ifndef AMREX_FE_INTEGRATOR_H
#define AMREX_FE_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_IntegratorBase.H>
#include <functional>

namespace amrex {

template<class T>
class FEIntegrator : public IntegratorBase<T>
{
private:
    using BaseT = IntegratorBase<T>;

    amrex::Vector<std::unique_ptr<T> > F_nodes;

    // Current (internal) state and time
    amrex::Vector<std::unique_ptr<T> > S_current;
    amrex::Real time_current;

    void initialize_stages (const T& S_data, const amrex::Real time)
    {
        // Create data for stage RHS
        IntegratorOps<T>::CreateLike(F_nodes, S_data);

        // Create and initialize data for current state
        IntegratorOps<T>::CreateLike(S_current, S_data, true);
        IntegratorOps<T>::Copy(*S_current[0], S_data);

        // Set the initial time
        time_current = time;
    }

public:
    FEIntegrator () {}

    FEIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize(S_data, time);
    }

    virtual ~FEIntegrator () {}

    void initialize (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize_stages(S_data, time);
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real dt) override
    {
        // Assume before step() that S_old is valid data at the current time ("time" argument)
        // So we initialize S_new by copying the old state.
        IntegratorOps<T>::Copy(S_new, S_old);

        // F = RHS(S, t)
        T& F = *F_nodes[0];
        BaseT::Rhs(F, S_new, time);

        // S_new += timestep * dS/dt
        IntegratorOps<T>::Saxpy(S_new, dt, F);

        // Call the post step hook
        BaseT::post_step_action(S_new, time + dt);

        // Return timestep
        return dt;
    }

    void evolve (T& S_out, const amrex::Real time_out) override
    {
        amrex::Real dt = BaseT::time_step;
        bool stop = false;

        for (int step_number = 0; step_number < BaseT::max_steps && !stop; ++step_number)
        {
            // Adjust step size to reach output time
            if (time_out - time_current < dt) {
                dt = time_out - time_current;
                stop = true;
            }

            // Call the time integrator step
            advance(*S_current[0], S_out, time_current, dt);

            // Update current state S_current = S_out
            IntegratorOps<T>::Copy(*S_current[0], S_out);

            // Update time
            time_current += dt;

            if (step_number == BaseT::max_steps - 1) {
                Error("Did not reach output time in max steps.");
            }
        }
    }

    virtual void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override
    {
        amrex::Error("Time interpolation not yet supported by the forward euler integrator.");
    }

    virtual void map_data (std::function<void(T&)> Map) override
    {
        for (auto& F : F_nodes) {
            Map(*F);
        }
    }

};

}

#endif
