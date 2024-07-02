#ifndef AMREX_RK_INTEGRATOR_H
#define AMREX_RK_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <functional>

namespace amrex {

enum struct ButcherTableauTypes {
    User = 0,
    ForwardEuler,
    Trapezoid,
    SSPRK3,
    RK4,
    NumTypes
};

template<class T>
class RKIntegrator : public IntegratorBase<T>
{
private:
    using BaseT = IntegratorBase<T>;

    // Butcher tableau identifiers
    ButcherTableauTypes tableau_type;

    // Number of RK stages
    int number_nodes;

    // A matrix from Butcher tableau
    amrex::Vector<amrex::Vector<amrex::Real> > tableau;

    // b vector from Butcher tableau
    amrex::Vector<amrex::Real> weights;

    // c vector from Butcher tableau
    amrex::Vector<amrex::Real> nodes;

    // RK embedded method b vector
    amrex::Vector<amrex::Real> extended_weights;

    // RK stage right-hand sides
    amrex::Vector<std::unique_ptr<T> > F_nodes;

    // Current (internal) state and time
    amrex::Vector<std::unique_ptr<T> > S_current;
    amrex::Real time_current;

    void initialize_preset_tableau ()
    {
        switch (tableau_type)
        {
            case ButcherTableauTypes::ForwardEuler:
                nodes = {0.0};
                tableau = {{0.0}};
                weights = {1.0};
                break;
            case ButcherTableauTypes::Trapezoid:
                nodes = {0.0,
                        1.0};
                tableau = {{0.0},
                        {1.0, 0.0}};
                weights = {0.5, 0.5};
                break;
            case ButcherTableauTypes::SSPRK3:
                nodes = {0.0,
                        1.0,
                        0.5};
                tableau = {{0.0},
                        {1.0, 0.0},
                        {0.25, 0.25, 0.0}};
                weights = {1./6., 1./6., 2./3.};
                break;
            case ButcherTableauTypes::RK4:
                nodes = {0.0,
                        0.5,
                        0.5,
                        1.0};
                tableau = {{0.0},
                        {0.5, 0.0},
                        {0.0, 0.5, 0.0},
                        {0.0, 0.0, 1.0, 0.0}};
                weights = {1./6., 1./3., 1./3., 1./6.};
                break;
            default:
                amrex::Error("Invalid RK Integrator tableau type");
                break;
        }

        number_nodes = weights.size();
    }

    void initialize_parameters ()
    {
        amrex::ParmParse pp("integration.rk");

        // Read an integrator type, if not recognized, then read weights/nodes/butcher tableau
        int _tableau_type = 0;
        pp.get("type", _tableau_type);
        tableau_type = static_cast<ButcherTableauTypes>(_tableau_type);

        // By default, define no extended weights
        extended_weights = {};

        if (tableau_type == ButcherTableauTypes::User)
        {
            // Read weights/nodes/butcher tableau
            pp.getarr("weights", weights);
            pp.queryarr("extended_weights", extended_weights);
            pp.getarr("nodes", nodes);

            amrex::Vector<amrex::Real> btable; // flattened into row major format
            pp.getarr("tableau", btable);

            // Sanity check the inputs
            if (weights.size() != nodes.size())
            {
                amrex::Error("integration.rk.weights should be the same length as integration.rk.nodes");
            } else {
                number_nodes = weights.size();
                const int nTableau = (number_nodes * (number_nodes + 1)) / 2; // includes diagonal
                if (btable.size() != nTableau)
                {
                    amrex::Error("integration.rk.tableau incorrect length - should include the Butcher Tableau diagonal.");
                }
            }

            // Fill tableau from the flattened entries
            int k = 0;
            for (int i = 0; i < number_nodes; ++i)
            {
                amrex::Vector<amrex::Real> stage_row;
                for (int j = 0; j <= i; ++j)
                {
                    stage_row.push_back(btable[k]);
                    ++k;
                }

                tableau.push_back(stage_row);
            }

            // Check that this is an explicit method
            for (const auto& astage : tableau)
            {
                if (astage.back() != 0.0)
                {
                    amrex::Error("RKIntegrator currently only supports explicit Butcher tableaus.");
                }
            }
        } else if (tableau_type > ButcherTableauTypes::User && tableau_type < ButcherTableauTypes::NumTypes)
        {
            initialize_preset_tableau();
        } else {
            amrex::Error("RKIntegrator received invalid input for integration.rk.type");
        }
    }

    void initialize_stages (const T& S_data, const amrex::Real time)
    {
        // Create data for stage RHS
        for (int i = 0; i < number_nodes; ++i)
        {
            IntegratorOps<T>::CreateLike(F_nodes, S_data);
        }

        // Create and initialize data for current state
        IntegratorOps<T>::CreateLike(S_current, S_data, true);
        IntegratorOps<T>::Copy(*S_current[0], S_data);

        // Set the initial time
        time_current = time;
    }

public:
    RKIntegrator () {}

    RKIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize(S_data, time);
    }

    void initialize (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize_parameters();
        initialize_stages(S_data, time);
    }

    virtual ~RKIntegrator () {}

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real dt) override
    {
        // Assume before step() that S_old is valid data at the current time ("time" argument)
        // And that if data is a MultiFab, both S_old and S_new contain ghost cells for evaluating a stencil based RHS
        // We need this from S_old. This is convenient for S_new to have so we can use it
        // as scratch space for stage values without creating a new scratch MultiFab with ghost cells.

        // Fill the RHS F_nodes at each stage
        for (int i = 0; i < number_nodes; ++i)
        {
            // Get current stage time, t = t_old + h * Ci
            amrex::Real stage_time = time + dt * nodes[i];

            // Fill S_new with the solution value for evaluating F at the current stage
            // Copy S_new = S_old
            IntegratorOps<T>::Copy(S_new, S_old);
            if (i > 0) {
                // Saxpy across the tableau row:
                // S_new += h * Aij * Fj
                // We should fuse these kernels ...
                for (int j = 0; j < i; ++j)
                {
                    IntegratorOps<T>::Saxpy(S_new, dt * tableau[i][j], *F_nodes[j]);
                }

                BaseT::post_stage_action(S_new, stage_time);
            }

            // Call the update hook for the stage state value
            BaseT::pre_rhs_action(S_new, stage_time);

            // Fill F[i], the RHS at the current stage
            // F[i] = RHS(y, t) at y = stage_value, t = stage_time
            BaseT::Rhs(*F_nodes[i], S_new, stage_time);
        }

        // Fill new State, starting with S_new = S_old.
        // Then Saxpy S_new += h * Wi * Fi for integration weights Wi
        // We should fuse these kernels ...
        IntegratorOps<T>::Copy(S_new, S_old);
        for (int i = 0; i < number_nodes; ++i)
        {
            IntegratorOps<T>::Saxpy(S_new, dt * weights[i], *F_nodes[i]);
        }

        BaseT::post_step_action(S_new, time + dt);

        // If we are working with an extended Butcher tableau, we can estimate the error here,
        // and then calculate an adaptive time step.

        // Save last completed step size for time_interpolate
        BaseT::previous_time_step = dt;

        // Return time step
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

            // Save last completed step size for time_interpolate
            BaseT::previous_time_step = dt;

            if (step_number == BaseT::max_steps - 1) {
                Error("Did not reach output time in max steps.");
            }
        }
    }

    void time_interpolate (const T& /* S_new */, const T& S_old, amrex::Real timestep_fraction, T& data) override
    {
        // data = S_old*(1-time_step_fraction) + S_new*(time_step_fraction)
        /*
        data.setVal(0);

        IntegratorOps<T>::Saxpy(data, 1-timestep_fraction, S_old);
        IntegratorOps<T>::Saxpy(data, timestep_fraction, S_new);
        */

        // currently we only do this for 4th order RK
        AMREX_ASSERT(number_nodes == 4);

        // fill data using MC Equation 39 at time + timestep_fraction * dt
        amrex::Real c = 0;

        // data = S_old
        IntegratorOps<T>::Copy(data, S_old);

        // data += (chi - 3/2 * chi^2 + 2/3 * chi^3) * k1
        c = timestep_fraction - 1.5 * std::pow(timestep_fraction, 2) + 2./3. * std::pow(timestep_fraction, 3);
        IntegratorOps<T>::Saxpy(data, c*BaseT::previous_time_step, *F_nodes[0]);

        // data += (chi^2 - 2/3 * chi^3) * k2
        c = std::pow(timestep_fraction, 2) - 2./3. * std::pow(timestep_fraction, 3);
        IntegratorOps<T>::Saxpy(data, c*BaseT::previous_time_step, *F_nodes[1]);

        // data += (chi^2 - 2/3 * chi^3) * k3
        c = std::pow(timestep_fraction, 2) - 2./3. * std::pow(timestep_fraction, 3);
        IntegratorOps<T>::Saxpy(data, c*BaseT::previous_time_step, *F_nodes[2]);

        // data += (-1/2 * chi^2 + 2/3 * chi^3) * k4
        c = -0.5 * std::pow(timestep_fraction, 2) + 2./3. * std::pow(timestep_fraction, 3);
        IntegratorOps<T>::Saxpy(data, c*BaseT::previous_time_step, *F_nodes[3]);

    }

    void map_data (std::function<void(T&)> Map) override
    {
        for (auto& F : F_nodes) {
            Map(*F);
        }
    }

};

}

#endif
