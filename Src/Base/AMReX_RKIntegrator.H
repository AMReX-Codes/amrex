#ifndef AMREX_RK_INTEGRATOR_H
#define AMREX_RK_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <functional>
#include <limits>

namespace amrex {

enum struct ButcherTableauTypes {
    User = 0,
    ForwardEuler,
    Trapezoid,
    SSPRK3,
    RK4,
    DormandPrince,
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

    // Order of the embedded error estimate (e.g. 4 for Dormand-Prince)
    int error_order;

    // Safety factor for adaptive step-size control
    amrex::Real safety_factor;

    // Min/max allowed step-size growth/shrink factors per step
    amrex::Real min_factor;
    amrex::Real max_factor;

    // RK stage right-hand sides
    amrex::Vector<std::unique_ptr<T> > F_nodes;

    // Scratch space for low-order solution in embedded error estimate
    amrex::Vector<std::unique_ptr<T> > S_error_vec;

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
            case ButcherTableauTypes::DormandPrince:
                // Dormand-Prince RK4(5) embedded pair
                nodes = {0.0,
                         1.0/5.0,
                         3.0/10.0,
                         4.0/5.0,
                         8.0/9.0,
                         1.0,
                         1.0};
                tableau = {{0.0},
                           {1.0/5.0, 0.0},
                           {3.0/40.0, 9.0/40.0, 0.0},
                           {44.0/45.0, -56.0/15.0, 32.0/9.0, 0.0},
                           {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0},
                           {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 0.0},
                           {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0}};
                // 5th-order weights
                weights = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0};
                // 4th-order weights for embedded error estimate
                extended_weights = {5179.0/57600.0, 0.0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0};
                error_order = 4;
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

        // Safety factor for adaptive step-size control (default 0.9)
        safety_factor = 0.9;
        pp.query("safety_factor", safety_factor);

        // Min/max allowed step-size growth/shrink factors per step
        min_factor = 0.1;
        max_factor = 5.0;
        pp.query("min_factor", min_factor);
        pp.query("max_factor", max_factor);

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

        // Allocate error scratch space for embedded methods
        if (!extended_weights.empty()) {
            IntegratorOps<T>::CreateLike(S_error_vec, S_data);
        }

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

        BaseT::previous_time_step = dt;

        // If we are working with an extended Butcher tableau, we can estimate the error here,
        // and then calculate an adaptive time step.
        amrex::Real dt_next = dt;
        if (!extended_weights.empty() && BaseT::use_adaptive_time_step)
        {
            // Compute low-order solution: S_error = S_old + dt * sum(extended_weights[i] * F[i])
            IntegratorOps<T>::Copy(*S_error_vec[0], S_old);
            for (int i = 0; i < number_nodes; ++i)
                IntegratorOps<T>::Saxpy(*S_error_vec[0], dt * extended_weights[i], *F_nodes[i]);

            // S_error = S_low - S_high (pointwise error estimate)
            IntegratorOps<T>::Saxpy(*S_error_vec[0], -1.0, S_new);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<bool>(BaseT::custom_error_norm),
                "Adaptive time stepping requires an error norm to be set via set_error_norm.");
            const amrex::Real scaled_err = BaseT::custom_error_norm(*S_error_vec[0], S_old, S_new,
                                                                    BaseT::abs_tol, BaseT::rel_tol);
            BaseT::last_scaled_error = scaled_err;

            // Step-size controller: exponent is 1/(error_order+1) per standard embedded RK theory
            amrex::Real factor = safety_factor * std::pow(1.0 / scaled_err, 1.0 / (error_order + 1));
            factor = amrex::min(amrex::max(factor, min_factor), max_factor);
            dt_next = dt * factor;

            if (scaled_err > 1.0)
            {
                // Step rejected: retry with a smaller dt without calling post_step_action
                return advance(S_old, S_new, time, dt_next);
            }

            // Step accepted. Set time_step to max before calling post_step_action so that
            // any CFL-based set_time_step call inside it can be detected and used to cap
            // the adaptive suggestion.
            BaseT::time_step = std::numeric_limits<amrex::Real>::max();
        }

        BaseT::post_step_action(S_new, time + dt);

        // Cap dt_next with any CFL limit set by post_step_action, then keep time_step in sync.
        dt_next = amrex::min(dt_next, BaseT::time_step);
        BaseT::time_step = dt_next;

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
            // protect against roundoff
            if ((time_out-time_current) < dt || almostEqual(time_out-time_current,dt,1000)) {
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
