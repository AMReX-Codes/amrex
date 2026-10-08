#ifndef AMREX_SUNDIALS_INTEGRATOR_H
#define AMREX_SUNDIALS_INTEGRATOR_H

#include <functional>
#include <utility>

#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <AMReX_NVector_MultiFab.H>
#include <AMReX_Sundials.H>

#include <nvector/nvector_manyvector.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>

// SUNDIALS 7.1 deprecated the stepper-specific functions in favor of
// ARKODE-wide ones.
#if defined(SUNDIALS_VERSION_MAJOR) && \
    ((SUNDIALS_VERSION_MAJOR > 7) || (SUNDIALS_VERSION_MAJOR == 7 && SUNDIALS_VERSION_MINOR >= 1))
#  define AMREX_SUNDIALS_ARKODE_API 1
#  define AMREX_ARKSTEP(f) ARKode##f
#  define AMREX_MRISTEP(f) ARKode##f
#else
#  define AMREX_SUNDIALS_ARKODE_API 0
#  define AMREX_ARKSTEP(f) ARKStep##f
#  define AMREX_MRISTEP(f) MRIStep##f
#endif

/**
 * \file AMReX_SundialsIntegrator.H
 *
 * Templates that adapt AMReX IntegratorBase to the SUNDIALS solvers (ARKStep /
 * MRIStep) plus the callback glue structures they require.
 */

namespace amrex {

/**
 * \brief User-supplied callbacks consumed by the AMReX/SUNDIALS bridge.
 *
 * Each functor must implement the signatures noted below and follow the SUNDIALS
 * convention of returning 0 on success and non-zero on failure.
 */
struct SundialsUserData {
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> f;   //!< ERK/DIRK RHS or MRI slow RHS.
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> fi;  //!< Implicit RHS for ImEx schemes.
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> fe;  //!< Explicit RHS for ImEx schemes.
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> ff;  //!< MRI fast-scale RHS.
    std::function<int(amrex::Real, N_Vector, void*)> post_stage;    //!< Hook invoked after each stage.
    std::function<int(amrex::Real, N_Vector, void*)> post_step;     //!< Hook invoked after each time step.
    std::function<int(amrex::Real, N_Vector, void*)> post_fast_stage; //!< Hook for MRI fast stages.
    std::function<int(amrex::Real, N_Vector, void*)> post_fast_step;  //!< Hook for MRI fast steps.
};

namespace SundialsUserFun {
    static int f (amrex::Real t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->f(t, y_data, y_rhs, user_data);
    }

    static int fi (amrex::Real t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->fi(t, y_data, y_rhs, user_data);
    }

    static int fe (amrex::Real t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->fe(t, y_data, y_rhs, user_data);
    }

    static int ff (amrex::Real t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->ff(t, y_data, y_rhs, user_data);
    }

    static int post_stage (amrex::Real t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->post_stage(t, y_data, user_data);
    }

    static int post_step (amrex::Real t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->post_step(t, y_data, user_data);
    }

    static int post_fast_stage (amrex::Real t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->post_fast_stage(t, y_data, user_data);
    }

    static int post_fast_step (amrex::Real t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->post_fast_step(t, y_data, user_data);
    }
}

/**
 * \brief IntegratorBase implementation powered by SUNDIALS ARKStep/MRIStep.
 *
 * \tparam T State container type (typically MultiFab) operated on by the user callbacks.
 */
template<class T>
class SundialsIntegrator : public IntegratorBase<T>
{
private:
    using BaseT = IntegratorBase<T>;

    // Method type: ERK, DIRK, IMEX-RK, EX-MRI, IM-MRI, IMEX-MRI
    std::string type = "ERK";

    // Use SUNDIALS default methods
    std::string method   = "DEFAULT"; // ERK, DIRK, or slow method with MRI
    std::string method_e = "DEFAULT"; // ERK in IMEX-RK
    std::string method_i = "DEFAULT"; // DIRK in IMEX-RK

    // Fast method type (ERK or DIRK) and method
    std::string fast_type   = "ERK";
    std::string fast_method = "DEFAULT";

    // Nonlinear solver
    std::string nonlinear_solver = "Newton";
    int max_nonlinear_iters = 0;

    std::string fast_nonlinear_solver = "Newton";
    int fast_max_nonlinear_iters = 0;

    // Linear solver
    std::string linear_solver = "GMRES";
    int max_linear_iters = 0;

    std::string fast_linear_solver = "GMRES";
    int fast_max_linear_iters = 0;

    // SUNDIALS package flags, set based on type
    bool use_ark = false;
    bool use_mri = false;

    // structure for interfacing with user-supplied functions
    SundialsUserData udata;

    // SUNDIALS context
    //
    // We should probably use context created by amrex:sundials::Initialize but
    // that context is not MPI-aware
    ::sundials::Context sunctx;

    // Single rate or slow time scale
    void *arkode_mem   = nullptr;
    SUNLinearSolver LS = nullptr;
    SUNNonlinearSolver NLS = nullptr;
    amrex::Real step_applied = 0;  // fixed step last given to SUNDIALS; 0: adaptive
    amrex::Real rtol_applied = -1; // tolerances last given to SUNDIALS; -1: none
    amrex::Real atol_applied = -1;

    // Fast time scale
    void *arkode_fast_mem            = nullptr;
    amrex::Real fast_step_applied    = 0; // fast fixed step last given to SUNDIALS; 0: adaptive
    amrex::Real fast_rtol_applied    = -1; // fast tolerances last given to SUNDIALS; -1: none
    amrex::Real fast_atol_applied    = -1;
    MRIStepInnerStepper fast_stepper = nullptr;
    SUNLinearSolver fast_LS          = nullptr;
    SUNNonlinearSolver fast_NLS      = nullptr;

    // Integrator stop time
    bool set_stop_time = false;
    amrex::Real stop_time = 0.0;

    // Max steps between returns
    amrex::Long max_num_steps = 0;

    //! Free the SUNDIALS objects created by initialize(). Safe to call twice.
    void free_memory ()
    {
        if (use_mri) {
            MRIStepInnerStepper_Free(&fast_stepper);
            AMREX_ARKSTEP(Free)(&arkode_fast_mem);
            AMREX_MRISTEP(Free)(&arkode_mem);
        } else if (use_ark) {
            AMREX_ARKSTEP(Free)(&arkode_mem);
        }
        SUNLinSolFree(LS);
        LS = nullptr;
        SUNLinSolFree(fast_LS);
        fast_LS = nullptr;
        SUNNonlinSolFree(NLS);
        NLS = nullptr;
        SUNNonlinSolFree(fast_NLS);
        fast_NLS = nullptr;
        use_ark = false;
        use_mri = false;
    }

    void initialize_parameters ()
    {
        amrex::ParmParse pp("integration.sundials");

        pp.query("type", type);
        pp.query("method", method);
        pp.query("method_e", method_e);
        pp.query("method_i", method_i);

        pp.query("fast_type", fast_type);
        pp.query("fast_method", fast_method);

        if (type == "ERK" || type == "DIRK" || type == "IMEX-RK") {
            use_ark = true;
        }
        else if (type == "EX-MRI" || type == "IM-MRI" || type == "IMEX-MRI") {
            use_mri = true;
        }
        else {
            std::string msg("Unknown method type: ");
            msg += type;
            amrex::Error(msg.c_str());
        }

        pp.query("nonlinear_solver", nonlinear_solver);
        pp.query("max_nonlinear_iters", max_nonlinear_iters);

        pp.query("fast_nonlinear_solver", fast_nonlinear_solver);
        pp.query("fast_max_nonlinear_iters", fast_max_nonlinear_iters);

        pp.query("linear_solver", linear_solver);
        pp.query("max_linear_iters", max_linear_iters);

        pp.query("fast_linear_solver", fast_linear_solver);
        pp.query("fast_max_linear_iters", fast_max_linear_iters);

        set_stop_time = pp.query("stop_time", stop_time);

        pp.query("max_num_steps", max_num_steps);
    }

    void SetupRK (amrex::Real time, N_Vector y_data)
    {
        step_applied = 0;
        rtol_applied = -1;
        atol_applied = -1;
        if (amrex::Verbose()) { amrex::Print() << "Using SUNDIALS time integrator\n"; }
        int flag = 0;

        // Create integrator and select method
        if (type == "ERK") {
            if (amrex::Verbose()) { amrex::Print() << "ERK method: " << method << "\n"; }
            arkode_mem = ARKStepCreate(SundialsUserFun::f, nullptr, time, y_data, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_mem != nullptr);
            if (method != "DEFAULT") {
                flag = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", method.c_str());
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }
        else if (type == "DIRK") {
            if (amrex::Verbose()) { amrex::Print() << "DIRK method: " << method << "\n"; }
            arkode_mem = ARKStepCreate(nullptr, SundialsUserFun::f, time, y_data, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_mem != nullptr);
            if (method != "DEFAULT") {
                flag = ARKStepSetTableName(arkode_mem, method.c_str(), "ARKODE_ERK_NONE");
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }
        else if (type == "IMEX-RK") {
            if (amrex::Verbose()) { amrex::Print() << "IMEX-RK method: " << method_i << " and "
                                    << method_e << "\n"; }
            arkode_mem = ARKStepCreate(SundialsUserFun::fe, SundialsUserFun::fi, time, y_data, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_mem != 0);
            if (method_e != "DEFAULT" && method_i != "DEFAULT")
            {
                flag = ARKStepSetTableName(arkode_mem, method_i.c_str(), method_e.c_str());
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Attach structure with user-supplied function wrappers
        flag = AMREX_ARKSTEP(SetUserData)(arkode_mem, &udata);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set integrator tolerances
        if (BaseT::use_adaptive_time_step || type == "DIRK" || type == "IMEX-RK") {
            if (amrex::Verbose()) {
                amrex::Print() << "Relative tolerance: " << BaseT::rel_tol << "\n";
                amrex::Print() << "Absolute tolerance: " << BaseT::abs_tol << "\n";
            }
            flag = AMREX_ARKSTEP(SStolerances)(arkode_mem, BaseT::rel_tol, BaseT::abs_tol);
            AMREX_ALWAYS_ASSERT(flag == 0);
            rtol_applied = BaseT::rel_tol;
            atol_applied = BaseT::abs_tol;
        }

        // Create and attach linear solver for implicit methods
        if (type == "DIRK" || type == "IMEX-RK") {
            if (amrex::Verbose()) {
                amrex::Print() << "Nonlinear solver: " << nonlinear_solver << "\n";
                amrex::Print() << "Max nonlinear iters: " << max_nonlinear_iters << "\n";
            }
            if (nonlinear_solver == "fixed-point") {
                NLS = SUNNonlinSol_FixedPoint(y_data, 0, sunctx);
                AMREX_ALWAYS_ASSERT(NLS != nullptr);
                flag = AMREX_ARKSTEP(SetNonlinearSolver)(arkode_mem, NLS);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
            flag = AMREX_ARKSTEP(SetMaxNonlinIters)(arkode_mem, max_nonlinear_iters);
            AMREX_ALWAYS_ASSERT(flag == 0);

            if (nonlinear_solver == "Newton") {
                if (amrex::Verbose()) {
                    amrex::Print() << "Linear solver: " << linear_solver << "\n";
                    amrex::Print() << "Max linear iters: " << max_linear_iters << "\n";
                }
                LS = SUNLinSol_SPGMR(y_data, SUN_PREC_NONE, max_linear_iters, sunctx);
                AMREX_ALWAYS_ASSERT(LS != nullptr);
                flag = AMREX_ARKSTEP(SetLinearSolver)(arkode_mem, LS, nullptr);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Set post stage and step function
        flag = AMREX_ARKSTEP(SetPostprocessStageFn)(arkode_mem, SundialsUserFun::post_stage);
        AMREX_ALWAYS_ASSERT(flag == 0);
        flag = AMREX_ARKSTEP(SetPostprocessStepFn)(arkode_mem, SundialsUserFun::post_step);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set a stop time
        if (set_stop_time) {
            if (amrex::Verbose()) { amrex::Print() << "Stop time: " << stop_time << "\n"; }
          flag = AMREX_ARKSTEP(SetStopTime)(arkode_mem, stop_time);
          AMREX_ALWAYS_ASSERT(flag == 0);
        }

        // Set max number of steps between returns
        flag = AMREX_ARKSTEP(SetMaxNumSteps)(arkode_mem, max_num_steps);
        AMREX_ALWAYS_ASSERT(flag == 0);
    }

    //! Pass slow tolerances to SUNDIALS when they are used and have changed.
    void apply_tolerances ()
    {
        if ((BaseT::use_adaptive_time_step || type == "DIRK" || type == "IMEX-RK" ||
             type == "IM-MRI" || type == "IMEX-MRI") &&
            (BaseT::rel_tol != rtol_applied || BaseT::abs_tol != atol_applied))
        {
#if AMREX_SUNDIALS_ARKODE_API
            int flag = ARKodeSStolerances(arkode_mem, BaseT::rel_tol, BaseT::abs_tol);
#else
            int flag = use_mri ? MRIStepSStolerances(arkode_mem, BaseT::rel_tol, BaseT::abs_tol)
                               : ARKStepSStolerances(arkode_mem, BaseT::rel_tol, BaseT::abs_tol);
#endif
            AMREX_ALWAYS_ASSERT(flag == 0);
            rtol_applied = BaseT::rel_tol;
            atol_applied = BaseT::abs_tol;
        }
    }

    //! Pass slow step size (0: adaptive) to SUNDIALS when it changes.
    void apply_step (amrex::Real h)
    {
        if (h != step_applied) {
#if AMREX_SUNDIALS_ARKODE_API
            int flag = ARKodeSetFixedStep(arkode_mem, h);
#else
            int flag = use_mri ? MRIStepSetFixedStep(arkode_mem, h)
                               : ARKStepSetFixedStep(arkode_mem, h);
#endif
            AMREX_ALWAYS_ASSERT(flag == 0);
            step_applied = h;
        }
    }

    //! Pass fast step size (0: adaptive) and tolerances to SUNDIALS when they change.
    void apply_fast_step ()
    {
        if ((BaseT::use_adaptive_fast_time_step || fast_type == "DIRK" || fast_type == "IMEX-RK") &&
            (BaseT::fast_rel_tol != fast_rtol_applied || BaseT::fast_abs_tol != fast_atol_applied))
        {
            int flag = AMREX_ARKSTEP(SStolerances)(arkode_fast_mem, BaseT::fast_rel_tol, BaseT::fast_abs_tol);
            AMREX_ALWAYS_ASSERT(flag == 0);
            fast_rtol_applied = BaseT::fast_rel_tol;
            fast_atol_applied = BaseT::fast_abs_tol;
        }
        amrex::Real h = BaseT::use_adaptive_fast_time_step ? amrex::Real(0) : BaseT::fast_time_step;
        if (h != fast_step_applied) {
            AMREX_ARKSTEP(SetFixedStep)(arkode_fast_mem, h);
            fast_step_applied = h;
        }
    }

    void SetupMRI (amrex::Real time, N_Vector y_data)
    {
        step_applied = 0;
        rtol_applied = -1;
        atol_applied = -1;
        fast_step_applied = 0;
        fast_rtol_applied = -1;
        fast_atol_applied = -1;
        if (amrex::Verbose()) { amrex::Print() << "Using SUNDIALS multirate time integrator\n"; }
        int flag = 0;

        // Create the fast integrator and select method
        if (fast_type == "ERK") {
            if (amrex::Verbose()) { amrex::Print() << "Fast ERK method: " << fast_method << "\n"; }
            arkode_fast_mem = ARKStepCreate(SundialsUserFun::ff, nullptr, time, y_data, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_fast_mem != nullptr);
            if (fast_method != "DEFAULT") {
                flag = ARKStepSetTableName(arkode_fast_mem, "ARKODE_DIRK_NONE", fast_method.c_str());
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }
        else if (fast_type == "DIRK") {
            if (amrex::Verbose()) { amrex::Print() << "Fast DIRK method: " << fast_method << "\n"; }
            arkode_fast_mem = ARKStepCreate(nullptr, SundialsUserFun::ff, time, y_data, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_fast_mem != nullptr);
            if (fast_method != "DEFAULT") {
                flag = ARKStepSetTableName(arkode_fast_mem, fast_method.c_str(), "ARKODE_ERK_NONE");
                AMREX_ALWAYS_ASSERT(flag == 0);
            }

            if (amrex::Verbose()) {
                amrex::Print() << "Fast nonlinear solver: " << fast_nonlinear_solver << "\n";
                amrex::Print() << "Fast max nonlinear iters: " << fast_max_nonlinear_iters << "\n";
            }
            if (fast_nonlinear_solver == "fixed-point") {
                fast_NLS = SUNNonlinSol_FixedPoint(y_data, 0, sunctx);
                AMREX_ALWAYS_ASSERT(fast_NLS != nullptr);
                flag = AMREX_ARKSTEP(SetNonlinearSolver)(arkode_fast_mem, fast_NLS);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
            flag = AMREX_ARKSTEP(SetMaxNonlinIters)(arkode_fast_mem, fast_max_nonlinear_iters);
            AMREX_ALWAYS_ASSERT(flag == 0);

            if (fast_nonlinear_solver == "Newton") {
                if (amrex::Verbose()) {
                    amrex::Print() << "Linear solver: " << fast_linear_solver << "\n";
                    amrex::Print() << "Max linear iters: " << fast_max_linear_iters << "\n";
                }
                fast_LS = SUNLinSol_SPGMR(y_data, SUN_PREC_NONE, fast_max_linear_iters, sunctx);
                AMREX_ALWAYS_ASSERT(fast_LS != nullptr);
                flag = AMREX_ARKSTEP(SetLinearSolver)(arkode_fast_mem, fast_LS, nullptr);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Attach structure with user-supplied function wrappers
        flag = AMREX_ARKSTEP(SetUserData)(arkode_fast_mem, &udata);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set integrator tolerances
        if (BaseT::use_adaptive_fast_time_step || fast_type == "DIRK" || fast_type == "IMEX-RK") {
            if (amrex::Verbose()) {
                amrex::Print() << "Fast relative tolerance: " << BaseT::fast_rel_tol << "\n";
                amrex::Print() << "Fast absolute tolerance: " << BaseT::fast_abs_tol << "\n";
            }
            flag = AMREX_ARKSTEP(SStolerances)(arkode_fast_mem, BaseT::fast_rel_tol, BaseT::fast_abs_tol);
            AMREX_ALWAYS_ASSERT(flag == 0);
            fast_rtol_applied = BaseT::fast_rel_tol;
            fast_atol_applied = BaseT::fast_abs_tol;
        }

        // Set post stage and step function
        flag = AMREX_ARKSTEP(SetPostprocessStageFn)(arkode_fast_mem, SundialsUserFun::post_fast_stage);
        AMREX_ALWAYS_ASSERT(flag == 0);
        flag = AMREX_ARKSTEP(SetPostprocessStepFn)(arkode_fast_mem, SundialsUserFun::post_fast_step);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set max number of steps between returns
        flag = AMREX_ARKSTEP(SetMaxNumSteps)(arkode_fast_mem, max_num_steps);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Wrap fast integrator as an inner stepper
#if defined(SUNDIALS_VERSION_MAJOR) && \
    ((SUNDIALS_VERSION_MAJOR > 7) || (SUNDIALS_VERSION_MAJOR == 7 && SUNDIALS_VERSION_MINOR >= 2))
        flag = ARKodeCreateMRIStepInnerStepper(arkode_fast_mem, &fast_stepper);
#else
        flag = ARKStepCreateMRIStepInnerStepper(arkode_fast_mem, &fast_stepper);
#endif
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Create slow integrator
        if (type == "EX-MRI") {
            if (amrex::Verbose()) { amrex::Print() << "EX-MRI method: " << method << "\n"; }
            arkode_mem = MRIStepCreate(SundialsUserFun::f, nullptr, time, y_data,
                                       fast_stepper, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_mem != nullptr);
        }
        else if (type == "IM-MRI") {
            if (amrex::Verbose()) { amrex::Print() << "IM-MRI method: " << method << "\n"; }
            arkode_mem = MRIStepCreate(nullptr, SundialsUserFun::f, time, y_data,
                                       fast_stepper, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_mem != nullptr);
        }
        else if (type == "IMEX-MRI") {
            if (amrex::Verbose()) { amrex::Print() << "IMEX-MRI method: " << method << "\n"; }
            arkode_mem = MRIStepCreate(SundialsUserFun::fe, SundialsUserFun::fi,
                                       time, y_data, fast_stepper, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_mem != nullptr);
        }

        // Set method
        if (method != "DEFAULT") {
            MRIStepCoupling MRIC = MRIStepCoupling_LoadTableByName(method.c_str());
            AMREX_ALWAYS_ASSERT(MRIC != nullptr);
            flag = MRIStepSetCoupling(arkode_mem, MRIC);
            AMREX_ALWAYS_ASSERT(flag == 0);
            MRIStepCoupling_Free(MRIC);
        }

        // Attach structure with user-supplied function wrappers
        flag = AMREX_MRISTEP(SetUserData)(arkode_mem, &udata);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set integrator tolerances
        if (BaseT::use_adaptive_time_step || type == "IM-MRI" || type == "IMEX-MRI") {
            if (amrex::Verbose()) {
                amrex::Print() << "Relative tolerance: " << BaseT::rel_tol << "\n";
                amrex::Print() << "Absolute tolerance: " << BaseT::abs_tol << "\n";
            }
            flag = AMREX_MRISTEP(SStolerances)(arkode_mem, BaseT::rel_tol, BaseT::abs_tol);
            AMREX_ALWAYS_ASSERT(flag == 0);
            rtol_applied = BaseT::rel_tol;
            atol_applied = BaseT::abs_tol;
        }

        // Create and attach linear solver
        if (type == "IM-MRI" || type == "IMEX-MRI") {
            if (amrex::Verbose()) {
                amrex::Print() << "Nonlinear solver: " << nonlinear_solver << "\n";
                amrex::Print() << "Max nonlinear iters: " << max_nonlinear_iters << "\n";
            }
            if (nonlinear_solver == "fixed-point") {
                NLS = SUNNonlinSol_FixedPoint(y_data, 0, sunctx);
                AMREX_ALWAYS_ASSERT(NLS != nullptr);
                flag = AMREX_MRISTEP(SetNonlinearSolver)(arkode_mem, NLS);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
            flag = AMREX_MRISTEP(SetMaxNonlinIters)(arkode_mem, max_nonlinear_iters);
            AMREX_ALWAYS_ASSERT(flag == 0);

            if (nonlinear_solver == "Newton") {
                if (amrex::Verbose()) {
                    amrex::Print() << "Linear solver: " << linear_solver << "\n";
                    amrex::Print() << "Max linear iters: " << max_linear_iters << "\n";
                }
                LS = SUNLinSol_SPGMR(y_data, SUN_PREC_NONE, max_linear_iters, sunctx);
                AMREX_ALWAYS_ASSERT(LS != nullptr);
                flag = AMREX_MRISTEP(SetLinearSolver)(arkode_mem, LS, nullptr);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Set post stage and step function
        flag = AMREX_MRISTEP(SetPostprocessStageFn)(arkode_mem, SundialsUserFun::post_stage);
        AMREX_ALWAYS_ASSERT(flag == 0);
        flag = AMREX_MRISTEP(SetPostprocessStepFn)(arkode_mem, SundialsUserFun::post_step);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set a stop time
        if (set_stop_time) {
            if (amrex::Verbose()) { amrex::Print() << "Stop time: " << stop_time << "\n"; }
          flag = AMREX_MRISTEP(SetStopTime)(arkode_mem, stop_time);
          AMREX_ALWAYS_ASSERT(flag == 0);
        }

        // Set max number of steps between returns
        flag = AMREX_MRISTEP(SetMaxNumSteps)(arkode_mem, max_num_steps);
        AMREX_ALWAYS_ASSERT(flag == 0);
    }

    // -------------------------------------
    // Vector<MultiFab> / N_Vector Utilities
    // -------------------------------------

    // Utility to unpack a SUNDIALS ManyVector into a vector of MultiFabs
    void unpack_vector (N_Vector y_data, amrex::Vector<amrex::MultiFab>& S_data)
    {
        const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
        S_data.resize(num_vecs);

        for(int i = 0; i < num_vecs; i++)
        {
            S_data.at(i) = amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),
                                           amrex::make_alias,
                                           0,
                                           amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
        }
    }

    // Utility to wrap vector of MultiFabs as a SUNDIALS ManyVector
    N_Vector wrap_data (amrex::Vector<amrex::MultiFab>& S_data)
    {
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_data[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        sunindextype NV_len = S_data.size();
        N_Vector* NV_array  = new N_Vector[NV_len];

        for (int i = 0; i < NV_len; ++i) {
            NV_array[i] = amrex::sundials::N_VMake_MultiFab(get_length(i),
                                                            &S_data[i], &sunctx);
        }

        N_Vector y_data = N_VNew_ManyVector(NV_len, NV_array, sunctx);

        delete[] NV_array;

        return y_data;
    }

    // Utility to wrap vector of MultiFabs as a SUNDIALS ManyVector
    N_Vector copy_data (const amrex::Vector<amrex::MultiFab>& S_data)
    {
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_data[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        sunindextype NV_len = S_data.size();
        N_Vector* NV_array  = new N_Vector[NV_len];

        for (int i = 0; i < NV_len; ++i) {
            NV_array[i] = amrex::sundials::N_VNew_MultiFab(get_length(i),
                                                           S_data[i].boxArray(),
                                                           S_data[i].DistributionMap(),
                                                           S_data[i].nComp(),
                                                           S_data[i].nGrow(),
                                                           &sunctx);

            MultiFab::Copy(*amrex::sundials::getMFptr(NV_array[i]),
                           S_data[i],
                           0,
                           0,
                           S_data[i].nComp(),
                           S_data[i].nGrow());
        }

        N_Vector y_data =  N_VNew_ManyVector(NV_len, NV_array, sunctx);

        delete[] NV_array;

        return y_data;
    }

    // -----------------------------
    // MultiFab / N_Vector Utilities
    // -----------------------------

    // Utility to unpack a SUNDIALS Vector into a MultiFab
    void unpack_vector (N_Vector y_data, amrex::MultiFab& S_data)
    {
        S_data = amrex::MultiFab(*amrex::sundials::getMFptr(y_data),
                                 amrex::make_alias,
                                 0,
                                 amrex::sundials::getMFptr(y_data)->nComp());
    }

    // Utility to wrap a MultiFab as a SUNDIALS Vector
    N_Vector wrap_data (amrex::MultiFab& S_data)
    {
        return amrex::sundials::N_VMake_MultiFab(S_data.nComp() * S_data.boxArray().numPts(),
                                                 &S_data, &sunctx);
    }

    // Utility to wrap a MultiFab as a SUNDIALS Vector
    N_Vector copy_data (const amrex::MultiFab& S_data)
    {
        N_Vector y_data = amrex::sundials::N_VNew_MultiFab(S_data.nComp() * S_data.boxArray().numPts(),
                                                           S_data.boxArray(),
                                                           S_data.DistributionMap(),
                                                           S_data.nComp(),
                                                           S_data.nGrow(),
                                                           &sunctx);

        MultiFab::Copy(*amrex::sundials::getMFptr(y_data),
                       S_data,
                       0,
                       0,
                       S_data.nComp(),
                       S_data.nGrow());

        return y_data;
    }

public:
    /**
     * \brief Construct an uninitialized integrator; call initialize() before use.
     */
    SundialsIntegrator () {}

    /**
     * \brief Construct and immediately configure the integrator with \p S_data at time \p time.
     *
     * \param S_data Initial state used to size and seed the solver.
     * \param time   Physical time associated with \p S_data.
     */
    SundialsIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize(S_data, time);
    }

    /**
     * \brief Configure (or reconfigure) the SUNDIALS integrator for the provided state.
     *
     * \param S_data State snapshot used to size vector storage.
     * \param time   Start time associated with \p S_data.
     */
    void initialize (const T& S_data, const amrex::Real time = 0.0)
    {
        free_memory(); // no-op on the first call
        {
            // Destroy the current context. Move assignment below does not free
            // it in SUNDIALS older than 7.6.
            ::sundials::Context old_sunctx = std::move(sunctx);
        }

        initialize_parameters();
#ifdef AMREX_USE_MPI
        MPI_Comm mpi_comm = ParallelContext::CommunicatorSub();
#endif
#if defined(SUNDIALS_VERSION_MAJOR) && (SUNDIALS_VERSION_MAJOR < 7)
#  ifdef AMREX_USE_MPI
        sunctx = ::sundials::Context(&mpi_comm);
#  else
        sunctx = ::sundials::Context(nullptr);
#  endif
#else
#  ifdef AMREX_USE_MPI
        sunctx = ::sundials::Context(mpi_comm);
#  else
        sunctx = ::sundials::Context(SUN_COMM_NULL);
#  endif
#endif

        // Right-hand side function wrappers
        udata.f = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                      void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::Rhs(S_rhs, state, rhs_time);

            return 0;
        };

        udata.fi = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                       void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::RhsIm(S_rhs, state, rhs_time);

            return 0;
        };

        udata.fe = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                       void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::RhsEx(S_rhs, state, rhs_time);

            return 0;
        };

        udata.ff = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                       void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::RhsFast(S_rhs, state, rhs_time);

            return 0;
        };

        udata.post_stage = [&](amrex::Real t, N_Vector y_data,
                               void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            BaseT::post_stage_action(state, t);

            return 0;
        };

        udata.post_step = [&](amrex::Real t, N_Vector y_data,
                              void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            BaseT::post_step_action(state, t);

            return 0;
        };

        udata.post_fast_stage = [&](amrex::Real t, N_Vector y_data,
                                    void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            BaseT::post_fast_stage_action(state, t);

            return 0;
        };

        udata.post_fast_step = [&](amrex::Real t, N_Vector y_data,
                                   void * /* user_data */) -> int {

            T state;
            unpack_vector(y_data, state);

            BaseT::post_fast_step_action(state, t);

            return 0;
        };

        N_Vector y_data = copy_data(S_data); // ideally just wrap and ignore const

        if (use_ark) {
            SetupRK(time, y_data);
        }
        else if (use_mri)
        {
            SetupMRI(time, y_data);
        }

        N_VDestroy(y_data);
    }

    /**
     * \brief Destroy the integrator, printing summary statistics when verbose.
     */
    virtual ~SundialsIntegrator () {
        // Print integrator statistics
        if (amrex::Verbose()) {
            if (type == "EX-MRI" || type == "IM-MRI" || type == "IMEX-MRI") {
                amrex::Print() << "Slow Time Integrator Stats\n";
                if (ParallelDescriptor::IOProcessor()) {
                    AMREX_MRISTEP(PrintAllStats)(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
                }
                amrex::Print() << "Fast Time Integrator Stats\n";
                if (ParallelDescriptor::IOProcessor()) {
                    AMREX_ARKSTEP(PrintAllStats)(arkode_fast_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
                }
            } else {
                amrex::Print() << "Time Integrator Stats\n";
                if (ParallelDescriptor::IOProcessor()) {
                    AMREX_ARKSTEP(PrintAllStats)(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
                }
            }
        }

        // Clean up allocated memory
        free_memory();
    }

    /**
     * \brief Take a single time step of size \p dt starting from \p S_old.
     *
     * \param S_old State at the beginning of the step (wrapped and passed to SUNDIALS).
     * \param S_new Output state receiving the step result.
     * \param time  Beginning-of-step time.
     * \param dt    Requested time step.
     * \return The step size actually taken (currently \p dt).
     */
    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real dt) override
    {
        amrex::Real tout = time + dt;
        amrex::Real tret;

        N_Vector y_old = wrap_data(S_old);
        N_Vector y_new = wrap_data(S_new);

        if (use_ark) {
            AMREX_ARKSTEP(Reset)(arkode_mem, time, y_old); // should probably resize
            apply_tolerances();
            AMREX_ARKSTEP(SetFixedStep)(arkode_mem, dt);
            step_applied = dt;
            int flag = AMREX_ARKSTEP(Evolve)(arkode_mem, tout, y_new, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
            AMREX_MRISTEP(Reset)(arkode_mem, time, y_old); // should probably resize -- need to resize inner stepper
            apply_tolerances();
            AMREX_MRISTEP(SetFixedStep)(arkode_mem, dt);
            step_applied = dt;
            apply_fast_step();
            int flag = AMREX_MRISTEP(Evolve)(arkode_mem, tout, y_new, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        } else {
            Error("SUNDIALS integrator type not specified.");
        }

        N_VDestroy(y_old);
        N_VDestroy(y_new);

        return dt;
    }

    /**
     * \brief Evolve the solution in \p S_out up to \p time_out using ARKStep/MRIStep.
     *
     * \param S_out    State providing the initial condition and receiving the final solution.
     * \param time_out Target time for the evolution.
     */
    void evolve (T& S_out, const amrex::Real time_out) override
    {
        int flag = 0;         // SUNDIALS return status
        amrex::Real time_ret; // SUNDIALS return time

        N_Vector y_out = wrap_data(S_out);

        if (use_ark) {
            apply_tolerances();
            apply_step(BaseT::use_adaptive_time_step ? amrex::Real(0) : BaseT::time_step);
            flag = AMREX_ARKSTEP(Evolve)(arkode_mem, time_out, y_out, &time_ret, ARK_NORMAL);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
#if defined(SUNDIALS_VERSION_MAJOR) && \
    ((SUNDIALS_VERSION_MAJOR < 7) || (SUNDIALS_VERSION_MAJOR == 7 && SUNDIALS_VERSION_MINOR < 2))
            if (BaseT::use_adaptive_time_step && step_applied == 0) {
                amrex::Abort("SundialsIntegrator: adaptive slow time steps with MRI methods need SUNDIALS 7.2 or later");
            }
#endif
            apply_tolerances();
            // Do not switch a fixed slow step back to adaptive; the fixed-step MRI table may lack an embedding.
            if (!BaseT::use_adaptive_time_step) {
                apply_step(BaseT::time_step);
            }
            apply_fast_step();
            flag = AMREX_MRISTEP(Evolve)(arkode_mem, time_out, y_out, &time_ret, ARK_NORMAL);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        } else {
            Error("SUNDIALS integrator type not specified.");
        }

        N_VDestroy(y_out);
    }

    /**
     * \brief Interpolate between SUNDIALS stages (not yet implemented for this integrator).
     */
    void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override {}

    /**
     * \brief Apply a user-supplied mapping to every MultiFab in the integrator (unused placeholder).
     */
    void map_data (std::function<void(T&)> /* Map */) override {}
};

}

#undef AMREX_SUNDIALS_ARKODE_API
#undef AMREX_ARKSTEP
#undef AMREX_MRISTEP

#endif
