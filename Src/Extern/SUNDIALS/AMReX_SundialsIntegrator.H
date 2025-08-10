#ifndef AMREX_SUNDIALS_INTEGRATOR_H
#define AMREX_SUNDIALS_INTEGRATOR_H

#include <functional>

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

namespace amrex {

struct SundialsUserData {
    // ERK or DIRK right-hand side function
    // EX-MRI or IM-MRI slow right-hand side function
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> f;

    // ImEx-RK right-hand side functions
    // ImEx-MRI slow right-hand side functions
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> fi;
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> fe;

    // MRI fast time scale right-hand side function
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> ff;

    // Post stage and step actions
    std::function<int(amrex::Real, N_Vector, void*)> post_stage;
    std::function<int(amrex::Real, N_Vector, void*)> post_step;

    // Post fast stage and step actions
    std::function<int(amrex::Real, N_Vector, void*)> post_fast_stage;
    std::function<int(amrex::Real, N_Vector, void*)> post_fast_step;
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

template<class T>
class SundialsIntegrator : public IntegratorBase<T>
{
private:
    using BaseT = IntegratorBase<T>;

    // Method type: ERK, DIRK, IMEX-RK, EX-MRI, IM-MRI, IMEX-MRI
    std::string type = "ERK";

    // Use SUNDIALS default methods
    std::string method   = "DEFAULT";
    std::string method_e = "DEFAULT";
    std::string method_i = "DEFAULT";

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

    // Fast time scale
    void *arkode_fast_mem            = nullptr;
    MRIStepInnerStepper fast_stepper = nullptr;
    SUNLinearSolver fast_LS          = nullptr;
    SUNNonlinearSolver fast_NLS      = nullptr;

    // Integrator stop time
    bool set_stop_time = false;
    amrex::Real stop_time = 0.0;

    // Max steps between returns
    amrex::Long max_num_steps = 0;

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
        flag = ARKStepSetUserData(arkode_mem, &udata);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set integrator tolerances
        if (BaseT::use_adaptive_time_step || type == "DIRK" || type == "IMEX-RK") {
            if (amrex::Verbose()) {
                amrex::Print() << "Relative tolerance: " << BaseT::rel_tol << "\n";
                amrex::Print() << "Absolute tolerance: " << BaseT::abs_tol << "\n";
            }
            flag = ARKStepSStolerances(arkode_mem, BaseT::rel_tol, BaseT::abs_tol);
            AMREX_ALWAYS_ASSERT(flag == 0);
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
                flag = ARKStepSetNonlinearSolver(arkode_mem, NLS);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
            flag = ARKStepSetMaxNonlinIters(arkode_mem, max_nonlinear_iters);
            AMREX_ALWAYS_ASSERT(flag == 0);

            if (nonlinear_solver == "Newton") {
                if (amrex::Verbose()) {
                    amrex::Print() << "Linear solver: " << linear_solver << "\n";
                    amrex::Print() << "Max linear iters: " << max_linear_iters << "\n";
                }
                LS = SUNLinSol_SPGMR(y_data, SUN_PREC_NONE, max_linear_iters, sunctx);
                AMREX_ALWAYS_ASSERT(LS != nullptr);
                flag = ARKStepSetLinearSolver(arkode_mem, LS, nullptr);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Set post stage and step function
        flag = ARKStepSetPostprocessStageFn(arkode_mem, SundialsUserFun::post_stage);
        AMREX_ALWAYS_ASSERT(flag == 0);
        flag = ARKStepSetPostprocessStepFn(arkode_mem, SundialsUserFun::post_step);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set a stop time
        if (set_stop_time) {
            if (amrex::Verbose()) { amrex::Print() << "Stop time: " << stop_time << "\n"; }
          flag = ARKStepSetStopTime(arkode_mem, stop_time);
          AMREX_ALWAYS_ASSERT(flag == 0);
        }

        // Set max number of steps between returns
        ARKStepSetMaxNumSteps(arkode_mem, max_num_steps);
        AMREX_ALWAYS_ASSERT(flag == 0);
    }

    void SetupMRI (amrex::Real time, N_Vector y_data)
    {
        if (amrex::Verbose()) { amrex::Print() << "Using SUNDIALS multirate time integrator\n"; }
        int flag = 0;

        // Create the fast integrator and select method
        if (fast_type == "ERK") {
            if (amrex::Verbose()) { amrex::Print() << "Fast ERK method: " << method << "\n"; }
            arkode_fast_mem = ARKStepCreate(SundialsUserFun::ff, nullptr, time, y_data, sunctx);
            AMREX_ALWAYS_ASSERT(arkode_fast_mem != nullptr);
            if (fast_method != "DEFAULT") {
                flag = ARKStepSetTableName(arkode_fast_mem, "ARKODE_DIRK_NONE", fast_method.c_str());
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }
        else if (fast_type == "DIRK") {
            if (amrex::Verbose()) { amrex::Print() << "Fast DIRK method: " << method << "\n"; }
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
                flag = ARKStepSetNonlinearSolver(arkode_mem, fast_NLS);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
            flag = ARKStepSetMaxNonlinIters(arkode_mem, fast_max_nonlinear_iters);
            AMREX_ALWAYS_ASSERT(flag == 0);

            if (fast_nonlinear_solver == "Newton") {
                if (amrex::Verbose()) {
                    amrex::Print() << "Linear solver: " << fast_linear_solver << "\n";
                    amrex::Print() << "Max linear iters: " << fast_max_linear_iters << "\n";
                }
                fast_LS = SUNLinSol_SPGMR(y_data, SUN_PREC_NONE, fast_max_linear_iters, sunctx);
                AMREX_ALWAYS_ASSERT(fast_LS != nullptr);
                flag = ARKStepSetLinearSolver(arkode_mem, fast_LS, nullptr);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Attach structure with user-supplied function wrappers
        flag = ARKStepSetUserData(arkode_fast_mem, &udata);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set integrator tolerances
        if (BaseT::use_adaptive_fast_time_step || fast_type == "DIRK" || fast_type == "IMEX-RK") {
            if (amrex::Verbose()) {
                amrex::Print() << "Fast relative tolerance: " << BaseT::fast_rel_tol << "\n";
                amrex::Print() << "Fast absolute tolerance: " << BaseT::fast_abs_tol << "\n";
            }
          flag = ARKStepSStolerances(arkode_fast_mem, BaseT::fast_rel_tol, BaseT::fast_abs_tol);
          AMREX_ALWAYS_ASSERT(flag == 0);
        }

        // Set post stage and step function
        flag = ARKStepSetPostprocessStageFn(arkode_fast_mem, SundialsUserFun::post_fast_stage);
        AMREX_ALWAYS_ASSERT(flag == 0);
        flag = ARKStepSetPostprocessStepFn(arkode_fast_mem, SundialsUserFun::post_fast_step);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set max number of steps between returns
        ARKStepSetMaxNumSteps(arkode_fast_mem, max_num_steps);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Wrap fast integrator as an inner stepper
        flag = ARKStepCreateMRIStepInnerStepper(arkode_fast_mem, &fast_stepper);
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
        flag = MRIStepSetUserData(arkode_mem, &udata);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set integrator tolerances
        if (BaseT::use_adaptive_time_step || fast_type == "IM-MRI" || fast_type == "IMEX-MRI") {
            if (amrex::Verbose()) {
                amrex::Print() << "Relative tolerance: " << BaseT::rel_tol << "\n";
                amrex::Print() << "Absolute tolerance: " << BaseT::abs_tol << "\n";
            }
            flag = MRIStepSStolerances(arkode_mem, BaseT::rel_tol, BaseT::abs_tol);
            AMREX_ALWAYS_ASSERT(flag == 0);
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
                flag = ARKStepSetNonlinearSolver(arkode_mem, NLS);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
            flag = ARKStepSetMaxNonlinIters(arkode_mem, max_nonlinear_iters);
            AMREX_ALWAYS_ASSERT(flag == 0);

            if (nonlinear_solver == "Newton") {
                if (amrex::Verbose()) {
                    amrex::Print() << "Linear solver: " << linear_solver << "\n";
                    amrex::Print() << "Max linear iters: " << max_linear_iters << "\n";
                }
                LS = SUNLinSol_SPGMR(y_data, SUN_PREC_NONE, max_linear_iters, sunctx);
                AMREX_ALWAYS_ASSERT(LS != nullptr);
                flag = ARKStepSetLinearSolver(arkode_mem, LS, nullptr);
                AMREX_ALWAYS_ASSERT(flag == 0);
            }
        }

        // Set post stage and step function
        flag = MRIStepSetPostprocessStageFn(arkode_mem, SundialsUserFun::post_stage);
        AMREX_ALWAYS_ASSERT(flag == 0);
        flag = MRIStepSetPostprocessStepFn(arkode_mem, SundialsUserFun::post_step);
        AMREX_ALWAYS_ASSERT(flag == 0);

        // Set a stop time
        if (set_stop_time) {
            if (amrex::Verbose()) { amrex::Print() << "Stop time: " << stop_time << "\n"; }
          flag = ARKStepSetStopTime(arkode_mem, stop_time);
          AMREX_ALWAYS_ASSERT(flag == 0);
        }

        // Set max number of steps between returns
        ARKStepSetMaxNumSteps(arkode_mem, max_num_steps);
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
    };

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
    };

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
    };

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
    };

    // Utility to wrap a MultiFab as a SUNDIALS Vector
    N_Vector wrap_data (amrex::MultiFab& S_data)
    {
        return amrex::sundials::N_VMake_MultiFab(S_data.nComp() * S_data.boxArray().numPts(),
                                                 &S_data, &sunctx);
    };

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
    };

public:
    SundialsIntegrator () {}

    SundialsIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize(S_data, time);
    }

    void initialize (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize_parameters();
        MPI_Comm mpi_comm = ParallelContext::CommunicatorSub();
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

            T S_data;
            unpack_vector(y_data, S_data);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::Rhs(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.fi = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                       void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::RhsIm(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.fe = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                       void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::RhsEx(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.ff = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs,
                       void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            T S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::RhsFast(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.post_stage = [&](amrex::Real time, N_Vector y_data,
                               void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            BaseT::post_stage_action(S_data, time);

            return 0;
        };

        udata.post_step = [&](amrex::Real time, N_Vector y_data,
                              void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            BaseT::post_step_action(S_data, time);

            return 0;
        };

        udata.post_fast_stage = [&](amrex::Real time, N_Vector y_data,
                                    void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            BaseT::post_fast_stage_action(S_data, time);

            return 0;
        };

        udata.post_fast_step = [&](amrex::Real time, N_Vector y_data,
                                   void * /* user_data */) -> int {

            T S_data;
            unpack_vector(y_data, S_data);

            BaseT::post_fast_step_action(S_data, time);

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

    virtual ~SundialsIntegrator () {
        // Print integrator statistics
        if (amrex::Verbose()) {
            if (type == "EX-MRI" || type == "IM-MRI" || type == "IMEX-MRI") {
                amrex::Print() << "Slow Time Integrator Stats\n";
                if (ParallelDescriptor::IOProcessor()) {
                    MRIStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
                }
                amrex::Print() << "Fast Time Integrator Stats\n";
                if (ParallelDescriptor::IOProcessor()) {
                    ARKStepPrintAllStats(arkode_fast_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
                }
            } else {
                amrex::Print() << "Time Integrator Stats\n";
                if (ParallelDescriptor::IOProcessor()) {
                    ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
                }
            }
        }

        // Clean up allocated memory
        SUNLinSolFree(LS);
        SUNLinSolFree(fast_LS);
        SUNNonlinSolFree(NLS);
        SUNNonlinSolFree(fast_NLS);
        MRIStepInnerStepper_Free(&fast_stepper);
        MRIStepFree(&arkode_fast_mem);
        ARKStepFree(&arkode_mem);
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real dt) override
    {
        amrex::Real tout = time + dt;
        amrex::Real tret;

        N_Vector y_old = wrap_data(S_old);
        N_Vector y_new = wrap_data(S_new);

        if (use_ark) {
            ARKStepReset(arkode_mem, time, y_old); // should probably resize
            ARKStepSetFixedStep(arkode_mem, dt);
            int flag = ARKStepEvolve(arkode_mem, tout, y_new, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
            MRIStepReset(arkode_mem, time, y_old); // should probably resize -- need to resize inner stepper
            MRIStepSetFixedStep(arkode_mem, dt);
            int flag = MRIStepEvolve(arkode_mem, tout, y_new, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        } else {
            Error("SUNDIALS integrator type not specified.");
        }

        N_VDestroy(y_old);
        N_VDestroy(y_new);

        return dt;
    }

    void evolve (T& S_out, const amrex::Real time_out) override
    {
        int flag = 0;         // SUNDIALS return status
        amrex::Real time_ret; // SUNDIALS return time

        N_Vector y_out = wrap_data(S_out);

        if (use_ark) {
            if (!BaseT::use_adaptive_time_step) {
                ARKStepSetFixedStep(arkode_mem, BaseT::time_step);
            }
            flag = ARKStepEvolve(arkode_mem, time_out, y_out, &time_ret, ARK_NORMAL);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
            if (!BaseT::use_adaptive_time_step) {
                MRIStepSetFixedStep(arkode_mem, BaseT::time_step);
            }
            if (!BaseT::use_adaptive_fast_time_step) {
                ARKStepSetFixedStep(arkode_fast_mem, BaseT::fast_time_step);
            }
            flag = MRIStepEvolve(arkode_mem, time_out, y_out, &time_ret, ARK_NORMAL);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        } else {
            Error("SUNDIALS integrator type not specified.");
        }

        N_VDestroy(y_out);
    }

    void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override {}

    void map_data (std::function<void(T&)> /* Map */) override {}
};

}

#endif
