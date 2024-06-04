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
#include <sunlinsol/sunlinsol_spgmr.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>

namespace amrex {

struct SundialsUserData {
    // ERK or DIRK right-hand side function
    // ExMRI or ImMRI slow right-hand side function
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> f;

    // ImEx-RK right-hand side functions
    // ImEx-MRI slow right-hand side functions
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> fi;
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> fe;

    // MRI fast time scale right-hand side function
    std::function<int(amrex::Real, N_Vector, N_Vector, void*)> ff;

    // Post step actions
    std::function<int(amrex::Real, N_Vector, void*)> post_step;
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

    static int post_step (amrex::Real t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->post_step(t, y_data, user_data);
    }
}

template<class T>
class SundialsIntegrator : public IntegratorBase<T>
{
private:
    using BaseT = IntegratorBase<T>;

    // method type and name
    std::string type     = "ERK";
    std::string method   = "DEFAULT";
    std::string method_e = "DEFAULT";
    std::string method_i = "DEFAULT";

    std::string fast_type   = "ERK";
    std::string fast_method = "DEFAULT";

    // SUNDIALS package flags, set based on type
    bool use_ark = false;
    bool use_mri = false;

    // structure for interfacing with user-supplied functions
    SundialsUserData udata;

    // SUNDIALS objects
    SUNContext sunctx = nullptr;

    void *arkode_mem   = nullptr;
    SUNLinearSolver LS = nullptr;

    void *arkode_fast_mem            = nullptr;
    MRIStepInnerStepper fast_stepper = nullptr;
    SUNLinearSolver fast_LS          = nullptr;

    // relative and absolute tolerances
    Real reltol = 1.0e-4;
    Real abstol = 1.0e-9;

    Real fast_reltol = 1.0e-4;
    Real fast_abstol = 1.0e-9;

    int NVar;               // NOTE: expects S_data to be a Vector<MultiFab>
    N_Vector* nv_many_arr;  /* vector array composed of cons, xmom, ymom, zmom component vectors */
    N_Vector nv_S;

    void initialize_parameters ()
    {
        amrex::ParmParse pp("integration.sundials");

        pp.query("type", type);
        pp.query("method", method);
        pp.query("method_e", method);
        pp.query("method_i", method);

        pp.query("fast_type", fast_type);
        pp.query("fast_method", fast_method);

        if (type == "ERK" || type == "DIRK" || type == "IMEX-RK") {
            use_ark = true;
        }
        else if (type == "EX-MRI" || type == "IM-MRI" || type == "IMEX-MRI") {
            use_mri = true;
        }
        else {
            std::string msg("Unknown strategy: ");
            msg += type;
            amrex::Error(msg.c_str());
        }
    }

    void SetupRK(amrex::Real time, N_Vector y_data)
    {
        // Create integrator and select method
        if (type == "ERK") {
            amrex::Print() << "SUNDIALS ERK time integrator\n";
            arkode_mem = ARKStepCreate(SundialsUserFun::f, nullptr, time, nv_S, sunctx);

            if (method != "DEFAULT") {
                amrex::Print() << "SUNDIALS ERK method " << method << "\n";
                ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", method.c_str());
            }
        }
        else if (type == "DIRK") {
            amrex::Print() << "SUNDIALS DIRK time integrator\n";
            arkode_mem = ARKStepCreate(nullptr, SundialsUserFun::f, time, nv_S, sunctx);

            if (method != "DEFAULT") {
                amrex::Print() << "SUNDIALS DIRK method " << method << "\n";
                ARKStepSetTableName(arkode_mem, method.c_str(), "ARKODE_ERK_NONE");
            }
        }
        else if (type == "IMEX-RK") {
            amrex::Print() << "SUNDIALS IMEX time integrator\n";
            arkode_mem = ARKStepCreate(SundialsUserFun::fe, SundialsUserFun::fi, time, nv_S, sunctx);

            if (method_e != "DEFAULT" && method_i != "DEFAULT")
            {
                amrex::Print() << "SUNDIALS IMEX method " << method_i << " and "
                               << method_e << "\n";
                ARKStepSetTableName(arkode_mem, method_i.c_str(), method_e.c_str());
            }
        }

        // Attach structure with user-supplied function wrappers
        ARKStepSetUserData(arkode_mem, &udata);

        // Set integrator tolerances
        ARKStepSStolerances(arkode_mem, reltol, abstol);

        // Create and attach linear solver for implicit methods
        if (type == "DIRK" || type == "IMEX-RK") {
            LS = SUNLinSol_SPGMR(nv_S, SUN_PREC_NONE, 0, sunctx);
            ARKStepSetLinearSolver(arkode_mem, LS, nullptr);
        }

        // Set post step function
        ARKStepSetPostprocessStepFn(arkode_mem, SundialsUserFun::post_step);
    }

    void SetupMRI(amrex::Real time, N_Vector y_data)
    {
        // Create the fast integrator and select method
        if (fast_type == "ERK") {
            amrex::Print() << "SUNDIALS ERK time integrator\n";
            arkode_fast_mem = ARKStepCreate(SundialsUserFun::ff, nullptr, time, nv_S, sunctx);

            if (method != "DEFAULT") {
                amrex::Print() << "SUNDIALS ERK method " << method << "\n";
                ARKStepSetTableName(arkode_fast_mem, "ARKODE_DIRK_NONE", fast_method.c_str());
            }
        }
        else if (fast_type == "DIRK") {
            amrex::Print() << "SUNDIALS DIRK time integrator\n";
            arkode_fast_mem = ARKStepCreate(nullptr, SundialsUserFun::ff, time, nv_S, sunctx);

            if (method != "DEFAULT") {
                amrex::Print() << "SUNDIALS DIRK method " << method << "\n";
                ARKStepSetTableName(arkode_fast_mem, fast_method.c_str(), "ARKODE_ERK_NONE");
            }

            fast_LS = SUNLinSol_SPGMR(nv_S, SUN_PREC_NONE, 0, sunctx);
            ARKStepSetLinearSolver(arkode_fast_mem, fast_LS, nullptr);
        }

        // Attach structure with user-supplied function wrappers
        ARKStepSetUserData(arkode_fast_mem, &udata);

        // Set integrator tolerances
        ARKStepSStolerances(arkode_fast_mem, fast_reltol, fast_abstol);

        // Wrap fast integrator as an inner stepper
        ARKStepCreateMRIStepInnerStepper(arkode_fast_mem, &fast_stepper);

        // Create slow integrator
        if (type == "EX-MRI") {
            amrex::Print() << "SUNDIALS ERK time integrator\n";
            arkode_mem = MRIStepCreate(SundialsUserFun::f, nullptr, time, nv_S,
                                       fast_stepper, sunctx);
        }
        else if (type == "IM-MRI") {
            amrex::Print() << "SUNDIALS DIRK time integrator\n";
            arkode_mem = MRIStepCreate(nullptr, SundialsUserFun::f, time, nv_S,
                                       fast_stepper, sunctx);
        }
        else if (type == "IMEX-MRI") {
            amrex::Print() << "SUNDIALS IMEX time integrator\n";
            arkode_mem = MRIStepCreate(SundialsUserFun::fe, SundialsUserFun::fi,
                                       time, nv_S, fast_stepper, sunctx);
        }

        // Set method
        if (method != "DEFAULT") {
            MRIStepCoupling MRIC = MRIStepCoupling_LoadTableByName(method.c_str());
            MRIStepSetCoupling(arkode_mem, MRIC);
            MRIStepCoupling_Free(MRIC);
        }

        // Attach structure with user-supplied function wrappers
        MRIStepSetUserData(arkode_mem, &udata);

        // Set integrator tolerances
        MRIStepSStolerances(arkode_mem, reltol, abstol);

        // Create and attach linear solver
        if (type == "IM-MRI" || type == "IMEX-MRI") {
            LS = SUNLinSol_SPGMR(nv_S, SUN_PREC_NONE, 0, sunctx);
            MRIStepSetLinearSolver(arkode_mem, LS, nullptr);
        }

        // Set post step function (only on slow integrator)
        MRIStepSetPostprocessStepFn(arkode_mem, SundialsUserFun::post_step);
    }

public:
    SundialsIntegrator () {}

    SundialsIntegrator (const T& S_data, const amrex::Real time = 0.0)
    {
        initialize(S_data, time);
    }

    void initialize (const T& S_data, const amrex::Real time = 0.0) override
    {
        initialize_parameters();
        MPI_Comm mpi_comm = ParallelContext::CommunicatorSub();
#if defined(SUNDIALS_VERSION_MAJOR) && (SUNDIALS_VERSION_MAJOR < 7)
        SUNContext_Create(&mpi_comm, &sunctx);
#else
#  ifdef AMREX_USE_MPI
        SUNContext_Create(mpi_comm, &sunctx);
#  else
        SUNContext_Create(SUN_COMM_NULL, &sunctx);
#  endif
#endif

        // Utility to unpack SUNDIALS vectors
        auto unpack_vector = [&](N_Vector y_data, amrex::Vector<amrex::MultiFab>& S_data) -> void {

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);

            for(int i = 0; i < num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),
                                             amrex::make_alias,
                                             0,
                                             amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
            }
        };

        // Right-hand side function wrappers
        udata.f = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {

            amrex::Vector<amrex::MultiFab> S_data;
            unpack_vector(y_data, S_data);

            amrex::Vector<amrex::MultiFab> S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::pre_rhs_action(S_data, rhs_time);
            BaseT::Rhs(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.fi = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {

            amrex::Vector<amrex::MultiFab> S_data;
            unpack_vector(y_data, S_data);

            amrex::Vector<amrex::MultiFab> S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::pre_rhs_action(S_data, rhs_time);
            BaseT::RhsIm(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.fe = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {

            amrex::Vector<amrex::MultiFab> S_data;
            unpack_vector(y_data, S_data);

            amrex::Vector<amrex::MultiFab> S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::pre_rhs_action(S_data, rhs_time);
            BaseT::RhsEx(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.ff = [&](amrex::Real rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {

            amrex::Vector<amrex::MultiFab> S_data;
            unpack_vector(y_data, S_data);

            amrex::Vector<amrex::MultiFab> S_rhs;
            unpack_vector(y_rhs, S_rhs);

            BaseT::pre_rhs_action(S_data, rhs_time);
            BaseT::RhsFast(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.post_step = [&](amrex::Real time, N_Vector y_data, void * /* user_data */) -> int {

            amrex::Vector<amrex::MultiFab> S_data;
            unpack_vector(y_data, S_data);

            BaseT::post_step_action(S_data, time);

            return 0;
        };

        // Create an N_Vector wrapper for the solution MultiFab
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_data[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        NVar        = S_data.size();      // NOTE: expects S_data to be a Vector<MultiFab>
        nv_many_arr = new N_Vector[NVar]; // vector array composed of cons, xmom, ymom, zmom component vectors */

        for (int i = 0; i < NVar; ++i) {
            nv_many_arr[i] = amrex::sundials::N_VNew_MultiFab(get_length(i),
                                                              S_data[i].boxArray(),
                                                              S_data[i].DistributionMap(),
                                                              S_data[i].nComp(),
                                                              S_data[i].nGrow());

            MultiFab::Copy(*amrex::sundials::getMFptr(nv_many_arr[i]),
                           S_data[i],
                           0,
                           0,
                           S_data[i].nComp(),
                           S_data[i].nGrow());
        }
        nv_S = N_VNew_ManyVector(NVar, nv_many_arr, sunctx);

        if (use_ark) {
            SetupRK(time, nv_S);
        }
        else if (use_mri)
        {
            SetupMRI(time, nv_S);
        }
    }

    virtual ~SundialsIntegrator () {
        // Clean up allocated memory
        for (int i = 0; i < NVar; ++i) {
            N_VDestroy(nv_many_arr[i]);
        }
        delete[] nv_many_arr;
        N_VDestroy(nv_S);
        SUNLinSolFree(LS);
        SUNLinSolFree(fast_LS);
        MRIStepInnerStepper_Free(&fast_stepper);
        MRIStepFree(&arkode_fast_mem);
        ARKStepFree(&arkode_mem);
        SUNContext_Free(&sunctx);
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step) override
    {
        amrex::Real tout = time + time_step;
        amrex::Real tret;

        // Copy the S_old to nv_S
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), S_old[i], 0, 0, S_old[i].nComp(), S_old[i].nGrow());
        }

        if (use_ark) {
            ARKStepReset(arkode_mem, time, nv_S);
            ARKStepSetFixedStep(arkode_mem, time_step);
            int flag = ARKStepEvolve(arkode_mem, tout, nv_S, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
            MRIStepReset(arkode_mem, time, nv_S);
            MRIStepSetFixedStep(arkode_mem, time_step);
            int flag = MRIStepEvolve(arkode_mem, tout, nv_S, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        } else {
            Error("SUNDIALS integrator type not specified.");
        }

        // Copy the result stored in nv_S to state_new
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(S_new[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, S_new[i].nComp(), S_new[i].nGrow());
        }

        return time_step;
    }

    void evolve (T& S_out, const amrex::Real t_out) override
    {
        amrex::Real tret;

        if (use_ark) {
            int flag = ARKStepEvolve(arkode_mem, t_out, nv_S, &tret, ARK_NORMAL);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
            int flag = MRIStepEvolve(arkode_mem, t_out, nv_S, &tret, ARK_NORMAL);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        } else {
            Error("SUNDIALS integrator type not specified.");
        }

        // Should be able to make an alias to S_out and avoid the copy
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(S_out[i],
                           *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)),
                           0,
                           0,
                           S_out[i].nComp(),
                           S_out[i].nGrow());
        }
    }

    void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override {}

    void map_data (std::function<void(T&)> /* Map */) override {}
};

}

#endif
