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

namespace amrex {

struct SundialsUserData {
    std::function<int(sunrealtype, N_Vector, N_Vector, void*)> f;
    std::function<int(sunrealtype, N_Vector, N_Vector, void*)> fe;
    std::function<int(sunrealtype, N_Vector, N_Vector, void*)> fi;
};

namespace SundialsUserFun {
    static int f (sunrealtype t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->f(t, y_data, y_rhs, user_data);
    }

    static int fe (sunrealtype t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->fe(t, y_data, y_rhs, user_data);
    }

    static int fi (sunrealtype t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->fi(t, y_data, y_rhs, user_data);
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

    // method type flags
    bool use_erk  = false;
    bool use_dirk = false;
    bool use_imex = false;
    bool use_mri  = false;

    // structure for interfacing with user-supplied functions
    SundialsUserData udata;

    // SUNDIALS objects
    SUNContext sunctx  = nullptr; // should get from singleton
    void *arkode_mem   = nullptr;
    SUNLinearSolver LS = nullptr;

    // relative and absolute tolerances
    Real reltol = 1.0e-4;
    Real abstol = 1.0e-9;

    int NVar;               // NOTE: expects S_data to be a Vector<MultiFab>
    N_Vector* nv_many_arr;  /* vector array composed of cons, xmom, ymom, zmom component vectors */
    N_Vector nv_S;

    void initialize_parameters ()
    {
        amrex::ParmParse pp("integration.sundials");

        pp.query("type", type);
        pp.query("method", method);

        if (type == "ERK") {
            use_erk = true;
        }
        else if (type == "DIRK") {
            use_dirk = true;
        }
        else if (type == "IMEX") {
            use_imex = true;
        }
        else {
            std::string msg("Unknown strategy: ");
            msg += type;
            amrex::Error(msg.c_str());
        }
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

        /* f routine to compute the ODE RHS function f(t,y). */
        udata.f = [&](sunrealtype rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {
            amrex::Vector<amrex::MultiFab> S_data;
            amrex::Vector<amrex::MultiFab> S_rhs;

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);
            S_rhs.resize(num_vecs);

            for(int i=0; i<num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),
                                             amrex::make_alias,
                                             0,
                                             amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());

                S_rhs.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i)),
                                            amrex::make_alias,
                                            0,
                                            amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i))->nComp());
            }

            BaseT::post_update(S_data, rhs_time);
            BaseT::rhs(S_rhs, S_data, rhs_time);

            return 0;
        };



        // Create an N_Vector wrapper for the solution MultiFab
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_data[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        NVar        = S_data.size();       // NOTE: expects S_data to be a Vector<MultiFab>
        nv_many_arr = new N_Vector[NVar]; // vector array composed of cons, xmom, ymom, zmom component vectors */

        for (int i = 0; i < NVar; ++i) {
            nv_many_arr[i] = amrex::sundials::N_VNew_MultiFab(get_length(i),
                                                              S_data[i].boxArray(),
                                                              S_data[i].DistributionMap(),
                                                              S_data[i].nComp(),
                                                              S_data[i].nGrow());
            MultiFab::Copy(*amrex::sundials::getMFptr(nv_many_arr[i]), S_data[i], 0, 0, S_data[i].nComp(), S_data[i].nGrow());
        }
        nv_S = N_VNew_ManyVector(NVar, nv_many_arr, sunctx);

        if (use_erk) {
            amrex::Print() << "SUNDIALS ERK time integrator\n";
            arkode_mem = ARKStepCreate(SundialsUserFun::f, nullptr, time, nv_S, sunctx);
        }
        else if (use_dirk) {
            amrex::Print() << "SUNDIALS DIRK time integrator\n";
            arkode_mem = ARKStepCreate(nullptr, SundialsUserFun::f, time, nv_S, sunctx);
        }
        else if (use_imex) {
            amrex::Print() << "SUNDIALS IMEX time integrator\n";
            arkode_mem = ARKStepCreate(SundialsUserFun::fe, SundialsUserFun::fi, time, nv_S, sunctx);
        }
        ARKStepSetUserData(arkode_mem, &udata);  /* Pass udata to user functions */
        ARKStepSStolerances(arkode_mem, reltol, abstol);

        if (method != "DEFAULT" || (method_e != "DEFAULT" && method_i != "DEFAULT"))
        {
            if (use_erk) {
                amrex::Print() << "SUNDIALS ERK method " << method << "\n";
                ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", method.c_str());
            }
            else if (use_dirk) {
                amrex::Print() << "SUNDIALS DIRK method " << method << "\n";
                ARKStepSetTableName(arkode_mem, method.c_str(), "ARKODE_ERK_NONE");
            }
            else {
                amrex::Print() << "SUNDIALS IMEX method " << method_i << " and "
                               << method_e << "\n";
                ARKStepSetTableName(arkode_mem, method_i.c_str(), method_e.c_str());
            }
        }

        // Create and attach linear solver
        if (use_dirk || use_imex) {
            LS = SUNLinSol_SPGMR(nv_S, SUN_PREC_NONE, 0, sunctx);
            ARKStepSetLinearSolver(arkode_mem, LS, nullptr);
        }
    }

    virtual ~SundialsIntegrator () {
        // Clean up allocated memory
        for (int i = 0; i < NVar; ++i) {
            N_VDestroy(nv_many_arr[i]);
        }
        delete[] nv_many_arr;
        N_VDestroy(nv_S);
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

        if (use_erk || use_dirk || use_imex) {
            ARKStepReset(arkode_mem, time, nv_S);
            ARKStepSetFixedStep(arkode_mem, time_step);
            int flag = ARKStepEvolve(arkode_mem, tout, nv_S, &tret, ARK_ONE_STEP);
            AMREX_ALWAYS_ASSERT(flag >= 0);
        }
        else if (use_mri) {
            Error("SUNDIALS integrator type not implemented, yet.");
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

    void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override {}

    void map_data (std::function<void(T&)> /* Map */) override {}
};

}

#endif
