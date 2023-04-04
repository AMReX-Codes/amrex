#ifndef AMREX_SUNDIALS_INTEGRATOR_H
#define AMREX_SUNDIALS_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <functional>
#include <arkode/arkode_erkstep.h>     /* prototypes for ERKStep fcts., consts   */
#include <arkode/arkode_arkstep.h>     /* prototypes for ARKStep fcts., consts   */
#include <arkode/arkode_mristep.h>     /* prototypes for MRIStep fcts., consts   */
#include <nvector/nvector_manyvector.h>/* manyvector N_Vector types, fcts. etc   */
#include <AMReX_NVector_MultiFab.H>    /* MultiFab N_Vector types, fcts., macros */
#include <AMReX_Sundials.H>            /* MultiFab N_Vector types, fcts., macros */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver        */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to FixedPoint SUNNonlinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype, etc   */

namespace amrex {

struct SundialsUserData {
    std::function<int(realtype, N_Vector, N_Vector, void*)> f0;
    std::function<int(realtype, N_Vector, N_Vector, void*)> f_fast;
    std::function<int(realtype, N_Vector, N_Vector, void*)> f;
    /* std::function<int(realtype, N_Vector*, int, void*)> StoreStage; */
    std::function<int(realtype, N_Vector, void*)> ProcessStage;
    std::function<int(realtype, N_Vector, void*)> PostStoreStage;
};

namespace SundialsUserFun {
    static int f0 (realtype t, N_Vector y, N_Vector ydot, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->f0(t, y, ydot, user_data);
    }

    static int f_fast (realtype t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->f_fast(t, y_data, y_rhs, user_data);
    }

    static int f (realtype t, N_Vector y_data, N_Vector y_rhs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->f(t, y_data, y_rhs, user_data);
    }

/*
    static int StoreStage (realtype t, N_Vector* f_data, int nvecs, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->StoreStage(t, f_data, nvecs, user_data);
    }
*/

    static int ProcessStage (realtype t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->ProcessStage(t, y_data, user_data);
    }

    static int PostStoreStage(realtype t, N_Vector y_data, void *user_data) {
        SundialsUserData* udata = static_cast<SundialsUserData*>(user_data);
        return udata->PostStoreStage(t, y_data, user_data);
    }
}

template<class T>
class SundialsIntegrator : public IntegratorBase<T>
{
private:
    amrex::Real timestep;
    typedef IntegratorBase<T> BaseT;

    bool use_erk_strategy;
    bool use_mri_strategy;
    bool use_mri_strategy_test;
    bool use_implicit_inner;

    SUNNonlinearSolver NLS;    /* empty nonlinear solver object */
    SUNLinearSolver LS;        /* empty linear solver object    */
    void *arkode_mem;          /* empty ARKode memory structure */
    SUNNonlinearSolver NLSf;   /* empty nonlinear solver object */
    SUNLinearSolver LSf;       /* empty linear solver object    */
    void *inner_mem;           /* empty ARKode memory structure */
    void *mristep_mem;         /* empty ARKode memory structure */
    MPI_Comm mpi_comm;         /* the MPI communicator          */
    SUNContext sunctx;         /* SUNDIALS Context object       */

    std::string mri_outer_method, mri_inner_method, erk_method;

    Real reltol;
    Real abstol;
    Real t;
    Real tout;
    Real hfixed;
    Real hfixed_mri;

    int NVar; // NOTE: expects S_data to be a Vector<MultiFab>
    N_Vector* nv_many_arr;     /* vector array composed of cons, xmom, ymom, zmom component vectors */
    N_Vector nv_S;
    N_Vector nv_stage_data;

    void initialize_parameters ()
    {
        use_erk_strategy=false;
        use_mri_strategy=false;
        use_mri_strategy_test=false;

        amrex::ParmParse pp("integration.sundials");

        std::string theStrategy;

        pp.get("strategy", theStrategy);

        if (theStrategy == "ERK")
        {
            use_erk_strategy=true;
            erk_method = "SSPRK3";
            amrex::ParmParse pp_erk("integration.sundials.erk");
            pp_erk.query("method", erk_method);
        }
        else if (theStrategy == "MRI")
        {
            use_mri_strategy=true;
        }
        else if (theStrategy == "MRITEST")
        {
            use_mri_strategy=true;
            use_mri_strategy_test=true;
        }
        else
        {
            std::string msg("Unknown strategy: ");
            msg += theStrategy;
            amrex::Error(msg.c_str());
        }

        if (theStrategy == "MRI" || theStrategy == "MRITEST")
        {
            use_implicit_inner = false;
            mri_outer_method = "KnothWolke3";
            mri_inner_method = "ForwardEuler";
            amrex::ParmParse pp_mri("integration.sundials.mri");
            pp_mri.query("implicit_inner", use_implicit_inner);
            pp_mri.query("outer_method", mri_outer_method);
            pp_mri.query("inner_method", mri_inner_method);
        }

        // SUNDIALS specific objects
        NLS = nullptr;            /* empty nonlinear solver object */
        LS = nullptr;             /* empty linear solver object    */
        arkode_mem = nullptr;     /* empty ARKode memory structure */
        NLSf = nullptr;           /* empty nonlinear solver object */
        LSf = nullptr;            /* empty linear solver object    */
        inner_mem = nullptr;      /* empty ARKode memory structure */
        mristep_mem = nullptr;    /* empty ARKode memory structure */

        // Arbitrary tolerances
        reltol          = 1e-4;
        abstol          = 1e-4;
    }

public:
    SundialsIntegrator () {}

    SundialsIntegrator (const T& /* S_data */)
    {
        initialize();
    }

    void initialize (const T& /* S_data */) override
    {
        initialize_parameters();
        mpi_comm = ParallelContext::CommunicatorSub();
        SUNContext_Create(&mpi_comm, &sunctx);
    }

    void initialize ()
    {
        initialize_parameters();
        mpi_comm = ParallelContext::CommunicatorSub();
        SUNContext_Create(&mpi_comm, &sunctx);
    }

    virtual ~SundialsIntegrator () {
        SUNContext_Free(&sunctx);
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step) override
    {
        if (use_mri_strategy) {
            return advance_mri(S_old, S_new, time, time_step);
        } else if (use_erk_strategy) {
            return advance_erk(S_old, S_new, time, time_step);
        } else {
            Error("SUNDIALS integrator backend not specified (ERK or MRI).");
        }

        return 0;
    }

    amrex::Real advance_erk (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        t               = time;
        tout            = time+time_step;
        hfixed          = time_step;
        timestep = time_step;

        // We use S_new as our working space, so first copy S_old to S_new
        IntegratorOps<T>::Copy(S_new, S_old);

        // Create an N_Vector wrapper for the solution MultiFab
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_new[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        /* Create manyvector for solution using S_new */
        NVar             = S_new.size();  // NOTE: expects S_new to be a Vector<MultiFab>
        nv_many_arr = new N_Vector[NVar]; // vector array composed of cons, xmom, ymom, zmom component vectors */

        for (int i = 0; i < NVar; ++i) {
            sunindextype length = get_length(i);
            N_Vector nvi   = amrex::sundials::N_VMake_MultiFab(length, &S_new[i]);
            nv_many_arr[i] = nvi;
        }

        nv_S = N_VNew_ManyVector(NVar, nv_many_arr, sunctx);
        nv_stage_data = N_VClone(nv_S);

        /* Create a temporary storage space for MRI */
        Vector<std::unique_ptr<T> > temp_storage;
        IntegratorOps<T>::CreateLike(temp_storage, S_old);
        T& state_store = *temp_storage.back();

        SundialsUserData udata;

        /* Begin Section: SUNDIALS FUNCTION HOOKS */
        /* f routine to compute the ODE RHS function f(t,y). */
        udata.f = [&](realtype rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {
            amrex::Vector<amrex::MultiFab> S_data;
            amrex::Vector<amrex::MultiFab> S_rhs;

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);
            S_rhs.resize(num_vecs);

            for(int i=0; i<num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
                S_rhs.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i))->nComp());
            }

            BaseT::post_update(S_data, rhs_time);
            BaseT::rhs(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.ProcessStage = [&](realtype rhs_time, N_Vector y_data, void * /* user_data */) -> int {
            amrex::Vector<amrex::MultiFab > S_data;

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);

            for (int i=0; i<num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
            }

            BaseT::post_update(S_data, rhs_time);

            return 0;
        };
        /* End Section: SUNDIALS FUNCTION HOOKS */

        /* Call ERKStepCreate to initialize the inner ARK timestepper module and
        specify the right-hand side function in y'=f(t,y), the initial time
        T0, and the initial dependent variable vector y. */
        arkode_mem = ERKStepCreate(SundialsUserFun::f, time, nv_S, sunctx);
        ERKStepSetUserData(arkode_mem, &udata);  /* Pass udata to user functions */
        ERKStepSetPostprocessStageFn(arkode_mem, SundialsUserFun::ProcessStage);
        /* Specify tolerances */
        ERKStepSStolerances(arkode_mem, reltol, abstol);
        ERKStepSetFixedStep(arkode_mem, hfixed);

        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(state_store[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_store[i].nComp(), state_store[i].nGrow());
            MultiFab::Copy(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_stage_data, i)), *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_store[i].nComp(), state_store[i].nGrow());
        }

        auto make_butcher_table = [&](std::string method) -> ARKodeButcherTable {
            ARKodeButcherTable B;
            if (method == "SSPRK3") {
                B = ARKodeButcherTable_Alloc(3, SUNFALSE);

                // 3rd order Strong Stability Preserving RK3
                B->A[1][0] = 1.0;
                B->A[2][0] = 0.25;
                B->A[2][1] = 0.25;
                B->b[0] = 1./6.;
                B->b[1] = 1./6.;
                B->b[2] = 2./3.;
                B->c[1] = 1.0;
                B->c[2] = 0.5;
                B->q=3;
            } else if (method == "Trapezoid") {
                B = ARKodeButcherTable_Alloc(2, SUNFALSE);

                // Trapezoidal rule
                B->A[1][0] = 1.0;
                B->b[0] = 0.5;
                B->b[1] = 0.5;
                B->c[1] = 1.0;
                B->q=2;
                B->p=0;
            } else if (method == "ForwardEuler") {
                B = ARKodeButcherTable_Alloc(1, SUNFALSE);

                // Forward Euler
                B->b[0] = 1.0;
                B->q=1;
                B->p=0;
            } else
                amrex::Error("ERK method not implemented");
            return B;
        };

        ARKodeButcherTable B = make_butcher_table(erk_method);

        //Set table
        ERKStepSetTable(arkode_mem, B);

        // Free the Butcher table
        ARKodeButcherTable_Free(B);

        // Use ERKStep to evolve state_old data (wrapped in nv_S) from t to tout=t+dt
        auto flag = ERKStepEvolve(arkode_mem, tout, nv_S, &t, ARK_NORMAL);
        AMREX_ALWAYS_ASSERT(flag >= 0);

        // Copy the result stored in nv_S to state_new
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(S_new[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, S_new[i].nComp(), S_new[i].nGrow());
        }

        // Clean up allocated memory
        N_VDestroyVectorArray(nv_many_arr, NVar);
        N_VDestroy(nv_S);
        N_VDestroy(nv_stage_data);

        ERKStepFree(&arkode_mem);

        // Return timestep
        return timestep;
    }

    amrex::Real advance_mri (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        int mri_time_step_ratio = BaseT::get_slow_fast_timestep_ratio();
        Real mri_fast_time_step = BaseT::get_fast_timestep();
        AMREX_ALWAYS_ASSERT(mri_time_step_ratio >= 1 || mri_fast_time_step >= 0.0);
        t               = time;
        tout            = time+time_step;
        hfixed          = time_step;
        hfixed_mri      = mri_fast_time_step >= 0.0 ? mri_fast_time_step : time_step / mri_time_step_ratio;
        timestep = time_step;

        // NOTE: hardcoded for now ...
        bool use_erk3 = !use_implicit_inner;
        bool use_linear = false;

        // We use S_new as our working space, so first copy S_old to S_new
        IntegratorOps<T>::Copy(S_new, S_old);

        // Create an N_Vector wrapper for the solution MultiFab
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_new[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        /* Create manyvector for solution using S_new */
        NVar             = S_new.size();  // NOTE: expects S_new to be a Vector<MultiFab>
        nv_many_arr = new N_Vector[NVar]; // vector array composed of cons, xmom, ymom, zmom component vectors */

        for (int i = 0; i < NVar; ++i) {
            sunindextype length = get_length(i);
            N_Vector nvi   = amrex::sundials::N_VMake_MultiFab(length, &S_new[i]);
            nv_many_arr[i] = nvi;
        }

        nv_S = N_VNew_ManyVector(NVar, nv_many_arr, sunctx);
        nv_stage_data = N_VClone(nv_S);

        // Copy the initial step data to nv_stage_data
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
           MultiFab* mf_y = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i));
           MultiFab* mf_stage = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_stage_data, i));
           MultiFab::Copy(*mf_stage, *mf_y, 0, 0, mf_y->nComp(), mf_y->nGrow());
        }

        /* Create a temporary storage space for MRI */
        Vector<std::unique_ptr<T> > temp_storage;
        IntegratorOps<T>::CreateLike(temp_storage, S_old);
        T& state_store = *temp_storage.back();

        SundialsUserData udata;

        /* Begin Section: SUNDIALS FUNCTION HOOKS */
        /* f0 routine to compute a zero-valued ODE RHS function f(t,y). */
        udata.f0 = [&](realtype /* rhs_time */, N_Vector /* y */, N_Vector ydot, void * /* user_data */) -> int {
            // Initialize ydot to zero and return
            N_VConst(0.0, ydot);
            return 0;
        };

        /* f routine to compute the ODE RHS function f(t,y). */
        udata.f_fast = [&](realtype rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {
            amrex::Vector<amrex::MultiFab> S_data;
            amrex::Vector<amrex::MultiFab> S_rhs;
            amrex::Vector<amrex::MultiFab> S_stage_data;

            N_VConst(0.0, y_rhs);

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);
            S_rhs.resize(num_vecs);
            S_stage_data.resize(num_vecs);

            for(int i=0; i<num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
                S_rhs.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i))->nComp());
                S_stage_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_stage_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_stage_data, i))->nComp());
            }

            // NOTE: we can optimize by calling a post_update_fast and only updating the variables the fast integration modifies
            BaseT::post_update(S_data, rhs_time);

            BaseT::fast_rhs(S_rhs, S_stage_data, S_data, rhs_time);

            return 0;
        };

        /* f routine to compute the ODE RHS function f(t,y). */
        udata.f = [&](realtype rhs_time, N_Vector y_data, N_Vector y_rhs, void * /* user_data */) -> int {
            amrex::Vector<amrex::MultiFab> S_data;
            amrex::Vector<amrex::MultiFab> S_rhs;

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);
            S_rhs.resize(num_vecs);

            for(int i=0; i<num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
                S_rhs.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i))->nComp());
            }

            BaseT::post_update(S_data, rhs_time);
            BaseT::rhs(S_rhs, S_data, rhs_time);

            return 0;
        };

        udata.ProcessStage = [&](realtype rhs_time, N_Vector y_data, void * /* user_data */) -> int {
            amrex::Vector<amrex::MultiFab > S_data;

            const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
            S_data.resize(num_vecs);

            for (int i=0; i<num_vecs; i++)
            {
                S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
            }

            BaseT::post_update(S_data, rhs_time);

            return 0;
        };

        udata.PostStoreStage = [&](realtype rhs_time, N_Vector y_data, void *user_data) -> int {
            udata.ProcessStage(rhs_time, y_data, user_data);

            for(int i=0; i<N_VGetNumSubvectors_ManyVector(y_data); i++)
            {
                MultiFab* mf_y = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i));
                MultiFab* mf_stage = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_stage_data, i));
                MultiFab::Copy(*mf_stage, *mf_y, 0, 0, mf_y->nComp(), mf_y->nGrow());
            }

            return 0;
        };
        /* End Section: SUNDIALS FUNCTION HOOKS */

        if(use_mri_strategy_test)
        {
            if(use_erk3)
                inner_mem = ARKStepCreate(SundialsUserFun::f0, nullptr, time, nv_S, sunctx); // explicit bc (explicit f, implicit f, time, data)
            else
                inner_mem = ARKStepCreate(nullptr, SundialsUserFun::f0, time, nv_S, sunctx); // implicit
        }
        else
        {
            if(use_erk3)
                inner_mem = ARKStepCreate(SundialsUserFun::f_fast, nullptr, time, nv_S, sunctx);
            else
                inner_mem = ARKStepCreate(nullptr, SundialsUserFun::f_fast, time, nv_S, sunctx);
        }

        ARKStepSetFixedStep(inner_mem, hfixed_mri);            // Specify fixed time step size

        ARKStepSetUserData(inner_mem, &udata);  /* Pass udata to user functions */

        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(state_store[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_store[i].nComp(), state_store[i].nGrow());
            MultiFab::Copy(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_stage_data, i)), *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_store[i].nComp(), state_store[i].nGrow());
        }

        auto make_butcher_table = [&](std::string method) -> ARKodeButcherTable {
            ARKodeButcherTable B;
            if (method == "KnothWolke3" || method == "Knoth-Wolke-3-3") {
                B = ARKodeButcherTable_Alloc(3, SUNFALSE);

                // 3rd order Knoth-Wolke method
                B->A[1][0] = 1.0/3.0;
                B->A[2][0] = -3.0/16.0;
                B->A[2][1] = 15.0/16.0;
                B->b[0] = 1./6.;
                B->b[1] = 3./10.;
                B->b[2] = 8./15.;
                B->c[1] = 1.0/3.0;
                B->c[2] = 3.0/4.0;
                B->q=3;
                B->p=0;
            } else if (method == "Trapezoid") {
                B = ARKodeButcherTable_Alloc(2, SUNFALSE);

                // Trapezoidal rule
                B->A[1][0] = 1.0;
                B->b[0] = 0.5;
                B->b[1] = 0.5;
                B->c[1] = 1.0;
                B->q=2;
                B->p=0;
            } else if (method == "ForwardEuler") {
                B = ARKodeButcherTable_Alloc(1, SUNFALSE);

                // Forward Euler
                B->b[0] = 1.0;
                B->q=1;
                B->p=0;
            } else
                amrex::Error("MRI method not implemented");
            return B;
        };

        ARKodeButcherTable B_outer = make_butcher_table(mri_outer_method);
        ARKodeButcherTable B_inner = make_butcher_table(mri_inner_method);

        if(use_erk3)
        {
            ARKStepSetTables(inner_mem, B_inner->q, B_inner->p, nullptr, B_inner);       // Specify Butcher table
        } else {
            ARKodeButcherTable_Free(B_inner);
            B_inner = ARKodeButcherTable_Alloc(2, SUNFALSE);

            B_inner->A[1][0] = 1.0;
            B_inner->A[2][0] = 1.0;
            B_inner->A[2][2] = 0.0;
            B_inner->b[0] = 0.5;
            B_inner->b[2] = 0.5;
            B_inner->c[1] = 1.0;
            B_inner->c[2] = 1.0;
            B_inner->q=2;
            ARKStepSetTables(inner_mem, B_inner->q, B_inner->p, B_inner, nullptr);       // Specify Butcher table
        }

        //Set table
        // Create fast time scale integrator from an ARKStep instance
        MRIStepInnerStepper inner_stepper = nullptr;
        ARKStepCreateMRIStepInnerStepper(inner_mem, &inner_stepper);

        // args: fast RHS, nullptr, initial time, initial state, fast time scale integrator, sundials context
        mristep_mem = MRIStepCreate(SundialsUserFun::f, nullptr, time, nv_S, inner_stepper, sunctx);

        MRIStepSetFixedStep(mristep_mem, hfixed);

        /* Specify tolerances */
        MRIStepSStolerances(mristep_mem, reltol, abstol);

        /* Initialize spgmr solver */
        LS = SUNLinSol_SPGMR(nv_S, PREC_NONE, 10, sunctx);
        NLS = SUNNonlinSol_FixedPoint(nv_S, 50, sunctx);

        if (use_implicit_inner) ARKStepSetNonlinearSolver(inner_mem, NLS);
        if(use_linear)
            MRIStepSetLinearSolver(mristep_mem, LS, nullptr);
        else
            MRIStepSetNonlinearSolver(mristep_mem, NLS);

        MRIStepSetUserData(mristep_mem, &udata);  /* Pass udata to user functions */
        MRIStepSetPostprocessStageFn(mristep_mem, SundialsUserFun::ProcessStage);

        MRIStepCoupling mri_coupling = MRIStepCoupling_MIStoMRI(B_outer, B_outer->q, B_outer->p);
        MRIStepSetCoupling(mristep_mem, mri_coupling);

        // Free the Butcher tables
        ARKodeButcherTable_Free(B_outer);
        ARKodeButcherTable_Free(B_inner);

        // Use MRIStep to evolve state_old data (wrapped in nv_S) from t to tout=t+dt
        auto flag = MRIStepEvolve(mristep_mem, tout, nv_S, &t, ARK_NORMAL);
        AMREX_ALWAYS_ASSERT(flag >= 0);

        // Copy the result stored in nv_S to state_new
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(S_new[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, S_new[i].nComp(), S_new[i].nGrow());
        }

        // Clean up allocated memory
        N_VDestroyVectorArray(nv_many_arr, NVar);
        N_VDestroy(nv_S);
        N_VDestroy(nv_stage_data);

        MRIStepFree(&mristep_mem);
        ARKStepFree(&inner_mem);

        SUNLinSolFree(LS);
        SUNNonlinSolFree(NLS);

        // Return timestep
        return timestep;
    }

    void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override {}

    void map_data (std::function<void(T&)> /* Map */) override {}

};

}

#endif
