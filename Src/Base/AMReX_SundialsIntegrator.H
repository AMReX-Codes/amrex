#ifndef AMREX_SUNDIALS_INTEGRATOR_H
#define AMREX_SUNDIALS_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntegratorBase.H>
#include <functional>
#include <arkode/arkode_erkstep.h>     /* prototypes for ERKStep fcts., consts */
#include <arkode/arkode_arkstep.h>     /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_mristep.h>     /* prototypes for MRIStep fcts., consts */
#include <nvector/nvector_manyvector.h>/* manyvector N_Vector types, fcts. etc */
#include <AMReX_NVector_MultiFab.H>    /* MultiFab N_Vector types, fcts., macros */
#include <AMReX_Sundials.H>    /* MultiFab N_Vector types, fcts., macros */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to FixedPoint SUNNonlinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype, etc */

namespace amrex {

//! The stepper strategies
enum Strategy { UNDEFINED = -1, NATIVE, ERK, MRI, MRITEST };

struct FastRhsData {
  std::function<void(amrex::Vector<amrex::MultiFab> &,
                     const amrex::Vector<amrex::MultiFab> &,
                     const amrex::Vector<amrex::MultiFab> &,
                     const Real)>  rhs_fun_fast;
  std::function<void(amrex::Vector<amrex::MultiFab> &,
                     const amrex::Vector<amrex::MultiFab> &,
                     const Real)>  rhs_fun;
  std::function<void(amrex::Vector<amrex::MultiFab> &,
                     const Real)>  post_update;
  void* inner_mem;
  amrex::Vector<amrex::MultiFab>* S_stage_data; // hold previous slow stage data
  N_Vector nv_stage_data; // hold previous slow stage data
};

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_fast(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int StoreStage(realtype t, N_Vector* f_data, int nvecs, void *user_data);
static int PostStoreStage(realtype t, N_Vector y_data, void *user_data);
static int ProcessStage(realtype t, N_Vector y_data, void *user_data);

typedef IntegratorBase<Vector<MultiFab> > IBaseVM;

template<class T>
class SundialsIntegrator : public IntegratorBase<T>
{
private:
    amrex::Real timestep;
    typedef IntegratorBase<T> BaseT;

    bool advance_erk;
    bool advance_mri;
    bool advance_mri_test;

    SUNNonlinearSolver NLS;    /* empty nonlinear solver object */
    SUNLinearSolver LS;    /* empty linear solver object */
    void *arkode_mem;      /* empty ARKode memory structure */
    SUNNonlinearSolver NLSf;    /* empty nonlinear solver object */
    SUNLinearSolver LSf;    /* empty linear solver object */
    void *inner_mem;      /* empty ARKode memory structure */
    void *mristep_mem;      /* empty ARKode memory structure */

    Real reltol;
    Real abstol;
    Real t;
    Real tout;
    Real hfixed;
    Real m;
    Real hfixed_mri;

    int NVar; // NOTE: expects S_data to be a Vector<MultiFab>
    N_Vector* nv_many_arr;              /* vector array composed of cons, xmom, ymom, zmom component vectors */
    N_Vector nv_S;
    N_Vector nv_store;

    void initialize_parameters ()
    {
        advance_erk=false;
        advance_mri=false;
        advance_mri_test=false;

        amrex::ParmParse pp("integration.sundials");

        std::string theStrategy;

        pp.get("strategy", theStrategy);

        if (theStrategy == "ERK")
        {
            advance_erk=true;
        }
        else if (theStrategy == "MRI")
        {
            advance_mri=true;
        }
        else if (theStrategy == "MRITEST")
        {
            advance_mri=true;
            advance_mri_test=true;
        }
        else
        {
            std::string msg("Unknown strategy: ");
            msg += theStrategy;
            amrex::Error(msg.c_str());
        }
    }

    void initialize_stages (const T& S_data)
    {
        ////STEP ONE
        // Create SUNDIALS specific objects
        NLS = NULL;    /* empty nonlinear solver object */
        LS = NULL;    /* empty linear solver object */
        arkode_mem = NULL;      /* empty ARKode memory structure */
        NLSf = NULL;    /* empty nonlinear solver object */
        LSf = NULL;    /* empty linear solver object */
        inner_mem = NULL;      /* empty ARKode memory structure */
        mristep_mem = NULL;      /* empty ARKode memory structure */

        // Create an N_Vector wrapper for the solution MultiFab
        auto get_length = [&](int index) -> sunindextype {
            auto* p_mf = &S_data[index];
            return p_mf->nComp() * (p_mf->boxArray()).numPts();
        };

        ////STEP TWO

        //Arbitrary tolerances
        reltol          = 1e-4;
        abstol          = 1e-4;
        m               = 2;

        ////STEP THREE
        /* Create manyvector for solution */
        NVar             = S_data.size(); // NOTE: expects S_data to be a Vector<MultiFab>
        nv_many_arr = new N_Vector[NVar];              /* vector array composed of cons, xmom, ymom, zmom component vectors */

        for (int i = 0; i < NVar; ++i) {
            sunindextype length = get_length(i);
            N_Vector nvi    = amrex::sundials::N_VMake_MultiFab(length, &S_data[i]);
            nv_many_arr[i] = nvi;
        }

        nv_S = N_VNew_ManyVector(NVar, nv_many_arr);
        nv_store = N_VClone(nv_S);
    }

public:
    SundialsIntegrator () {}

    SundialsIntegrator (const T& S_data)
    {
        initialize(S_data);
    }

    void initialize (const T& S_data)
    {
        initialize_parameters();
        initialize_stages(S_data);
    }

    virtual ~SundialsIntegrator () {
        ////STEP THIRTEEN
        N_VDestroyVectorArray(NVar, nv_many_arr);
        N_VDestroy(nv_S);
        N_VDestroy(nv_store);

        ////STEP FOURTEEN
        if(advance_mri) {
            MRIStepFree(&mristep_mem);
        }

        ARKStepFree(&inner_mem);
        ERKStepFree(&arkode_mem);

        if(advance_mri)
        {
            SUNLinSolFree(LS);    // Free nonlinear solvers
            SUNNonlinSolFree(NLS);    // Free nonlinear solvers
        }
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step)
    {
        t               = time;
        tout            = time+time_step;
        hfixed          = time_step;
        hfixed_mri      = time_step / m;
        timestep = time_step;

        // NOTE: hardcoded for now ...
        bool use_erk3 = true;
        bool use_linear = false;

        Vector<T> temp_storage;
        IntegratorOps<T>::CreateLike(temp_storage, S_old);
        T& state_store = temp_storage.back();

        ////STEP FOUR
        if(advance_mri_test)
        {
            if(use_erk3)
                inner_mem = ARKStepCreate(f0, NULL, time, nv_S); // explicit bc (explicit f, implicit f, time, data)
            else
                inner_mem = ARKStepCreate(NULL, f0, time, nv_S); // implicit
        }
        else
        {
            if(use_erk3)
                inner_mem = ARKStepCreate(f_fast, NULL, time, nv_S);
            else
                inner_mem = ARKStepCreate(NULL, f_fast, time, nv_S);
        }

        FastRhsData fast_userdata;
        fast_userdata.rhs_fun_fast = BaseT::get_fast_rhs();
        //For the sundials solve, use state_new as temporary data;
        fast_userdata.S_stage_data = &state_store;
        fast_userdata.nv_stage_data = nv_store;
        fast_userdata.inner_mem = inner_mem;

        /* BEGIN: SUNDIALS FUNCTION HOOKS */
        // f0 routine to compute a zero-valued ODE RHS function f(t,y).
        static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data)
        {
        // Initialize ydot to zero and return
        N_VConst(0.0, ydot);
        return 0;
        }

        /* f routine to compute the ODE RHS function f(t,y). */
        static int f_fast(realtype t, N_Vector y_data, N_Vector y_rhs, void *user_data)
        {
        FastRhsData* fast_userdata = (FastRhsData*) user_data;
        amrex::Vector<amrex::MultiFab> S_data;
        amrex::Vector<amrex::MultiFab> S_rhs;
        amrex::Vector<amrex::MultiFab> S_stage_data;
        auto call_post_update = IBaseVM::get_post_update();
        auto call_rhs = IBaseVM::get_rhs();

        N_VConst(0.0, y_rhs);

        const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
        S_data.resize(num_vecs);
        S_rhs.resize(num_vecs);
        S_stage_data.resize(num_vecs);

        for(int i=0; i<num_vecs; i++)
        {
            S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
            S_rhs.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i))->nComp());
            S_stage_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(fast_userdata->nv_stage_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(fast_userdata->nv_stage_data, i))->nComp());
        }

        // NOTE: we can optimize by calling a post_update_fast and only updating the variables the fast integration modifies
        call_post_update(S_data, t);
        call_post_update(S_stage_data, t);

        //Call rhs_fun_fast lambda stored in userdata which uses erf_fast_rhs
        fast_userdata->rhs_fun_fast(S_rhs, S_stage_data, S_data, t);

        return 0;
        }

        static int StoreStage(realtype t, N_Vector* f_data, int nvecs, void *user_data)
        {

        //Note that user_data may be different between mristep and inner stepper
        FastRhsData* fast_userdata = (FastRhsData*) user_data;
        void* inner_mem = fast_userdata->inner_mem;

        N_Vector y_data;
        Real tcur;
        ARKStepGetCurrentState(inner_mem, &y_data);
        ARKStepGetCurrentTime(inner_mem, &tcur);

        //This segfaults if this function is called as PreInner
        const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);

        for(int i=0; i<N_VGetNumSubvectors_ManyVector(y_data); i++)
        {
            MultiFab* mf_y = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i));
            MultiFab* mf_stage = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(fast_userdata->nv_stage_data, i));
            MultiFab::Copy(*mf_stage, *mf_y, 0, 0, mf_y->nComp(), mf_y->nGrow());
        }

        return 0;
        }

        static int PostStoreStage(realtype t, N_Vector y_data, void *user_data)
        {

        FastRhsData* fast_userdata = (FastRhsData*) user_data;

        const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);

        ProcessStage(t, y_data, user_data);

        for(int i=0; i<N_VGetNumSubvectors_ManyVector(y_data); i++)
        {
            MultiFab* mf_y = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i));
            MultiFab* mf_stage = amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(fast_userdata->nv_stage_data, i));
            MultiFab::Copy(*mf_stage, *mf_y, 0, 0, mf_y->nComp(), mf_y->nGrow());
        }

        return 0;
        }

        /* f routine to compute the ODE RHS function f(t,y). */
        static int f(realtype t, N_Vector y_data, N_Vector y_rhs, void *user_data)
        {
        FastRhsData* fast_userdata = (FastRhsData*) user_data;
        amrex::Vector<amrex::MultiFab> S_data;
        amrex::Vector<amrex::MultiFab> S_rhs;
        auto call_post_update = IBaseVM::get_post_update();
        auto call_rhs = IBaseVM::get_rhs();

        const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
        S_data.resize(num_vecs);
        S_rhs.resize(num_vecs);

        for(int i=0; i<num_vecs; i++)
        {
            S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
            S_rhs.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_rhs, i))->nComp());
        }

        call_post_update(S_data, t);
        call_rhs(S_rhs, S_data, t);

        return 0;
        }

        static int ProcessStage(realtype t, N_Vector y_data, void *user_data)
        {
        FastRhsData* fast_userdata = (FastRhsData*) user_data;
        amrex::Vector<amrex::MultiFab > S_data;
        auto call_post_update = BaseT::get_post_update();

        const int num_vecs = N_VGetNumSubvectors_ManyVector(y_data);
        S_data.resize(num_vecs);

        for(int i=0; i<num_vecs; i++)
        {
            S_data.at(i)=amrex::MultiFab(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i)),amrex::make_alias,0,amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(y_data, i))->nComp());
        }

        call_post_update(S_data, t);

        return 0;
        }
        /* END: SUNDIALS FUNCTION HOOKS */

        ////STEP FIVE
        ARKStepSetFixedStep(inner_mem, hfixed_mri);            // Specify fixed time step size

        /* Call ERKStepCreate to initialize the inner ARK timestepper module and
        specify the right-hand side function in y'=f(t,y), the initial time
        T0, and the initial dependent variable vector y. */
        arkode_mem = ERKStepCreate(f, time, nv_S); // NOTE: unused if we are doing multirate
        ERKStepSetUserData(arkode_mem, (void *) &fast_userdata);  /* Pass udata to user functions */
        ERKStepSetPostprocessStageFn(arkode_mem, ProcessStage);
        /* Specify tolerances */
        ERKStepSStolerances(arkode_mem, reltol, abstol);
        ERKStepSetFixedStep(arkode_mem, hfixed);

        ARKStepSetUserData(inner_mem, (void *) &fast_userdata);  /* Pass udata to user functions */

        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(state_store[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_store[i].nComp(), state_store[i].nGrow());
            MultiFab::Copy(*amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_store, i)), *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, state_store[i].nComp(), state_store[i].nGrow());
        }

        ARKodeButcherTable B = ARKodeButcherTable_Alloc(3, SUNFALSE);
        ARKodeButcherTable B2 = ARKodeButcherTable_Alloc(2, SUNFALSE);

        if(use_erk3)
        {
            if(advance_erk)
            {
                // Use SSP-RK3
                B->A[1][0] = 1.0;
                B->A[2][0] = 0.25;
                B->A[2][1] = 0.25;
                B->b[0] = 1./6.;
                B->b[1] = 1./6.;
                B->b[2] = 2./3.;
                B->c[1] = 1.0;
                B->c[2] = 0.5;
                B->q=3;
            }

            if(advance_mri)
            {
                B->A[1][0] = 1.0;
                B->A[2][0] = 1.0;
                B->A[2][2] = 0.0;
                B->b[0] = 0.5;
                B->b[2] = 0.5;
                B->c[1] = 1.0;
                B->c[2] = 1.0;
                B->q=2;
                B->p=0;
                B2->A[1][0] = 1.0;
                B2->b[0] = 0.5;
                B2->b[1] = 0.5;
                B2->c[1] = 1.0;
                B2->q=2;
                B2->p=0;
            }

            if(advance_mri) {
                if(advance_mri_test)
                    ARKStepSetTables(inner_mem, B->q, B->p, NULL, B);       // Specify Butcher table
                else
                    ARKStepSetTables(inner_mem, B2->q, B2->p, NULL, B2);       // Specify Butcher table
            }
        } else {
            B->A[1][0] = 1.0;
            B->A[2][0] = 1.0;
            B->A[2][2] = 0.0;
            B->b[0] = 0.5;
            B->b[2] = 0.5;
            B->c[1] = 1.0;
            B->c[2] = 1.0;
            B->q=2;
            ARKStepSetTables(inner_mem, B->q, B->p, B, NULL);       // Specify Butcher table
        }

        //Set table
        //    ERKStepSetTable(arkode_mem, B);
        if(advance_mri)
        {
            ////STEP SIX
            mristep_mem = MRIStepCreate(f, time, nv_S, MRISTEP_ARKSTEP, inner_mem);

            ////STEP SEVEN
            MRIStepSetFixedStep(mristep_mem, hfixed);

            ////STEP 8.1
            /* Specify tolerances */
            MRIStepSStolerances(mristep_mem, reltol, abstol);
            ////STEP 8.2
            /* Initialize spgmr solver */
            LS = SUNLinSol_SPGMR(nv_S, PREC_NONE, 10);
            NLS = SUNNonlinSol_FixedPoint(nv_S, 50);

            ////STEP 8.3
            //    ARKStepSetNonlinearSolver(inner_mem, NLS); // NOTE: uncomment for implicit
            if(use_linear)
                MRIStepSetLinearSolver(mristep_mem, LS, NULL);
            else
                MRIStepSetNonlinearSolver(mristep_mem, NLS);

            ////STEP NINE
            MRIStepSetUserData(mristep_mem, (void *) &fast_userdata);  /* Pass udata to user functions */

            MRIStepSetPostInnerFn(mristep_mem, PostStoreStage);
            //    ARKStepSetPostprocessStepFn(inner_mem, PostStoreStage);
            //    MRIStepSetPreInnerFn(mristep_mem, StoreStage);
            MRIStepSetPostprocessStageFn(mristep_mem, ProcessStage);
        }

        //Set table
        ERKStepSetTable(arkode_mem, B);

        if(advance_mri) {
            if(advance_mri_test)
                MRIStepSetTable(mristep_mem, B->q, B);
            else
                MRIStepSetTable(mristep_mem, B2->q, B2);
        }

        // Free the Butcher table
        ARKodeButcherTable_Free(B);

        if(advance_erk)
        {
            // Use ERKStep to evolve state_old data (wrapped in nv_S) from t to tout=t+dt
            auto flag = ERKStepEvolve(arkode_mem, tout, nv_S, &t, ARK_NORMAL);
            // check flag status w/ helper function using Assert macro 0 = success, <0 = hard error, >0 = not-catastrophic warning
        }

        ////STEP ELEVEN
        if(advance_mri)
        {
            // Use MRIStep to evolve state_old data (wrapped in nv_S) from t to tout=t+dt
            auto flag = MRIStepEvolve(mristep_mem, tout, nv_S, &t, ARK_NORMAL);
            // check flag status w/ helper function
        }

        // Copy the result stored in nv_S to state_new
        for(int i=0; i<N_VGetNumSubvectors_ManyVector(nv_S); i++)
        {
            MultiFab::Copy(S_new[i], *amrex::sundials::getMFptr(N_VGetSubvector_ManyVector(nv_S, i)), 0, 0, S_new[i].nComp(), S_new[i].nGrow());
        }

        // Return timestep
        return timestep;
    }

    void time_interpolate (const T& S_new, const T& S_old, amrex::Real timestep_fraction, T& data)
    {
    }

    void map_data (std::function<void(T&)> Map)
    {
    }

};

}

#endif
