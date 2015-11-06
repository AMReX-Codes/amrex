#ifndef _ADR_F_H_
#define _ADR_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(FILCC,filcc)
                (const Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                 const int * domlo, const int * domhi,
                 const Real * dx_crse, const Real * xlo,
                 const int * bc);

BL_FORT_PROC_DECL(F_NETWORK_INIT,f_network_init)();

BL_FORT_PROC_DECL(GET_NUM_SPEC,get_num_spec)(int* nspec);

BL_FORT_PROC_DECL(GET_SPEC_NAMES,get_spec_names)(int* spec_names, int* ispec, int* len);

BL_FORT_PROC_DECL(GET_METHOD_PARAMS,get_method_params)(int* HYP_GROW);

BL_FORT_PROC_DECL(SET_METHOD_PARAMS,set_method_params)
    (const int& dm,
     const int& Density, const int& Xvel, 
     const int& FirstAdv, const int& FirstSpec, 
     const int& NumAdv, const int& normalize_species);

BL_FORT_PROC_DECL(SET_PROBLEM_PARAMS,set_problem_params)
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Outflow_value, 
     const int& Symmetry_value, const int& coord_type);

BL_FORT_PROC_DECL(INITDATA,initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const Real dx[], const Real xlo[], const Real xhi[]);

BL_FORT_PROC_DECL(STATE_ERROR,state_error)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const Real dx[], const Real xlo[], const Real problo[],
     const Real* time, const int* level);

BL_FORT_PROC_DECL(DENFILL,denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(XVELFILL,xvelfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(YVELFILL,yvelfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(ZVELFILL,zvelfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(SPECFILL,specfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

#ifdef REACTIONS
BL_FORT_PROC_DECL(REACTFILL,reactfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

BL_FORT_PROC_DECL(ESTDT,estdt)
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const Real dx[], Real* dt);

BL_FORT_PROC_DECL(ADVECT,advect)
    (const Real* time, 
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     D_DECL(const BL_FORT_FAB_ARG(udgnvx),
            const BL_FORT_FAB_ARG(udgnvy),
            const BL_FORT_FAB_ARG(udgnvz)),
     const BL_FORT_FAB_ARG(src),
     D_DECL(const BL_FORT_FAB_ARG(xflux),
            const BL_FORT_FAB_ARG(yflux),
            const BL_FORT_FAB_ARG(zflux)),
     const Real dx[], const Real* dt,
     const int& verbose);

BL_FORT_PROC_DECL(SUMMASS,summass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real dx[], Real* mass);

BL_FORT_PROC_DECL(ENFORCE_NONNEGATIVE_SPECIES,enforce_nonnegative_species)
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[]);

BL_FORT_PROC_DECL(EXT_SRC,ext_src)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(ext_src),
     const Real* dx, const Real* time, const Real* dt);

BL_FORT_PROC_DECL(REACT_STATE,react_state)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(Sold),
     const BL_FORT_FAB_ARG(Snew),
     const Real& time, const Real& dt_react);

BL_FORT_PROC_DECL(FILL_DIFF_COEFF,fill_diff_coeff)
    (const int lo[], const int hi[],
     D_DECL(const BL_FORT_FAB_ARG(xcoeffs),
            const BL_FORT_FAB_ARG(ycoeffs),
            const BL_FORT_FAB_ARG(zcoeffs)),
     const Real* dx);

BL_FORT_PROC_DECL(APPLY_METRIC,apply_metric)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     D_DECL(const BL_FORT_FAB_ARG(xedge),
            const BL_FORT_FAB_ARG(yedge),
            const BL_FORT_FAB_ARG(zedge)), 
            const Real* dx, const int* coord_type);

BL_FORT_PROC_DECL(UNWEIGHT_CC,unweight_cc)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const Real* dx, const int* coord_type);

#endif
