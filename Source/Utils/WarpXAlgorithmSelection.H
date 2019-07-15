#ifndef UTILS_WARPXALGORITHMSELECTION_H_
#define UTILS_WARPXALGORITHMSELECTION_H_

#include <AMReX_ParmParse.H>
#include <string>

struct MaxwellSolverAlgo {
    // These numbers corresponds to the algorithm code in WarpX's
    // `warpx_push_bvec` and `warpx_push_evec_f`
    enum {
        Yee = 0,
        CKC = 1
    };
};

struct ParticlePusherAlgo {
    // These numbers correspond to the algorithm code in WarpX's
    // `warpx_particle_pusher`
    enum {
        Boris = 0,
        Vay = 1
    };
};

struct CurrentDepositionAlgo {
    enum {
        // These numbers corresponds to the algorithm code in PICSAR's
        // `depose_jxjyjz_generic` and `depose_jxjyjz_generic_2d`
         Direct = 3,
         DirectVectorized = 2,
         EsirkepovNonOptimized = 1,
         Esirkepov = 0
    };
};

struct ChargeDepositionAlgo {
    // These numbers corresponds to the algorithm code in WarpX's
    // `warpx_charge_deposition` function
    enum {
         Vectorized = 0,
         Standard = 1
    };
};

struct GatheringAlgo {
    // These numbers corresponds to the algorithm code in PICSAR's
    // `geteb3d_energy_conserving_generic` function
    enum {
         Vectorized = 0,
         Standard = 1
    };
};

int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key );

#endif // UTILS_WARPXALGORITHMSELECTION_H_
