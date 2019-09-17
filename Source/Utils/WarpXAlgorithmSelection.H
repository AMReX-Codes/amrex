#ifndef UTILS_WARPXALGORITHMSELECTION_H_
#define UTILS_WARPXALGORITHMSELECTION_H_

#include <AMReX_ParmParse.H>
#include <string>

struct MaxwellSolverAlgo {
    enum {
        Yee = 0,
        CKC = 1
    };
};

struct ParticlePusherAlgo {
    enum {
        Boris = 0,
        Vay = 1
    };
};

struct CurrentDepositionAlgo {
    enum {
         Esirkepov = 0,
         Direct = 1
    };
};

struct ChargeDepositionAlgo {
    // Only the Standard algorithm is implemented
    enum {
         Standard = 0
    };
};

struct GatheringAlgo {
    // Only the Standard algorithm is implemented
    enum {
         Standard = 0
    };
};

int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key );

#endif // UTILS_WARPXALGORITHMSELECTION_H_
