/* Copyright 2019 David Grote, Luca Fedeli, Remi Lehe
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
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
        Vay = 1,
        HigueraCary = 2
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
    enum {
         EnergyConserving = 0,
         MomentumConserving
    };
};

/** Strategy to compute weights for use in load balance.
 */
struct LoadBalanceCostsUpdateAlgo {
    enum {
        Timers = 0,   //!< load balance according to in-code timer-based weights (i.e., with  `costs`)
        Heuristic = 1 /**< load balance according to weights computed from number of cells
                             and number of particles per box (i.e., with `costs_heuristic`)*/
    };
};

int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key );

#endif // UTILS_WARPXALGORITHMSELECTION_H_
