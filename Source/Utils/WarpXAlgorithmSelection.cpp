/* Copyright 2019-2020 Axel Huebl, David Grote, Luca Fedeli
 * Remi Lehe, Weiqun Zhang, Yinjian Zhao
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXAlgorithmSelection.H"

#include <algorithm>
#include <cstring>
#include <map>


// Define dictionary with correspondance between user-input strings,
// and corresponding integer for use inside the code

const std::map<std::string, int> maxwell_solver_algo_to_int = {
    {"yee",     MaxwellSolverAlgo::Yee },
#ifndef WARPX_DIM_RZ // Not available in RZ
    {"ckc",     MaxwellSolverAlgo::CKC },
#endif
    {"default", MaxwellSolverAlgo::Yee }
};

const std::map<std::string, int> particle_pusher_algo_to_int = {
    {"boris",   ParticlePusherAlgo::Boris },
    {"vay",     ParticlePusherAlgo::Vay },
    {"higuera", ParticlePusherAlgo::HigueraCary },
    {"default", ParticlePusherAlgo::Boris }
};

const std::map<std::string, int> current_deposition_algo_to_int = {
    {"esirkepov", CurrentDepositionAlgo::Esirkepov },
    {"direct",    CurrentDepositionAlgo::Direct },
#ifdef WARPX_USE_PSATD
    {"default",   CurrentDepositionAlgo::Direct }
#else
    {"default",   CurrentDepositionAlgo::Esirkepov }
#endif
};

const std::map<std::string, int> charge_deposition_algo_to_int = {
    {"standard",   ChargeDepositionAlgo::Standard },
    {"default",    ChargeDepositionAlgo::Standard }
};

const std::map<std::string, int> gathering_algo_to_int = {
    {"energy-conserving",   GatheringAlgo::EnergyConserving },
    {"momentum-conserving", GatheringAlgo::MomentumConserving },
    {"default",             GatheringAlgo::EnergyConserving }
};

const std::map<std::string, int> load_balance_costs_update_algo_to_int = {
    {"timers",    LoadBalanceCostsUpdateAlgo::Timers },
    {"heuristic", LoadBalanceCostsUpdateAlgo::Heuristic },
    {"default",   LoadBalanceCostsUpdateAlgo::Timers }
};


int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key ){

    // Read user input ; use "default" if it is not found
    std::string algo = "default";
    pp.query( pp_search_key, algo );
    // Convert to lower case
    std::transform(algo.begin(), algo.end(), algo.begin(), ::tolower);

    // Pick the right dictionary
    std::map<std::string, int> algo_to_int;
    if (0 == std::strcmp(pp_search_key, "maxwell_fdtd_solver")) {
        algo_to_int = maxwell_solver_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "particle_pusher")) {
        algo_to_int = particle_pusher_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "current_deposition")) {
        algo_to_int = current_deposition_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "charge_deposition")) {
        algo_to_int = charge_deposition_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "field_gathering")) {
        algo_to_int = gathering_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "load_balance_costs_update")) {
        algo_to_int = load_balance_costs_update_algo_to_int;
    } else {
        std::string pp_search_string = pp_search_key;
        amrex::Abort("Unknown algorithm type: " + pp_search_string);
    }

    // Check if the user-input is a valid key for the dictionary
    if (algo_to_int.count(algo) == 0) {
        // Not a valid key ; print error message
        std::string pp_search_string = pp_search_key;
        std::string error_message = "Invalid string for algo." + pp_search_string
            + ": " + algo + ".\nThe valid values are:\n";
        for ( const auto &valid_pair : algo_to_int ) {
            if (valid_pair.first != "default"){
                error_message += " - " + valid_pair.first + "\n";
            }
        }
        amrex::Abort(error_message);
    }

    // If the input is a valid key, return the value
    return algo_to_int[algo];
}
