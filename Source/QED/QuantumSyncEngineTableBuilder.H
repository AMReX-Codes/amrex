/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_quantum_sync_engine_table_builder_h_
#define WARPX_quantum_sync_engine_table_builder_h_

#include "QedWrapperCommons.H"
#include "QuantumSyncEngineInnards.H"

//This includes only the definition of a simple datastructure
//used to control the Quantum Synchrotron engine.
#include <quantum_sync_engine_ctrl.h>

/**
 * A class which computes the lookup tables for the Quantum Synchrotron engine.
 */
class QuantumSynchrotronEngineTableBuilder{
public:
      /**
       * Computes the tables.
       * @param[in] ctrl control parameters to generate the tables
       * @param[out] innards structure holding both a copy of ctrl and lookup tables data
       */
      void compute_table
         (picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real> ctrl,
         QuantumSynchrotronEngineInnards& innards) const;
};

#endif //WARPX_quantum_sync_engine_table_builder_h_
