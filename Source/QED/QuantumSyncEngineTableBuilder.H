#ifndef WARPX_quantum_sync_engine_table_builder_h_
#define WARPX_quantum_sync_engine_table_builder_h_

#include "QedWrapperCommons.H"
#include "QuantumSyncEngineInnards.H"

#include <quantum_sync_engine_ctrl.h>

class QuantumSynchrotronEngineTableBuilder{
public:
      void compute_table
         (picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real> ctrl,
         QuantumSynchrotronEngineInnards& innards) const;
};

#endif //WARPX_quantum_sync_engine_table_builder_h_
