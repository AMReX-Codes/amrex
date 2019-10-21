#ifndef WARPX_breit_wheeler_engine_table_builder_h_
#define WARPX_breit_wheeler_engine_table_builder_h_

#include "QedWrapperCommons.H"
#include "BreitWheelerEngineInnards.H"

#include <breit_wheeler_engine_ctrl.h>

class BreitWheelerEngineTableBuilder{
   public:
      void compute_table
         (picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real> ctrl,
         BreitWheelerEngineInnards& innards) const;
};

#endif //WARPX_breit_wheeler_engine_table_builder_h_

