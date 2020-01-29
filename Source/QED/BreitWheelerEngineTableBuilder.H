/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_breit_wheeler_engine_table_builder_h_
#define WARPX_breit_wheeler_engine_table_builder_h_

#include "QedWrapperCommons.H"
#include "BreitWheelerEngineInnards.H"

//This includes only the definition of a simple datastructure
//used to control the Breit Wheeler engine.
#include <breit_wheeler_engine_ctrl.h>

/**
 * A class which computes the lookup tables for the Breit Wheeler engine.
 */
class BreitWheelerEngineTableBuilder{
   public:
      /**
       * Computes the tables.
       * @param[in] ctrl control parameters to generate the tables
       * @param[out] innards structure holding both a copy of ctrl and lookup tables data
       */
      void compute_table
         (picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real> ctrl,
         BreitWheelerEngineInnards& innards) const;
};

#endif //WARPX_breit_wheeler_engine_table_builder_h_
