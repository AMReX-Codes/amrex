/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BreitWheelerEngineTableBuilder.H"

//Include the full Breit Wheeler engine with table generation support
//(after some consistency tests). This requires to have a recent version
// of the Boost library.
#ifdef PXRMP_CORE_ONLY
    #error The Table Builder is incompatible with PXRMP_CORE_ONLY
#endif

#ifdef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
    #warning breit_wheeler_engine.hpp should not have been included before reaching this point.
#endif
#include <breit_wheeler_engine.hpp>
//_______________________________________________

//Some handy aliases
using PicsarBreitWheelerEngine = picsar::multi_physics::
    breit_wheeler_engine<amrex::Real, QedUtils::DummyStruct>;

using PicsarBreitWheelerCtrl =
    picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real>;
//_______________________________________________

void
BreitWheelerEngineTableBuilder::compute_table
    (PicsarBreitWheelerCtrl ctrl,
     BreitWheelerEngineInnards& innards) const
{
    PicsarBreitWheelerEngine bw_engine(
        std::move(QedUtils::DummyStruct()), 1.0, ctrl);

    bw_engine.compute_dN_dt_lookup_table();
    bw_engine.compute_cumulative_pair_table();

    auto bw_innards_picsar = bw_engine.export_innards();

    //Copy data in a GPU-friendly data-structure
    innards.ctrl = bw_innards_picsar.bw_ctrl;
    innards.TTfunc_coords.assign(bw_innards_picsar.TTfunc_table_coords_ptr,
        bw_innards_picsar.TTfunc_table_coords_ptr +
        bw_innards_picsar.TTfunc_table_coords_how_many);
    innards.TTfunc_data.assign(bw_innards_picsar.TTfunc_table_data_ptr,
        bw_innards_picsar.TTfunc_table_data_ptr +
        bw_innards_picsar.TTfunc_table_data_how_many);
    innards.cum_distrib_coords_1.assign(
        bw_innards_picsar.cum_distrib_table_coords_1_ptr,
        bw_innards_picsar.cum_distrib_table_coords_1_ptr +
        bw_innards_picsar.cum_distrib_table_coords_1_how_many);
    innards.cum_distrib_coords_2.assign(
        bw_innards_picsar.cum_distrib_table_coords_2_ptr,
        bw_innards_picsar.cum_distrib_table_coords_2_ptr +
        bw_innards_picsar.cum_distrib_table_coords_2_how_many);
    innards.cum_distrib_data.assign(
        bw_innards_picsar.cum_distrib_table_data_ptr,
        bw_innards_picsar.cum_distrib_table_data_ptr +
        bw_innards_picsar.cum_distrib_table_data_how_many);
    //____
}
