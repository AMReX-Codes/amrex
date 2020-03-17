/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "QuantumSyncEngineTableBuilder.H"

//Include the full Quantum Synchrotron engine with table generation support
//(after some consistency tests). This requires to have a recent version
// of the Boost library.
#ifdef PXRMP_CORE_ONLY
    #error The Table Builder is incompatible with PXRMP_CORE_ONLY
#endif

#ifdef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
    #warning quantum_sync_engine.hpp should not have been included before reaching this point.
#endif
#include <quantum_sync_engine.hpp>
//_______________________________________________

//Some handy aliases
using PicsarQuantumSynchrotronEngine = picsar::multi_physics::
    quantum_synchrotron_engine<amrex::Real, QedUtils::DummyStruct>;

using PicsarQuantumSynchrotronCtrl =
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real>;
//_______________________________________________

void
QuantumSynchrotronEngineTableBuilder::compute_table
    (PicsarQuantumSynchrotronCtrl ctrl,
     QuantumSynchrotronEngineInnards& innards) const
{
    PicsarQuantumSynchrotronEngine qs_engine(
        std::move(QedUtils::DummyStruct()), 1.0, ctrl);

    qs_engine.compute_dN_dt_lookup_table();
    qs_engine.compute_cumulative_phot_em_table();

    auto qs_innards_picsar = qs_engine.export_innards();

    //Copy data in a GPU-friendly data-structure
    innards.ctrl = qs_innards_picsar.qs_ctrl;
    innards.KKfunc_coords.assign(qs_innards_picsar.KKfunc_table_coords_ptr,
        qs_innards_picsar.KKfunc_table_coords_ptr +
        qs_innards_picsar.KKfunc_table_coords_how_many);
    innards.KKfunc_data.assign(qs_innards_picsar.KKfunc_table_data_ptr,
        qs_innards_picsar.KKfunc_table_data_ptr +
        qs_innards_picsar.KKfunc_table_data_how_many);
        innards.cum_distrib_coords_1.assign(
    qs_innards_picsar.cum_distrib_table_coords_1_ptr,
        qs_innards_picsar.cum_distrib_table_coords_1_ptr +
        qs_innards_picsar.cum_distrib_table_coords_1_how_many);
    innards.cum_distrib_coords_2.assign(
        qs_innards_picsar.cum_distrib_table_coords_2_ptr,
        qs_innards_picsar.cum_distrib_table_coords_2_ptr +
        qs_innards_picsar.cum_distrib_table_coords_2_how_many);
    innards.cum_distrib_data.assign(
        qs_innards_picsar.cum_distrib_table_data_ptr,
        qs_innards_picsar.cum_distrib_table_data_ptr +
        qs_innards_picsar.cum_distrib_table_data_how_many);
    //____
}
