#include "QuantumSyncEngineWrapper.H"
//This file provides a wrapper aroud the quantum_sync engine
//provided by the PICSAR library

// Factory class =============================

QuantumSynchrotronEngine::QuantumSynchrotronEngine (){}

//Builds the functor to evolve the optical depth
QuantumSynchrotronGetOpticalDepth
QuantumSynchrotronEngine::build_optical_depth_functor ()
{
    return QuantumSynchrotronGetOpticalDepth();
}

//Builds the functor to evolve the optical depth
QuantumSynchrotronEvolveOpticalDepth QuantumSynchrotronEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return QuantumSynchrotronEvolveOpticalDepth(&innards);
}

bool QuantumSynchrotronEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

/* \brief Reads lookup tables from 'file' on disk */
void QuantumSynchrotronEngine::read_lookup_tables (std::string file)
{
    std::ifstream ifile(file, std::ios::in | std::ios::binary);

    //Header (control parameters)
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_min),
        sizeof(innards.ctrl.chi_part_min));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_tdndt_min),
        sizeof(innards.ctrl.chi_part_tdndt_min));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_tdndt_max),
        sizeof(innards.ctrl.chi_part_tdndt_max));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_tdndt_how_many),
        sizeof(innards.ctrl.chi_part_tdndt_how_many));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_tem_min),
        sizeof(innards.ctrl.chi_part_tem_min));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_tem_max),
        sizeof(innards.ctrl.chi_part_tem_max));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_part_tem_how_many),
        sizeof(innards.ctrl.chi_part_tem_how_many));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.prob_tem_how_many),
        sizeof(innards.ctrl.prob_tem_how_many));
    //_______

    //Data
    size_t size_buf = sizeof(amrex::Real)*innards.ctrl.chi_part_tdndt_how_many;
    char* data_buf = new char(size_buf);
    ifile.read(data_buf, size_buf);
    innards.KKfunc_coords.assign((amrex::Real*)data_buf,
        (amrex::Real*)data_buf + size_buf);
    ifile.read(data_buf, size_buf);
    innards.KKfunc_data.assign((amrex::Real*)data_buf,
        (amrex::Real*)data_buf + size_buf);
    delete[] data_buf;
    //_______

    ifile.close();

    lookup_tables_initialized = true;
}

#ifdef WARPX_QED_TABLE_GEN

//Initializes the Lookup tables using the default settings
//provided by the library
void QuantumSynchrotronEngine::compute_lookup_tables_default ()
{
    //A control parameters structure
    //with the default values provided by the library
    WarpXQuantumSynchrotronWrapperCtrl ctrl_default;

    computes_lookup_tables(ctrl_default);

    lookup_tables_initialized = true;
}

// Computes the Lookup tables using user-defined settings
void QuantumSynchrotronEngine::compute_custom_lookup_tables (
    WarpXQuantumSynchrotronWrapperCtrl ctrl)
{
    computes_lookup_tables(ctrl);

    lookup_tables_initialized = true;
}


/* \brief Writes lookup tables on disk in 'file'
 *  return false if it fails. */
bool QuantumSynchrotronEngine::write_lookup_tables (
        std::string file) const
{
    if(!lookup_tables_initialized)
        return false;

    std::ofstream of(file, std::ios::out | std::ios::binary);

    //Header (control parameters)
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_min),
        sizeof(innards.ctrl.chi_part_min));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_tdndt_min),
        sizeof(innards.ctrl.chi_part_tdndt_min));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_tdndt_max),
        sizeof(innards.ctrl.chi_part_tdndt_max));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_tdndt_how_many),
        sizeof(innards.ctrl.chi_part_tdndt_how_many));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_tem_min),
        sizeof(innards.ctrl.chi_part_tem_max));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_tem_max),
        sizeof(innards.ctrl.chi_part_tem_max));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_part_tem_how_many),
        sizeof(innards.ctrl.chi_part_tem_how_many));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.prob_tem_how_many),
        sizeof(innards.ctrl.prob_tem_how_many));
    //_______

    //Data
    of.write(reinterpret_cast<const char*>(innards.KKfunc_coords.dataPtr()),
        sizeof(amrex::Real)*innards.KKfunc_coords.size());
    of.write(reinterpret_cast<const char*>(innards.KKfunc_data.dataPtr()),
        sizeof(amrex::Real)*innards.KKfunc_data.size());
    // TODO: add other table
    //_______

    of.close();

    return true;
}

    //Private function which actually computes the lookup tables
void QuantumSynchrotronEngine::computes_lookup_tables (
    WarpXQuantumSynchrotronWrapperCtrl ctrl)
{
    //Lambda is not actually used if S.I. units are enabled
    WarpXQuantumSynchrotronWrapper qs_engine(
        std::move(QedUtils::DummyStruct()), 1.0, ctrl);

    qs_engine.compute_dN_dt_lookup_table();
    //qs_engine.compute_cumulative_pair_table();

    auto qs_innards_picsar = qs_engine.export_innards();

    //Copy data in a GPU-friendly data-structure
    innards.ctrl = qs_innards_picsar.qs_ctrl;
    innards.KKfunc_coords.assign(qs_innards_picsar.KKfunc_table_coords_ptr,
        qs_innards_picsar.KKfunc_table_coords_ptr +
        qs_innards_picsar.KKfunc_table_coords_how_many);
    innards.KKfunc_data.assign(qs_innards_picsar.KKfunc_table_data_ptr,
        qs_innards_picsar.KKfunc_table_data_ptr +
        qs_innards_picsar.KKfunc_table_data_how_many);
    //____
}
#endif

//============================================
