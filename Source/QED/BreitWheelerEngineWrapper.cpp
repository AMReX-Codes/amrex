#include "BreitWheelerEngineWrapper.H"
//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

// Factory class =============================

BreitWheelerEngine::BreitWheelerEngine (){}

//Builds the functor to initialize the optical depth
BreitWheelerGetOpticalDepth BreitWheelerEngine::build_optical_depth_functor ()
{
    return BreitWheelerGetOpticalDepth();
}

//Builds the functor to evolve the optical depth
BreitWheelerEvolveOpticalDepth BreitWheelerEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerEvolveOpticalDepth(&innards);
}

//Builds the functor to generate the pairs
BreitWheelerGeneratePairs BreitWheelerEngine::build_pair_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerGeneratePairs(&innards);
}

bool BreitWheelerEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

/* \brief Reads lookup tables from 'file' on disk */
void BreitWheelerEngine::read_lookup_tables (std::string file)
{
    std::ifstream ifile(file, std::ios::in | std::ios::binary);

    //Header (control parameters)
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_min),
        sizeof(innards.ctrl.chi_phot_min));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_tdndt_min),
        sizeof(innards.ctrl.chi_phot_tdndt_min));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_tdndt_max),
        sizeof(innards.ctrl.chi_phot_tdndt_max));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_tdndt_how_many),
        sizeof(innards.ctrl.chi_phot_tdndt_how_many));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_tpair_min),
        sizeof(innards.ctrl.chi_phot_tpair_min));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_tpair_max),
        sizeof(innards.ctrl.chi_phot_tpair_max));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_phot_tpair_how_many),
        sizeof(innards.ctrl.chi_phot_tpair_how_many));
    ifile.read(reinterpret_cast<char*>(&innards.ctrl.chi_frac_tpair_how_many),
        sizeof(innards.ctrl.chi_frac_tpair_how_many));
    //_______

    //Data
    size_t size_buf = sizeof(amrex::Real)*innards.ctrl.chi_phot_tdndt_how_many;
    char* data_buf = new char(size_buf);
    ifile.read(data_buf, size_buf);
    innards.TTfunc_coords.assign((amrex::Real*)data_buf,
        (amrex::Real*)data_buf + size_buf);
    ifile.read(data_buf, size_buf);
    innards.TTfunc_data.assign((amrex::Real*)data_buf,
        (amrex::Real*)data_buf + size_buf);
    delete[] data_buf;
    //_______

    ifile.close();

    lookup_tables_initialized = true;
}


#ifdef WARPX_QED_TABLE_GEN
//Initializes the Lookup tables using the default settings
//provided by the library
void BreitWheelerEngine::compute_lookup_tables_default ()
{
    //A control parameters structure
    //with the default values provided by the library
    WarpXBreitWheelerWrapperCtrl ctrl_default;

    computes_lookup_tables(ctrl_default);

    lookup_tables_initialized = true;
}

// Computes the Lookup tables using user-defined settings
void BreitWheelerEngine::compute_custom_lookup_tables (
    WarpXBreitWheelerWrapperCtrl ctrl)
{
    computes_lookup_tables(ctrl);

    lookup_tables_initialized = true;
}


/* \brief Writes lookup tables on disk in 'file'
 *  return false if it fails. */
bool BreitWheelerEngine::write_lookup_tables (
        std::string file) const
{
    if(!lookup_tables_initialized)
        return false;

    std::ofstream of(file, std::ios::out | std::ios::binary);

    //Header (control parameters)
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_min),
        sizeof(innards.ctrl.chi_phot_min));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_tdndt_min),
        sizeof(innards.ctrl.chi_phot_tdndt_min));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_tdndt_max),
        sizeof(innards.ctrl.chi_phot_tdndt_max));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_tdndt_how_many),
        sizeof(innards.ctrl.chi_phot_tdndt_how_many));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_tpair_min),
        sizeof(innards.ctrl.chi_phot_tpair_min));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_tpair_max),
        sizeof(innards.ctrl.chi_phot_tpair_max));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_phot_tpair_how_many),
        sizeof(innards.ctrl.chi_phot_tpair_how_many));
    of.write(reinterpret_cast<const char*>(&innards.ctrl.chi_frac_tpair_how_many),
        sizeof(innards.ctrl.chi_frac_tpair_how_many));
    //_______

    //Data
    of.write(reinterpret_cast<const char*>(innards.TTfunc_coords.dataPtr()),
        sizeof(amrex::Real)*innards.TTfunc_coords.size());
    of.write(reinterpret_cast<const char*>(innards.TTfunc_data.dataPtr()),
        sizeof(amrex::Real)*innards.TTfunc_data.size());
    // TODO: add other table
    //_______

    of.close();

    return true;
}

//Private function which actually computes the lookup tables
void BreitWheelerEngine::computes_lookup_tables (
    WarpXBreitWheelerWrapperCtrl ctrl)
{
    //Lambda is not actually used if S.I. units are enabled
    WarpXBreitWheelerWrapper bw_engine(
        std::move(QedDummyStruct()), 1.0, ctrl);

    bw_engine.compute_dN_dt_lookup_table();
    //bw_engine.compute_cumulative_pair_table();

    auto bw_innards_picsar = bw_engine.export_innards();

    //Copy data in a GPU-friendly data-structure
    innards.ctrl = bw_innards_picsar.bw_ctrl;
    innards.TTfunc_coords.assign(bw_innards_picsar.TTfunc_table_coords_ptr,
        bw_innards_picsar.TTfunc_table_coords_ptr +
        bw_innards_picsar.TTfunc_table_coords_how_many);
    innards.TTfunc_data.assign(bw_innards_picsar.TTfunc_table_data_ptr,
        bw_innards_picsar.TTfunc_table_data_ptr +
        bw_innards_picsar.TTfunc_table_data_how_many);
    //____
}

#endif

//============================================
