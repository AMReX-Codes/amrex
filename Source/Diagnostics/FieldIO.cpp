/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote, Maxence Thevenet,
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang, Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldIO.H"
#include "WarpX.H"

#include <AMReX_FillPatchUtil_F.H>
#include <AMReX_Interpolater.H>

#ifdef WARPX_USE_OPENPMD
#include <openPMD/openPMD.hpp>
#endif
#ifdef WARPX_USE_PSATD
#include <SpectralSolver.H>
#endif


using namespace amrex;

namespace
{
    // Return true if any element in elems is in vect, false otherwise
    static bool is_in_vector(std::vector<std::string> vect,
                             std::vector<std::string> elems)
    {
        bool value = false;
        for (std::string elem : elems){
            if (std::find(vect.begin(), vect.end(), elem) != vect.end())
                value = true;
        }
        return value;
    }
}

#ifdef WARPX_USE_OPENPMD
/** \brief For a given field that is to be written to an openPMD file,
 * set the metadata that indicates the physical unit.
 */
void
setOpenPMDUnit( openPMD::Mesh mesh, const std::string field_name )
{
    if (field_name[0] == 'E'){  // Electric field
        mesh.setUnitDimension({
            {openPMD::UnitDimension::L,  1},
            {openPMD::UnitDimension::M,  1},
            {openPMD::UnitDimension::T, -3},
            {openPMD::UnitDimension::I, -1},
        });
    } else if (field_name[0] == 'B'){ // Magnetic field
        mesh.setUnitDimension({
            {openPMD::UnitDimension::M,  1},
            {openPMD::UnitDimension::I, -1},
            {openPMD::UnitDimension::T, -2}
        });
    } else if (field_name[0] == 'j'){ // current
        mesh.setUnitDimension({
            {openPMD::UnitDimension::L, -2},
            {openPMD::UnitDimension::I,  1},
        });
    } else if (field_name.substr(0,3) == "rho"){ // charge density
        mesh.setUnitDimension({
            {openPMD::UnitDimension::L, -3},
            {openPMD::UnitDimension::I,  1},
            {openPMD::UnitDimension::T,  1},
        });
    }
}


/** \brief
 * Convert an IntVect to a std::vector<std::uint64_t>
 * and reverse the order of the elements
 * (used for compatibility with the openPMD API)
 */
std::vector<std::uint64_t>
getReversedVec( const IntVect& v )
{
  // Convert the IntVect v to and std::vector u
  std::vector<std::uint64_t> u = {
    AMREX_D_DECL(
                 static_cast<std::uint64_t>(v[0]),
                 static_cast<std::uint64_t>(v[1]),
                 static_cast<std::uint64_t>(v[2])
                 )
  };
  // Reverse the order of elements, if v corresponds to the indices of a
  // Fortran-order array (like an AMReX FArrayBox)
  // but u is intended to be used with a C-order API (like openPMD)
  std::reverse( u.begin(), u.end() );

  return u;
}

/** \brief
 * Convert Real* pointer to a std::vector<double>,
 * and reverse the order of the elements
 * (used for compatibility with the openPMD API)
 */
std::vector<double>
getReversedVec( const Real* v )
{
  // Convert Real* v to and std::vector u
  std::vector<double> u = {
    AMREX_D_DECL(
                 static_cast<double>(v[0]),
                 static_cast<double>(v[1]),
                 static_cast<double>(v[2])
                 )
  };
  // Reverse the order of elements, if v corresponds to the indices of a
  // Fortran-order array (like an AMReX FArrayBox)
  // but u is intended to be used with a C-order API (like openPMD)
  std::reverse( u.begin(), u.end() );

  return u;
}

/** \brief Write the `ncomp` components of `mf` (with names `varnames`)
 * into a file `filename` in openPMD format.
 **/
void
WriteOpenPMDFields( const std::string& filename,
                  const std::vector<std::string>& varnames,
                  const MultiFab& mf, const Geometry& geom,
                  const int iteration, const double time )
{
  WARPX_PROFILE("WriteOpenPMDFields()");

  const int ncomp = mf.nComp();

  // Create a few vectors that store info on the global domain
  // Swap the indices for each of them, since AMReX data is Fortran order
  // and since the openPMD API assumes contiguous C order
  // - Size of the box, in integer number of cells
  const Box& global_box = geom.Domain();
  auto global_size = getReversedVec(global_box.size());
  // - Grid spacing
  std::vector<double> grid_spacing = getReversedVec(geom.CellSize());
  // - Global offset
  std::vector<double> global_offset = getReversedVec(geom.ProbLo());
  // - AxisLabels
#if AMREX_SPACEDIM==3
  std::vector<std::string> axis_labels{"x", "y", "z"};
#else
  std::vector<std::string> axis_labels{"x", "z"};
#endif

  // Prepare the type of dataset that will be written
  openPMD::Datatype datatype = openPMD::determineDatatype<Real>();
  auto dataset = openPMD::Dataset(datatype, global_size);

  // Create new file and store the time/iteration info
  auto series = [filename](){
      if( ParallelDescriptor::NProcs() > 1 )
          return openPMD::Series( filename,
                                  openPMD::AccessType::CREATE,
                                  ParallelDescriptor::Communicator() );
      else
          return openPMD::Series( filename,
                                  openPMD::AccessType::CREATE );
  }();

  auto series_iteration = series.iterations[iteration];
  series_iteration.setTime( time );

  // Loop through the different components, i.e. different fields stored in mf
  for (int icomp=0; icomp<ncomp; icomp++){

    // Check if this field is a vector or a scalar, and extract the field name
    const std::string& varname = varnames[icomp];
    std::string field_name = varname;
    std::string comp_name = openPMD::MeshRecordComponent::SCALAR;
    bool is_vector = false;
    for (const char* vector_field: {"E", "B", "j"}){
        for (const char* comp: {"x", "y", "z"}){
            if (varname[0] == *vector_field && varname[1] == *comp ){
                is_vector = true;
                field_name = varname[0] + varname.substr(2); // Strip component
                comp_name = varname[1];
            }
        }
    }

    // Setup the mesh record accordingly
    auto mesh = series_iteration.meshes[field_name];
    mesh.setDataOrder(openPMD::Mesh::DataOrder::F); // MultiFab: Fortran order
    mesh.setAxisLabels( axis_labels );
    mesh.setGridSpacing( grid_spacing );
    mesh.setGridGlobalOffset( global_offset );
    setOpenPMDUnit( mesh, field_name );

    // Create a new mesh record component, and store the associated metadata
    auto mesh_comp = mesh[comp_name];
    mesh_comp.resetDataset( dataset );
    // Cell-centered data: position is at 0.5 of a cell size.
    mesh_comp.setPosition(std::vector<double>{AMREX_D_DECL(0.5, 0.5, 0.5)});

    // Loop through the multifab, and store each box as a chunk,
    // in the openPMD file.
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

      const FArrayBox& fab = mf[mfi];
      const Box& local_box = fab.box();

      // Determine the offset and size of this chunk
      IntVect box_offset = local_box.smallEnd() - global_box.smallEnd();
      auto chunk_offset = getReversedVec(box_offset);
      auto chunk_size = getReversedVec(local_box.size());

      // Write local data
      Real const * local_data = fab.dataPtr(icomp);
      mesh_comp.storeChunk(openPMD::shareRaw(local_data),
                           chunk_offset, chunk_size);
    }
  }
  // Flush data to disk after looping over all components
  series.flush();
}
#endif // WARPX_USE_OPENPMD

#ifdef WARPX_DIM_RZ
void
ConstructTotalRZField(std::array< std::unique_ptr<MultiFab>, 3 >& mf_total,
                      const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field)
{
    // Sum over the real components, giving quantity at theta=0
    MultiFab::Copy(*mf_total[0], *vector_field[0], 0, 0, 1, vector_field[0]->nGrowVect());
    MultiFab::Copy(*mf_total[1], *vector_field[1], 0, 0, 1, vector_field[1]->nGrowVect());
    MultiFab::Copy(*mf_total[2], *vector_field[2], 0, 0, 1, vector_field[2]->nGrowVect());
    for (int ic=1 ; ic < vector_field[0]->nComp() ; ic += 2) {
        MultiFab::Add(*mf_total[0], *vector_field[0], ic, 0, 1, vector_field[0]->nGrowVect());
        MultiFab::Add(*mf_total[1], *vector_field[1], ic, 0, 1, vector_field[1]->nGrowVect());
        MultiFab::Add(*mf_total[2], *vector_field[2], ic, 0, 1, vector_field[2]->nGrowVect());
    }
}
#endif

void
PackPlotDataPtrs (Vector<const MultiFab*>& pmf,
                  const std::array<std::unique_ptr<MultiFab>,3>& data)
{
    BL_ASSERT(pmf.size() == AMREX_SPACEDIM);
#if (AMREX_SPACEDIM == 3)
    pmf[0] = data[0].get();
    pmf[1] = data[1].get();
    pmf[2] = data[2].get();
#elif (AMREX_SPACEDIM == 2)
    pmf[0] = data[0].get();
    pmf[1] = data[2].get();
#endif
}

/** \brief Takes an array of 3 MultiFab `vector_field`
 * (representing the x, y, z components of a vector),
 * averages it to the cell center, and stores the
 * resulting MultiFab in mf_avg (in the components dcomp to dcomp+2)
 */
void
AverageAndPackVectorField( MultiFab& mf_avg,
                           const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
                           const DistributionMapping& dm,
                           const int dcomp, const int ngrow )
{
#ifndef WARPX_DIM_RZ
    (void)dm;
#endif
    // The object below is temporary, and is needed because
    // `average_edge_to_cellcenter` requires fields to be passed as Vector
    Vector<const MultiFab*> srcmf(AMREX_SPACEDIM);

    // Check the type of staggering of the 3-component `vector_field`
    // and average accordingly:
    // - Fully cell-centered field (no average needed; simply copy)
    if ( vector_field[0]->is_cell_centered() ){

        MultiFab::Copy( mf_avg, *vector_field[0], 0, dcomp  , 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[1], 0, dcomp+1, 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[2], 0, dcomp+2, 1, ngrow);

        // - Fully nodal
    } else if ( vector_field[0]->is_nodal() ){

        amrex::average_node_to_cellcenter( mf_avg, dcomp  ,
                                          *vector_field[0], 0, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                          *vector_field[1], 0, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+2,
                                          *vector_field[2], 0, 1, ngrow);

        // - Face centered, in the same way as B on a Yee grid
    } else if ( vector_field[0]->is_nodal(0) ){

        // Note that average_face_to_cellcenter operates only on the number of
        // arrays equal to the number of dimensions. So, for 2D, PackPlotDataPtrs
        // packs in the x and z (or r and z) arrays, which are then cell averaged.
        // The Copy code then copies the z from the 2nd to the 3rd field,
        // and copies over directly the y (or theta) component (which is
        // already cell centered).
        if (vector_field[0]->nComp() > 1) {
#ifdef WARPX_DIM_RZ
            // When there are more than one components, the total
            // fields needs to be constructed in temporary MultiFabs.
            // Note that mf_total is declared in the same way as
            // vector_field so that it can be passed into PackPlotDataPtrs.
            std::array<std::unique_ptr<MultiFab>,3> mf_total;
            mf_total[0].reset(new MultiFab(vector_field[0]->boxArray(), dm, 1, vector_field[0]->nGrowVect()));
            mf_total[1].reset(new MultiFab(vector_field[1]->boxArray(), dm, 1, vector_field[1]->nGrowVect()));
            mf_total[2].reset(new MultiFab(vector_field[2]->boxArray(), dm, 1, vector_field[2]->nGrowVect()));
            ConstructTotalRZField(mf_total, vector_field);
            PackPlotDataPtrs(srcmf, mf_total);
            amrex::average_face_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            MultiFab::Copy( mf_avg, *mf_total[1], 0, dcomp+1, 1, ngrow);
#else
           amrex::Abort("AverageAndPackVectorField not implemented for ncomp > 1");
#endif
        } else {
            PackPlotDataPtrs(srcmf, vector_field);
            amrex::average_face_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
#if (AMREX_SPACEDIM == 2)
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            MultiFab::Copy( mf_avg, *vector_field[1], 0, dcomp+1, 1, ngrow);
#endif
        }

        // - Edge centered, in the same way as E on a Yee grid
    } else if ( !vector_field[0]->is_nodal(0) ){

        // See comment above, though here, the y (or theta) component
        // has node centering.
        if (vector_field[0]->nComp() > 1) {
#ifdef WARPX_DIM_RZ
            // When there are more than one components, the total
            // fields needs to be constructed in temporary MultiFabs
            // Note that mf_total is declared in the same way as
            // vector_field so that it can be passed into PackPlotDataPtrs.
            std::array<std::unique_ptr<MultiFab>,3> mf_total;
            mf_total[0].reset(new MultiFab(vector_field[0]->boxArray(), dm, 1, vector_field[0]->nGrowVect()));
            mf_total[1].reset(new MultiFab(vector_field[1]->boxArray(), dm, 1, vector_field[1]->nGrowVect()));
            mf_total[2].reset(new MultiFab(vector_field[2]->boxArray(), dm, 1, vector_field[2]->nGrowVect()));
            ConstructTotalRZField(mf_total, vector_field);
            PackPlotDataPtrs(srcmf, mf_total);
            amrex::average_edge_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                              *mf_total[1], 0, 1, ngrow);
#else
           amrex::Abort("AverageAndPackVectorField not implemented for ncomp > 1");
#endif
        } else {
            PackPlotDataPtrs(srcmf, vector_field);
            amrex::average_edge_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
#if (AMREX_SPACEDIM == 2)
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                              *vector_field[1], 0, 1, ngrow);
#endif
        }

    } else {
        amrex::Abort("Unknown staggering.");
    }
}

/** \brief Takes all of the components of the three fields and
 * averages and packs them into the MultiFab mf_avg.
 */
void
AverageAndPackVectorFieldComponents (MultiFab& mf_avg,
                                     const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
                                     const DistributionMapping& dm,
                                     int& dcomp, const int ngrow )
{
    if (vector_field[0]->nComp() > 1) {
        std::array<std::unique_ptr<MultiFab>,3> vector_field_component;
        for (int icomp=0 ; icomp < vector_field[0]->nComp() ; icomp++) {
            vector_field_component[0].reset(new MultiFab(*vector_field[0], amrex::make_alias, icomp, 1));
            vector_field_component[1].reset(new MultiFab(*vector_field[1], amrex::make_alias, icomp, 1));
            vector_field_component[2].reset(new MultiFab(*vector_field[2], amrex::make_alias, icomp, 1));
            AverageAndPackVectorField(mf_avg, vector_field_component, dm, dcomp, ngrow);
            dcomp += 3;
        }
    }
}

/** \brief Take a MultiFab `scalar_field`
 * averages it to the cell center, and stores the
 * resulting MultiFab in mf_avg (in the components dcomp)
 */
void
AverageAndPackScalarField( MultiFab& mf_avg,
                           const MultiFab & scalar_field,
                           const int dcomp, const int ngrow )
{
    // Check the type of staggering of the 3-component `vector_field`
    // and average accordingly:
    // - Fully cell-centered field (no average needed; simply copy)
    if ( scalar_field.is_cell_centered() ){
        MultiFab::Copy( mf_avg, scalar_field, 0, dcomp, 1, ngrow);
        // - Fully nodal
    } else if ( scalar_field.is_nodal() ){

        amrex::average_node_to_cellcenter( mf_avg, dcomp, scalar_field, 0, 1, ngrow);

    } else {
        amrex::Abort("Unknown staggering.");
    }
}

/** \brief Takes the specified component of the scalar and
 * averages and packs it into the MultiFab mf_avg.
 */
void
AverageAndPackScalarFieldComponent (MultiFab& mf_avg,
                                    const MultiFab& scalar_field,
                                    const int icomp,
                                    const int dcomp, const int ngrow )
{
    MultiFab scalar_field_component(scalar_field, amrex::make_alias, icomp, 1);
    AverageAndPackScalarField(mf_avg, scalar_field_component, dcomp, ngrow);
}

/** \brief Generate mode variable name
 */
std::string
ComponentName(std::string fieldname, int mode, std::string type)
{
    if (type == "real") {
        return fieldname + "_" + std::to_string(mode) + "_" + "real";
    } else if (type == "imag") {
        return fieldname + "_" + std::to_string(mode) + "_" + "imag";
    } else {
        AMREX_ALWAYS_ASSERT( false );
    }
    // This should never be done
    return "";
}

/* \brief Copy vector field component data into the MultiFab that will be written out
 */
void
CopyVectorFieldComponentsToMultiFab (int lev, amrex::Vector<MultiFab>& mf_avg, MultiFab& mf_tmp,
                                     int icomp, int& dcomp, int ngrow,
                                     std::string fieldname, Vector<std::string>& varnames)
{
    if (mf_tmp.nComp() > 3) {
        if (lev==0) varnames.push_back(ComponentName(fieldname, 0, "real"));
        MultiFab::Copy( mf_avg[lev], mf_tmp, 3+icomp, dcomp++, 1, ngrow);
        int const nmodes = mf_tmp.nComp()/6;
        for (int mode=1 ; mode < nmodes ; mode++) {
            if (lev==0) varnames.push_back(ComponentName(fieldname, mode, "real"));
            MultiFab::Copy( mf_avg[lev], mf_tmp, 3*2*mode+icomp, dcomp++, 1, ngrow);
            if (lev==0) varnames.push_back(ComponentName(fieldname, mode, "imag"));
            MultiFab::Copy( mf_avg[lev], mf_tmp, 3*2*mode+3+icomp, dcomp++, 1, ngrow);
        }
    }
}

/* \brief Copy scalar field component data into the MultiFab that will be written out
 */
void
CopyScalarFieldComponentsToMultiFab (int lev, amrex::Vector<MultiFab>& mf_avg, MultiFab& mf_tmp,
                                     int& dcomp, int ngrow, int n_rz_azimuthal_modes,
                                     std::string fieldname, Vector<std::string>& varnames)
{
    if (n_rz_azimuthal_modes > 1) {
        if (lev==0) varnames.push_back(ComponentName(fieldname, 0, "real"));
        AverageAndPackScalarFieldComponent(mf_avg[lev], mf_tmp, 0, dcomp++, ngrow);
        for (int mode=1 ; mode < n_rz_azimuthal_modes ; mode++) {
            if (lev==0) varnames.push_back(ComponentName(fieldname, mode, "real"));
            AverageAndPackScalarFieldComponent(mf_avg[lev], mf_tmp, 2*mode-1, dcomp++, ngrow);
            if (lev==0) varnames.push_back(ComponentName(fieldname, mode, "imag"));
            AverageAndPackScalarFieldComponent(mf_avg[lev], mf_tmp, 2*mode  , dcomp++, ngrow);
        }
    }
}

/** \brief Add variable names to the list.
 */
void
AddToVarNames (Vector<std::string>& varnames,
               std::string name, std::string suffix) {
    auto coords = {"x", "y", "z"};
    for(auto coord:coords) varnames.push_back(name+coord+suffix);
}

/** \brief Add RZ variable names to the list.
 */
void
AddToVarNamesRZ (Vector<std::string>& varnames,
                 std::string name, std::string suffix) {
    auto coords = {"r", "theta", "z"};
    for(auto coord:coords) varnames.push_back(name+coord+suffix);
}

/** \brief Write the different fields that are meant for output,
 * into the vector of MultiFab `mf_avg` (one MultiFab per level)
 * after averaging them to the cell centers.
 */
void
WarpX::AverageAndPackFields ( Vector<std::string>& varnames,
                              amrex::Vector<MultiFab>& mf_avg, const int ngrow) const
{
    // Count how many different fields should be written (ncomp)
    MultiFab* cost = WarpX::getCosts(0);
    int ncomp = fields_to_plot.size()
        + static_cast<int>(plot_finepatch)*6
        + static_cast<int>(plot_crsepatch)*6
        + static_cast<int>(cost != nullptr and plot_costs);

    // Add in the RZ modes
    if (n_rz_azimuthal_modes > 1) {
        for (std::string field : fields_to_plot) {
            if (is_in_vector({"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz", "rho", "F"}, {field})) {
                ncomp += 2*n_rz_azimuthal_modes - 1;
            }
        }
        if (plot_finepatch) {
            ncomp += 6*(2*n_rz_azimuthal_modes - 1);
        }
    }

    int nvecs = 3;
    if (n_rz_azimuthal_modes > 1) {
        nvecs += 3*(2*n_rz_azimuthal_modes - 1);
    }

    // Loop over levels of refinement
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Allocate pointers to the `ncomp` fields that will be added
        mf_avg.push_back( MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

        // For E, B and J, if at least one component is requested,
        // build cell-centered temporary MultiFab with 3 comps
        MultiFab mf_tmp_E, mf_tmp_B, mf_tmp_J;
        // Build mf_tmp_E is at least one component of E is requested
        if (is_in_vector(fields_to_plot, {"Ex", "Ey", "Ez"} )){
            // Allocate temp MultiFab with 3 components
            mf_tmp_E = MultiFab(grids[lev], dmap[lev], nvecs, ngrow);
            // Fill MultiFab mf_tmp_E with averaged E
            AverageAndPackVectorField(mf_tmp_E, Efield_aux[lev], dmap[lev], 0, ngrow);
            int dcomp = 3;
            AverageAndPackVectorFieldComponents(mf_tmp_E, Efield_aux[lev], dmap[lev], dcomp, ngrow);
        }
        // Same for B
        if (is_in_vector(fields_to_plot, {"Bx", "By", "Bz"} )){
            mf_tmp_B = MultiFab(grids[lev], dmap[lev], nvecs, ngrow);
            AverageAndPackVectorField(mf_tmp_B, Bfield_aux[lev], dmap[lev], 0, ngrow);
            int dcomp = 3;
            AverageAndPackVectorFieldComponents(mf_tmp_B, Bfield_aux[lev], dmap[lev], dcomp, ngrow);
        }
        // Same for J
        if (is_in_vector(fields_to_plot, {"jx", "jy", "jz"} )){
            mf_tmp_J = MultiFab(grids[lev], dmap[lev], nvecs, ngrow);
            AverageAndPackVectorField(mf_tmp_J, current_fp[lev], dmap[lev], 0, ngrow);
            int dcomp = 3;
            AverageAndPackVectorFieldComponents(mf_tmp_J, current_fp[lev], dmap[lev], dcomp, ngrow);
        }

        int dcomp;
        // Go through the different fields in fields_to_plot, pack them into
        // mf_avg[lev] add the corresponding names to `varnames`.
        // plot_fine_patch and plot_coarse_patch are treated separately
        // (after this for loop).
        dcomp = 0;
        for (int ifield=0; ifield<fields_to_plot.size(); ifield++){
            std::string fieldname = fields_to_plot[ifield];
            if(lev==0) varnames.push_back(fieldname);
            if        (fieldname == "Ex"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_E, 0, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_E, 0, dcomp, ngrow, "Er", varnames);
            } else if (fieldname == "Ey"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_E, 1, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_E, 1, dcomp, ngrow, "Etheta", varnames);
            } else if (fieldname == "Ez"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_E, 2, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_E, 2, dcomp, ngrow, "Ez", varnames);
            } else if (fieldname == "Bx"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_B, 0, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_B, 0, dcomp, ngrow, "Br", varnames);
            } else if (fieldname == "By"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_B, 1, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_B, 1, dcomp, ngrow, "Btheta", varnames);
            } else if (fieldname == "Bz"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_B, 2, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_B, 2, dcomp, ngrow, "Bz", varnames);
            } else if (fieldname == "jx"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_J, 0, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_J, 0, dcomp, ngrow, "jr", varnames);
            } else if (fieldname == "jy"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_J, 1, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_J, 1, dcomp, ngrow, "jtheta", varnames);
            } else if (fieldname == "jz"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_J, 2, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_J, 2, dcomp, ngrow, "jz", varnames);
            } else if (fieldname == "rho"){
                AverageAndPackScalarField( mf_avg[lev], *rho_fp[lev], dcomp++, ngrow );
                CopyScalarFieldComponentsToMultiFab(lev, mf_avg, *rho_fp[lev], dcomp, ngrow, n_rz_azimuthal_modes,
                                                    fieldname, varnames);
            } else if (fieldname == "F"){
                AverageAndPackScalarField( mf_avg[lev], *F_fp[lev], dcomp++, ngrow );
                CopyScalarFieldComponentsToMultiFab(lev, mf_avg, *F_fp[lev], dcomp, ngrow, n_rz_azimuthal_modes,
                                                    fieldname, varnames);
            } else if (fieldname == "part_per_cell") {
                MultiFab temp_dat(grids[lev],mf_avg[lev].DistributionMap(),1,0);
                temp_dat.setVal(0);
                // MultiFab containing number of particles in each cell
                mypc->Increment(temp_dat, lev);
                AverageAndPackScalarField( mf_avg[lev], temp_dat, dcomp++, ngrow );
            } else if (fieldname == "part_per_grid"){
                const Vector<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);
                // MultiFab containing number of particles per grid
                // (stored as constant for all cells in each grid)
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(mf_avg[lev]); mfi.isValid(); ++mfi) {
                    (mf_avg[lev])[mfi].setVal(static_cast<Real>(npart_in_grid[mfi.index()]));
                }
                dcomp++;
            } else if (fieldname == "part_per_proc"){
                const Vector<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);
                // MultiFab containing number of particles per process
                // (stored as constant for all cells in each grid)
                long n_per_proc = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:n_per_proc)
#endif
                for (MFIter mfi(mf_avg[lev]); mfi.isValid(); ++mfi) {
                    n_per_proc += npart_in_grid[mfi.index()];
                }
                mf_avg[lev].setVal(static_cast<Real>(n_per_proc),dcomp++,1);
            } else if (fieldname == "proc_number"){
                // MultiFab containing the Processor ID
#ifdef _OPENMP
#pragma omp parallel
#endif
                mf_avg[lev].setVal(static_cast<Real>(ParallelDescriptor::MyProc()),dcomp++,1);
            } else if (fieldname == "divB"){
                if (do_nodal) amrex::Abort("TODO: do_nodal && plot divb");
                ComputeDivB(mf_avg[lev], dcomp++,
                            {Bfield_aux[lev][0].get(),
                                    Bfield_aux[lev][1].get(),
                                    Bfield_aux[lev][2].get()},
                            WarpX::CellSize(lev) );
            } else if (fieldname == "divE"){
                if (do_nodal) amrex::Abort("TODO: do_nodal && plot divE");
                const BoxArray& ba = amrex::convert(boxArray(lev),IntVect::TheUnitVector());
                MultiFab divE( ba, DistributionMap(lev), 1, 0 );
#ifdef WARPX_USE_PSATD
                spectral_solver_fp[lev]->ComputeSpectralDivE( Efield_aux[lev], divE );
#else
                ComputeDivE( divE, 0, {Efield_aux[lev][0].get(), Efield_aux[lev][1].get(),
                             Efield_aux[lev][2].get()}, WarpX::CellSize(lev) );
#endif
                AverageAndPackScalarField( mf_avg[lev], divE, dcomp++, ngrow );
            } else {
                amrex::Abort("unknown field in fields_to_plot: " + fieldname);
            }
        }
        if (plot_finepatch)
        {
            AverageAndPackVectorField( mf_avg[lev], Efield_fp[lev], dmap[lev], dcomp, ngrow );
            dcomp += 3;
            AverageAndPackVectorFieldComponents(mf_avg[lev], Efield_fp[lev], dmap[lev], dcomp, ngrow);
            if (lev == 0) {
                AddToVarNames(varnames, "E", "_fp");
                if (n_rz_azimuthal_modes > 1) {
                    AddToVarNamesRZ(varnames, "E", ComponentName("_fp", 0, "real"));
                    for (int mode=1 ; mode < n_rz_azimuthal_modes ; mode++) {
                        AddToVarNamesRZ(varnames, "E", ComponentName("_fp", mode, "real"));
                        AddToVarNamesRZ(varnames, "E", ComponentName("_fp", mode, "imag"));
                    }
                }
            }
            AverageAndPackVectorField( mf_avg[lev], Bfield_fp[lev], dmap[lev], dcomp, ngrow );
            dcomp += 3;
            AverageAndPackVectorFieldComponents(mf_avg[lev], Bfield_fp[lev], dmap[lev], dcomp, ngrow);
            if (lev == 0) {
                AddToVarNames(varnames, "B", "_fp");
                if (n_rz_azimuthal_modes > 1) {
                    AddToVarNamesRZ(varnames, "B", ComponentName("_fp", 0, "real"));
                    for (int mode=1 ; mode < n_rz_azimuthal_modes ; mode++) {
                        AddToVarNamesRZ(varnames, "B", ComponentName("_fp", mode, "real"));
                        AddToVarNamesRZ(varnames, "B", ComponentName("_fp", mode, "imag"));
                    }
                }
            }
        }

        if (plot_crsepatch)
        {
            if (lev == 0)
            {
                mf_avg[lev].setVal(0.0, dcomp, 3, ngrow);
            }
            else
            {
                if (do_nodal) amrex::Abort("TODO: do_nodal && plot_crsepatch");
                std::array<std::unique_ptr<MultiFab>, 3> E = getInterpolatedE(lev);
                AverageAndPackVectorField( mf_avg[lev], E, dmap[lev], dcomp, ngrow );

            }
            if (lev == 0) AddToVarNames(varnames, "E", "_cp");
            dcomp += 3;

            // now the magnetic field
            if (lev == 0)
            {
                mf_avg[lev].setVal(0.0, dcomp, 3, ngrow);
            }
            else
            {
                if (do_nodal) amrex::Abort("TODO: do_nodal && plot_crsepatch");
                std::array<std::unique_ptr<MultiFab>, 3> B = getInterpolatedB(lev);
                AverageAndPackVectorField( mf_avg[lev], B, dmap[lev], dcomp, ngrow );
            }
            if (lev == 0) AddToVarNames(varnames, "B", "_cp");
            dcomp += 3;
        }

        if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            if (costs[0] != nullptr and plot_costs)
            {
                AverageAndPackScalarField( mf_avg[lev], *costs[lev], dcomp, ngrow );
                if(lev==0) varnames.push_back("costs");
                dcomp += 1;
            }
        }
        BL_ASSERT(dcomp == ncomp);
    } // end loop over levels of refinement

};

/** \brief Reduce the size of all the fields in `source_mf`
  * by `coarse_ratio` and store the results in `coarse_mf`.
  * Calculate the corresponding coarse Geometry from `source_geom`
  * and store the results in `coarse_geom` */
void
coarsenCellCenteredFields(
    Vector<MultiFab>& coarse_mf, Vector<Geometry>& coarse_geom,
    const Vector<MultiFab>& source_mf, const Vector<Geometry>& source_geom,
    int coarse_ratio, int finest_level )
{
    // Check that the Vectors to be filled have an initial size of 0
    AMREX_ALWAYS_ASSERT( coarse_mf.size()==0 );
    AMREX_ALWAYS_ASSERT( coarse_geom.size()==0 );

    // Fill the vectors with one element per level
    int ncomp = source_mf[0].nComp();
    for (int lev=0; lev<=finest_level; lev++) {
        AMREX_ALWAYS_ASSERT( source_mf[lev].is_cell_centered() );

        coarse_geom.push_back(amrex::coarsen( source_geom[lev], IntVect(coarse_ratio)));

        BoxArray small_ba = amrex::coarsen(source_mf[lev].boxArray(), coarse_ratio);
        coarse_mf.push_back( MultiFab(small_ba, source_mf[lev].DistributionMap(), ncomp, 0) );
        average_down(source_mf[lev], coarse_mf[lev], 0, ncomp, IntVect(coarse_ratio));
    }
};


/** \brief Write the data from MultiFab `F` into the file `filename`
 *  as a raw field (i.e. no interpolation to cell centers).
 *  Write guard cells if `plot_guards` is True.
 */
void
WriteRawField( const MultiFab& F, const DistributionMapping& dm,
                const std::string& filename,
                const std::string& level_prefix,
                const std::string& field_name,
                const int lev, const bool plot_guards )
{
    std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                            filename, level_prefix, field_name);

    if (plot_guards) {
        // Dump original MultiFab F
        VisMF::Write(F, prefix);
    } else {
        // Copy original MultiFab into one that does not have guard cells
        MultiFab tmpF( F.boxArray(), dm, F.nComp(), 0);
        MultiFab::Copy(tmpF, F, 0, 0, F.nComp(), 0);
        VisMF::Write(tmpF, prefix);
    }

}

/** \brief Write a multifab of the same shape as `F` but filled with 0.
 *  (The shape includes guard cells if `plot_guards` is True.)
 *  This is mainly needed because the yt reader requires all levels of the
 *  coarse/fine patch to be written, but WarpX does not have data for
 *  the coarse patch of level 0 (meaningless).
 */
void
WriteZeroRawField( const MultiFab& F, const DistributionMapping& dm,
                const std::string& filename,
                const std::string& level_prefix,
                const std::string& field_name,
                const int lev, const int ng )
{
    std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                            filename, level_prefix, field_name);

    MultiFab tmpF(F.boxArray(), dm, F.nComp(), ng);
    tmpF.setVal(0.);
    VisMF::Write(tmpF, prefix);
}

/** \brief Write the coarse scalar multifab `F_cp` to the file `filename`
 *  *after* sampling/interpolating its value on the fine grid corresponding
 *  to `F_fp`. This is mainly needed because the yt reader requires the
 *  coarse and fine patch to have the same shape.
 */
void
WriteCoarseScalar( const std::string field_name,
    const std::unique_ptr<MultiFab>& F_cp,
    const std::unique_ptr<MultiFab>& F_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards,
    const int icomp )
{
    int ng = 0;
    if (plot_guards) ng = F_fp->nGrow();

    if (lev == 0) {
        // No coarse field for level 0: instead write a MultiFab
        // filled with 0, with the same number of cells as the _fp field
        WriteZeroRawField( *F_fp, dm, filename, level_prefix, field_name+"_cp", lev, ng );
    } else {
        // Create an alias to the component `icomp` of F_cp
        MultiFab F_comp(*F_cp, amrex::make_alias, icomp, 1);
        // Interpolate coarse data onto fine grid
        const int r_ratio = WarpX::GetInstance().refRatio(lev-1)[0];
        const Real* dx = WarpX::GetInstance().Geom(lev-1).CellSize();
        auto F = getInterpolatedScalar( F_comp, *F_fp, dm, r_ratio, dx, ng );
        // Write interpolated raw data
        WriteRawField( *F, dm, filename, level_prefix, field_name+"_cp", lev, plot_guards );
    }
}

/** \brief Write the coarse vector multifab `F*_cp` to the file `filename`
 *  *after* sampling/interpolating its value on the fine grid corresponding
 *  to `F*_fp`. This is mainly needed because the yt reader requires the
 *  coarse and fine patch to have the same shape.
 */
void
WriteCoarseVector( const std::string field_name,
    const std::unique_ptr<MultiFab>& Fx_cp,
    const std::unique_ptr<MultiFab>& Fy_cp,
    const std::unique_ptr<MultiFab>& Fz_cp,
    const std::unique_ptr<MultiFab>& Fx_fp,
    const std::unique_ptr<MultiFab>& Fy_fp,
    const std::unique_ptr<MultiFab>& Fz_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards )
{
    int ng = 0;
    if (plot_guards) ng = Fx_fp->nGrow();

    if (lev == 0) {
        // No coarse field for level 0: instead write a MultiFab
        // filled with 0, with the same number of cells as the _fp field
        WriteZeroRawField( *Fx_fp, dm, filename, level_prefix, field_name+"x_cp", lev, ng );
        WriteZeroRawField( *Fy_fp, dm, filename, level_prefix, field_name+"y_cp", lev, ng );
        WriteZeroRawField( *Fz_fp, dm, filename, level_prefix, field_name+"z_cp", lev, ng );
    } else {
        // Interpolate coarse data onto fine grid
        const int r_ratio = WarpX::GetInstance().refRatio(lev-1)[0];
        const Real* dx = WarpX::GetInstance().Geom(lev-1).CellSize();
        auto F = getInterpolatedVector( Fx_cp, Fy_cp, Fz_cp, Fx_fp, Fy_fp, Fz_fp,
                                    dm, r_ratio, dx, ng );
        // Write interpolated raw data
        WriteRawField( *F[0], dm, filename, level_prefix, field_name+"x_cp", lev, plot_guards );
        WriteRawField( *F[1], dm, filename, level_prefix, field_name+"y_cp", lev, plot_guards );
        WriteRawField( *F[2], dm, filename, level_prefix, field_name+"z_cp", lev, plot_guards );
    }
}

/** \brief Samples/Interpolates the coarse scalar multifab `F_cp` on the
  * fine grid associated with the fine multifab `F_fp`.
  */
std::unique_ptr<MultiFab>
getInterpolatedScalar(
    const MultiFab& F_cp, const MultiFab& F_fp,
    const DistributionMapping& dm, const int r_ratio,
    const Real* /*dx*/, const int ngrow )
{
    // Prepare the structure that will contain the returned fields
    std::unique_ptr<MultiFab> interpolated_F;
    interpolated_F.reset( new MultiFab(F_fp.boxArray(), dm, 1, ngrow) );
    interpolated_F->setVal(0.);

    // Loop through the boxes and interpolate the values from the _cp data
#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        FArrayBox ffab; // Temporary array ; contains interpolated fields
        for (MFIter mfi(*interpolated_F); mfi.isValid(); ++mfi)
        {
            Box finebx = mfi.fabbox();
            finebx.coarsen(r_ratio).refine(r_ratio); // so that finebx is coarsenable

            const FArrayBox& cfab = (F_cp)[mfi];
            ffab.resize(finebx);

            // - Fully nodal
            if ( F_fp.is_nodal() ){
                IntVect refinement_vector{AMREX_D_DECL(r_ratio, r_ratio, r_ratio)};
                node_bilinear_interp.interp(cfab, 0, ffab, 0, 1,
                        finebx, refinement_vector, {}, {}, {}, 0, 0, RunOn::Cpu);
            } else {
                amrex::Abort("Unknown field staggering.");
            }

            // Add temporary array to the returned structure
            const Box& bx = (*interpolated_F)[mfi].box();
            (*interpolated_F)[mfi].plus(ffab, bx, bx, 0, 0, 1);
        }
    }
    return interpolated_F;
}

/** \brief Samples/Interpolates the coarse vector multifab `F*_cp` on the
  * fine grid associated with the fine multifab `F*_fp`.
  */
std::array<std::unique_ptr<MultiFab>, 3>
getInterpolatedVector(
    const std::unique_ptr<MultiFab>& Fx_cp,
    const std::unique_ptr<MultiFab>& Fy_cp,
    const std::unique_ptr<MultiFab>& Fz_cp,
    const std::unique_ptr<MultiFab>& Fx_fp,
    const std::unique_ptr<MultiFab>& Fy_fp,
    const std::unique_ptr<MultiFab>& Fz_fp,
    const DistributionMapping& dm, const int r_ratio,
    const Real* dx, const int ngrow )
{

    // Prepare the structure that will contain the returned fields
    std::array<std::unique_ptr<MultiFab>, 3> interpolated_F;
    interpolated_F[0].reset( new MultiFab(Fx_fp->boxArray(), dm, 1, ngrow) );
    interpolated_F[1].reset( new MultiFab(Fy_fp->boxArray(), dm, 1, ngrow) );
    interpolated_F[2].reset( new MultiFab(Fz_fp->boxArray(), dm, 1, ngrow) );
    for (int i=0; i<3; i++) interpolated_F[i]->setVal(0.);

    // Loop through the boxes and interpolate the values from the _cp data
    const int use_limiter = 0;
#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,3> ffab; // Temporary array ; contains interpolated fields
        for (MFIter mfi(*interpolated_F[0]); mfi.isValid(); ++mfi)
        {
            Box ccbx = mfi.fabbox();
            ccbx.enclosedCells();
            ccbx.coarsen(r_ratio).refine(r_ratio); // so that ccbx is coarsenable

            const FArrayBox& cxfab = (*Fx_cp)[mfi];
            const FArrayBox& cyfab = (*Fy_cp)[mfi];
            const FArrayBox& czfab = (*Fz_cp)[mfi];
            ffab[0].resize(amrex::convert(ccbx,(*Fx_fp)[mfi].box().type()));
            ffab[1].resize(amrex::convert(ccbx,(*Fy_fp)[mfi].box().type()));
            ffab[2].resize(amrex::convert(ccbx,(*Fz_fp)[mfi].box().type()));

            // - Face centered, in the same way as B on a Yee grid
            if ( (*Fx_fp)[mfi].box().type() == IntVect{AMREX_D_DECL(1,0,0)} ){
#if (AMREX_SPACEDIM == 3)
                amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                             BL_TO_FORTRAN_ANYD(ffab[0]),
                                             BL_TO_FORTRAN_ANYD(ffab[1]),
                                             BL_TO_FORTRAN_ANYD(ffab[2]),
                                             BL_TO_FORTRAN_ANYD(cxfab),
                                             BL_TO_FORTRAN_ANYD(cyfab),
                                             BL_TO_FORTRAN_ANYD(czfab),
                                             dx, &r_ratio, &use_limiter);
#else
                amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                             BL_TO_FORTRAN_ANYD(ffab[0]),
                                             BL_TO_FORTRAN_ANYD(ffab[2]),
                                             BL_TO_FORTRAN_ANYD(cxfab),
                                             BL_TO_FORTRAN_ANYD(czfab),
                                             dx, &r_ratio, &use_limiter);
                amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                       BL_TO_FORTRAN_ANYD(ffab[1]),
                                       BL_TO_FORTRAN_ANYD(cyfab),
                                       &r_ratio, &use_limiter);
#endif
            // - Edge centered, in the same way as E on a Yee grid
            } else if ( (*Fx_fp)[mfi].box().type() == IntVect{AMREX_D_DECL(0,1,1)} ){
#if (AMREX_SPACEDIM == 3)
                amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                    BL_TO_FORTRAN_ANYD(ffab[0]),
                                    BL_TO_FORTRAN_ANYD(ffab[1]),
                                    BL_TO_FORTRAN_ANYD(ffab[2]),
                                    BL_TO_FORTRAN_ANYD(cxfab),
                                    BL_TO_FORTRAN_ANYD(cyfab),
                                    BL_TO_FORTRAN_ANYD(czfab),
                                    &r_ratio, &use_limiter);
#else
                amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                    BL_TO_FORTRAN_ANYD(ffab[0]),
                                    BL_TO_FORTRAN_ANYD(ffab[2]),
                                    BL_TO_FORTRAN_ANYD(cxfab),
                                    BL_TO_FORTRAN_ANYD(czfab),
                                    &r_ratio,&use_limiter);
                amrex_interp_nd_efield(ccbx.loVect(), ccbx.hiVect(),
                                       BL_TO_FORTRAN_ANYD(ffab[1]),
                                       BL_TO_FORTRAN_ANYD(cyfab),
                                       &r_ratio);
#endif
            } else {
                amrex::Abort("Unknown field staggering.");
            }

            // Add temporary array to the returned structure
            for (int i = 0; i < 3; ++i) {
                const Box& bx = (*interpolated_F[i])[mfi].box();
                (*interpolated_F[i])[mfi].plus(ffab[i], bx, bx, 0, 0, 1);
            }
        }
    }
    return interpolated_F;
}

std::array<std::unique_ptr<MultiFab>, 3> WarpX::getInterpolatedE(int lev) const
{

    const int ngrow = 0;
    const DistributionMapping& dm = DistributionMap(lev);
    const Real* dx = Geom(lev-1).CellSize();
    const int r_ratio = refRatio(lev-1)[0];

    return getInterpolatedVector(
        Efield_cp[lev][0], Efield_cp[lev][1], Efield_cp[lev][2],
        Efield_fp[lev][0], Efield_fp[lev][1], Efield_fp[lev][2],
        dm, r_ratio, dx, ngrow );
}

std::array<std::unique_ptr<MultiFab>, 3> WarpX::getInterpolatedB(int lev) const
{
    const int ngrow = 0;
    const DistributionMapping& dm = DistributionMap(lev);
    const Real* dx = Geom(lev-1).CellSize();
    const int r_ratio = refRatio(lev-1)[0];

    return getInterpolatedVector(
        Bfield_cp[lev][0], Bfield_cp[lev][1], Bfield_cp[lev][2],
        Bfield_fp[lev][0], Bfield_fp[lev][1], Bfield_fp[lev][2],
        dm, r_ratio, dx, ngrow );
}
