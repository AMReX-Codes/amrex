#ifndef WARPX_WRITER_H
#define WARPX_WRITER_H

#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_AmrParticles.H>

#include "pseudo_warpx.H"
#include "AMReX_PlotFileUtilOPENPMD.H"

class AMReX_warpxWriter: public amrex::openpmd_api::AMReX_openPMDWriter
{
public:
  AMReX_warpxWriter(std::vector<bool> fieldPMLdirections,
                    openPMD::IterationEncoding ie =  openPMD::IterationEncoding::variableBased,
                    std::string filetype = "default",
                    std::string openpmdOptions = "{}")
    :AMReX_openPMDWriter("{}", ie, std::move(filetype), std::move(openpmdOptions)),
     m_fieldPMLdirections(std::move(fieldPMLdirections))
  {}

  //
  //  fields related funcitions
  //
  void SetupMeshComp (openPMD::Mesh& mesh,
                      const amrex::Geometry& full_geom,
                      amrex::MultiFab const& mf,
                      const amrex::openpmd_api::AMReX_VarNameParser& varName) const override;

  void SetupFields (openPMD::Container< openPMD::Mesh >& meshes,
                    amrex::Geometry& full_geom) const override;

  [[nodiscard]] inline std::vector< std::string > getParticlePositionComponentLabels() const override;

  virtual void SetConstantMassCharge(openPMD::ParticleSpecies& currSpecies,
                                     const unsigned long long& np,
                                     amrex::ParticleReal charge,
                                     amrex::ParticleReal mass) const;

  void SetParticleSpecieAttributes(openPMD::ParticleSpecies& currSpecies) const override;

private:
  std::vector<bool> m_fieldPMLdirections;

  [[nodiscard]] static std::vector< std::string > getFieldAxisLabels ( const amrex::openpmd_api::AMReX_VarNameParser& varName );
};


void AMReX_warpxWriter::SetupMeshComp (openPMD::Mesh& mesh,
                                       const amrex::Geometry& full_geom,
                                       amrex::MultiFab const& mf,
                                       const amrex::openpmd_api::AMReX_VarNameParser& varName) const
{
  BL_PROFILE("SetupMeshComp(WarpX)");

  auto mesh_comp = mesh[varName.m_CompName];
  amrex::Box const & global_box = full_geom.Domain();
  auto global_size = amrex::openpmd_api::helper::getReversedVec(global_box.size() );

  // - Grid spacing
  std::vector<double> const grid_spacing = amrex::openpmd_api::helper::getReversedVec(full_geom.CellSize());
  mesh.setGridSpacing(grid_spacing);

  // - Global offset
  std::vector<double> const global_offset = amrex::openpmd_api::helper::getReversedVec(full_geom.ProbLo());
  mesh.setGridGlobalOffset(global_offset);

  // - AxisLabels
  std::vector<std::string> axis_labels = getFieldAxisLabels(varName);
  mesh.setAxisLabels(axis_labels);

#if defined(WARPX_DIM_RZ)
  auto & warpx = WarpX::GetInstance();
  if (curr.m_ThetaMode) {
    global_size.emplace(global_size.begin(), warpx.ncomps);
  }
#endif

  // Prepare the type of dataset that will be written
  openPMD::Datatype const datatype = openPMD::determineDatatype<amrex::Real>();
  auto const dataset = openPMD::Dataset(datatype, global_size);
  mesh.setDataOrder(openPMD::Mesh::DataOrder::C);

  if (varName.m_ThetaMode) {
    mesh.setGeometry("thetaMode");
    mesh.setGeometryParameters("m=" + std::to_string(WarpX::n_rz_azimuthal_modes) + ";imag=+");
  }

  mesh.setAttribute("fieldSmoothing", "none");
  mesh_comp.resetDataset(dataset);

  amrex::openpmd_api::helper::setOpenPMDUnit( mesh, varName.m_FieldName );

  auto relative_cell_pos = amrex::openpmd_api::helper::getRelativeCellPosition(mf);     // AMReX Fortran index order
  std::reverse( relative_cell_pos.begin(), relative_cell_pos.end() ); // now in C order
  mesh_comp.setPosition( relative_cell_pos );
}

void AMReX_warpxWriter::SetupFields (openPMD::Container< openPMD::Mesh >& meshes,
                                      amrex::Geometry& full_geom) const
{
  // meta data for ED-PIC extension
  auto const period = full_geom.periodicity();

  std::vector<std::string> fieldBoundary(6, "reflecting");
  std::vector<std::string> particleBoundary(6, "absorbing");
  fieldBoundary.resize(AMREX_SPACEDIM * 2);
  particleBoundary.resize(AMREX_SPACEDIM * 2);

#if AMREX_SPACEDIM != 3
  fieldBoundary.resize(4);
  particleBoundary.resize(4);
#endif
  if (m_fieldPMLdirections.size() >= fieldBoundary.size())
  {
      for (int i = 0; i < int(fieldBoundary.size() / 2); ++i) {
          if (m_fieldPMLdirections.at(i)) {
              fieldBoundary.at(i) = "open";
          }
      }
  }

  for (int i = 0; i < int(fieldBoundary.size() / 2); ++i) {
    if (period.isPeriodic(i)) {
      fieldBoundary.at(2 * i) = "periodic";
      fieldBoundary.at(2 * i + 1) = "periodic";
      particleBoundary.at(2 * i) = "periodic";
      particleBoundary.at(2 * i + 1) = "periodic";
    }
  }

  meshes.setAttribute("fieldSolver", []() {
    switch (WarpX::electromagnetic_solver_id) {
    case WarpX::ElectromagneticSolverAlgo::Yee :
        { return "Yee"; }
    case WarpX::ElectromagneticSolverAlgo::CKC :
        { return "CK"; }
    case WarpX::ElectromagneticSolverAlgo::PSATD :
        { return "PSATD"; }
    default:
        { return "other"; }
    }
  }())
    ;

  meshes.setAttribute("fieldBoundary", fieldBoundary);
  meshes.setAttribute("particleBoundary", particleBoundary);

  meshes.setAttribute("currentSmoothing", []() {
    if (WarpX::use_filter) { return "Binomial"; }
    else { return "none"; }
  }());

  if (WarpX::use_filter) {
    meshes.setAttribute("currentSmoothingParameters", []() {
      std::stringstream ss;
      ss << "period=1;compensator=false";
#if (AMREX_SPACEDIM >= 2)
      ss << ";numPasses_x=" << WarpX::filter_npass_each_dir[0];
#endif
#if (AMREX_SPACEDIM == 3)
      ss << ";numPasses_y=" << WarpX::filter_npass_each_dir[1];
      ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[2];
#endif
#if defined(WARPX_DIM_3D)
      ss << ";numPasses_y=" << WarpX::filter_npass_each_dir[1];
      ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[2];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
      ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[1];
#elif defined(WARPX_DIM_1D_Z)
      ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[0];
#endif
      std::string currentSmoothingParameters = ss.str();
      return currentSmoothingParameters;
    }());
  }

  meshes.setAttribute("chargeCorrection", []() {
    if (WarpX::do_dive_cleaning) { return "hyperbolic"; } // TODO or "spectral" or something? double-check
    else { return "none"; }
  }());

  if (WarpX::do_dive_cleaning) {
    meshes.setAttribute("chargeCorrectionParameters", "period=1");
  }
}



std::vector< std::string > AMReX_warpxWriter::getFieldAxisLabels ( const amrex::openpmd_api::AMReX_VarNameParser& /*varName*/ )
  {
    using vs = std::vector< std::string >;
    // Fortran order of the index labels for the AMReX FArrayBox
#if defined(WARPX_DIM_1D_Z)
    vs const axisLabels{"z"};  // z varies fastest in memory
#elif defined(WARPX_DIM_XZ)
    vs const axisLabels{"x", "z"};  // x varies fastest in memory
#elif defined(WARPX_DIM_RZ)
    // when we write individual modes of a field (default)
    vs const circAxisLabels{"r", "z"};  // r varies fastest in memory
    // if we just write reconstructed 2D fields at theta=0
    vs const cartAxisLabels{"x", "z"};  // x varies fastest in memory
    vs const axisLabels = varName.m_ThetaMode ? circAxisLabels : cartAxisLabels;
#elif defined(WARPX_DIM_3D)
    vs const axisLabels{"x", "y", "z"};  // x varies fastest in memory
#else
    //error Unknown WarpX dimensionality.
    vs const axisLabels{"w", "k", "i"};  // default as this is just a test
#endif
    // revert to C order (fastest varying index last)
    return {axisLabels.rbegin(), axisLabels.rend()};
  }

  //
  // ptl related warpx functions
  //
inline std::vector< std::string >
AMReX_warpxWriter::getParticlePositionComponentLabels() const
  {
    using vs = std::vector< std::string >;
#if defined(WARPX_DIM_1D_Z)
    vs positionComponents{"z"};
#elif defined(WARPX_DIM_XZ)
    vs positionComponents{"x", "z"};
#elif defined(WARPX_DIM_RZ)
    // note: although we internally store particle positions
    //       for AMReX in r,z and a theta attribute, we
    //       actually need them for algorithms (e.g. push)
    //       and I/O in Cartesian.
    //       Other attributes like momentum are consequently
    //       stored in x,y,z internally.
    vs positionComponents{"x", "y", "z"};
#elif (WARPX_DIM_3D)
    vs positionComponents{"x", "y", "z"};
#else
    //#   error Unknown WarpX dimensionality.
    // default:
    vs positionComponents{"kk", "jj", "ii"};
#endif
    return positionComponents;
  }

  // TODO should it be here?
void AMReX_warpxWriter::SetConstantMassCharge(openPMD::ParticleSpecies& currSpecies,
                                              const unsigned long long& np,
                                              amrex::ParticleReal const charge,
                                              amrex::ParticleReal const mass) const
  {
    auto realType = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});
    auto const * const scalar = openPMD::RecordComponent::SCALAR;

    currSpecies["charge"][scalar].resetDataset( realType );
    currSpecies["mass"][scalar].resetDataset( realType );

    currSpecies["charge"][scalar].makeConstant( charge );
    currSpecies["mass"][scalar].makeConstant( mass );

    currSpecies["charge"].setUnitDimension( amrex::openpmd_api::helper::getUnitDimension("charge") );
    currSpecies["mass"].setUnitDimension( amrex::openpmd_api::helper::getUnitDimension("mass") );

    currSpecies["charge"].setAttribute( "macroWeighted", 0U );
    currSpecies["charge"].setAttribute( "weightingPower", 1.0 );
    currSpecies["mass"].setAttribute( "macroWeighted", 0U );
    currSpecies["mass"].setAttribute( "weightingPower", 1.0 );
  }



  // should be called after
  //    - setupPos  so records "position" and "id" are defined
  //    - setConstants so record "positonOffset" is defined
void AMReX_warpxWriter::SetParticleSpecieAttributes(openPMD::ParticleSpecies& currSpecies) const
  {

    // ED-PIC extension
    currSpecies["id"].setAttribute( "macroWeighted", 0U );
    currSpecies["id"].setAttribute( "weightingPower", 0.0 );
    currSpecies["position"].setAttribute( "macroWeighted", 0U );
    currSpecies["position"].setAttribute( "weightingPower", 0.0 );

    currSpecies["positionOffset"].setAttribute( "macroWeighted", 0U );
    currSpecies["positionOffset"].setAttribute( "weightingPower", 0.0 );

      // meta data for ED-PIC extension
      currSpecies.setAttribute( "particleShape", double( WarpX::noz ) );

      currSpecies.setAttribute( "particleShapes", [](){
        return std::vector< double >{
          double(WarpX::nox),
#if AMREX_SPACEDIM==3
          double(WarpX::noy),
#endif
          double(WarpX::noz)
        };
      }() );

      currSpecies.setAttribute( "particlePush", [](){
        switch( WarpX::particle_pusher_algo ) {
        case WarpX::ParticlePusherAlgo::Boris : { return "Boris"; }
        case WarpX::ParticlePusherAlgo::Vay : { return "Vay"; }
        case WarpX::ParticlePusherAlgo::HigueraCary : { return "HigueraCary"; }
        default: { return "other"; }
        }
      }() );
      currSpecies.setAttribute( "particleInterpolation", [](){
        switch( WarpX::field_gathering_algo ) {
        case WarpX::GatheringAlgo::EnergyConserving : { return "energyConserving"; }
        case WarpX::GatheringAlgo::MomentumConserving : { return "momentumConserving"; }
        default: { return "other"; }
        }
      }() );
      currSpecies.setAttribute( "particleSmoothing", "none" );
      currSpecies.setAttribute( "currentDeposition", [](){
        switch( WarpX::current_deposition_algo ) {
        case WarpX::CurrentDepositionAlgo::Esirkepov : { return "Esirkepov"; }
        default: { return "directMorseNielson"; }
        }
      }() );
  }


#endif
