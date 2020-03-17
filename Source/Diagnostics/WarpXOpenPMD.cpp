/* Copyright 2019-2020 Axel Huebl, Junmin Gu
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXOpenPMD.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "FieldIO.H"  // for getReversedVec

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <tuple>
#include <utility>
#include <iostream>


namespace detail
{

    /** Convert AMReX AoS .id and .cpu to a globally unique particle ID
     */
    union GlobalID {
        struct { int id; int cpu; }; //! MPI-rank local ID (and rank/cpu)
        uint64_t global_id; //! global ID that is unique in the whole simulation
    };
    static_assert(sizeof(int) * 2u <= sizeof(uint64_t), "int size might cause collisions in global IDs");

#ifdef WARPX_USE_OPENPMD
    /** Unclutter a real_names to openPMD record
     *
     * @param fullName name as in real_names variable
     * @return pair of openPMD record and component name
     */
    inline std::pair< std::string, std::string >
    name2openPMD( std::string const& fullName )
    {
        std::string record_name = fullName;
        std::string component_name = openPMD::RecordComponent::SCALAR;
        std::size_t startComp = fullName.find_last_of("_");

        if( startComp != std::string::npos ) {  // non-scalar
            record_name = fullName.substr(0, startComp);
            component_name = fullName.substr(startComp + 1u);
        }
        return make_pair(record_name, component_name);
    }

    /** Get the openPMD physical dimensionality of a record
     *
     * @param record_name name of the openPMD record
     * @return map with base quantities and power scaling
     */
    std::map< openPMD::UnitDimension, double >
    getUnitDimension( std::string const & record_name )
    {

        if( record_name == "position" ) return {
            {openPMD::UnitDimension::L,  1.}
        };
        else if( record_name == "positionOffset" ) return {
            {openPMD::UnitDimension::L,  1.}
        };
        else if( record_name == "momentum" ) return {
            {openPMD::UnitDimension::L,  1.},
            {openPMD::UnitDimension::M,  1.},
            {openPMD::UnitDimension::T, -1.}
        };
        else if( record_name == "charge" ) return {
            {openPMD::UnitDimension::T,  1.},
            {openPMD::UnitDimension::I,  1.}
        };
        else if( record_name == "mass" ) return {
            {openPMD::UnitDimension::M,  1.}
        };
        else if( record_name == "E" ) return {
            {openPMD::UnitDimension::L,  1.},
            {openPMD::UnitDimension::M,  1.},
            {openPMD::UnitDimension::T, -3.},
            {openPMD::UnitDimension::I, -1.},
        };
        else if( record_name == "B" ) return {
            {openPMD::UnitDimension::M,  1.},
            {openPMD::UnitDimension::I, -1.},
            {openPMD::UnitDimension::T, -2.}
        };
        else return {};
    }
#endif // WARPX_USE_OPENPMD
}

#ifdef WARPX_USE_OPENPMD
WarpXOpenPMDPlot::WarpXOpenPMDPlot(bool oneFilePerTS,
    std::string openPMDFileType, std::vector<bool> fieldPMLdirections)
  :m_Series(nullptr),
   m_OneFilePerTS(oneFilePerTS),
   m_OpenPMDFileType(std::move(openPMDFileType)),
   m_fieldPMLdirections(std::move(fieldPMLdirections))
{
  // pick first available backend if default is chosen
  if( m_OpenPMDFileType == "default" )
#if openPMD_HAVE_ADIOS2==1
    m_OpenPMDFileType = "bp";
#elif openPMD_HAVE_ADIOS1==1
    m_OpenPMDFileType = "bp";
#elif openPMD_HAVE_HDF5==1
    m_OpenPMDFileType = "h5";
#else
    m_OpenPMDFileType = "json";
#endif
}

WarpXOpenPMDPlot::~WarpXOpenPMDPlot()
{
  if( m_Series )
  {
    m_Series->flush();
    m_Series.reset( nullptr );
  }
}


//
//
//
void WarpXOpenPMDPlot::GetFileName(std::string& filename)
{
  std::string dir = "diags/";

  filename = dir;
  filename.append(m_OpenPMDFileType).append("/simData");
  //
  // OpenPMD supports timestepped names
  //
  if (m_OneFilePerTS)
      filename = filename.append("_%07T");
  filename.append(".").append(m_OpenPMDFileType);
}


void WarpXOpenPMDPlot::SetStep(int ts)
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ts >= 0 , "openPMD iterations are unsigned");

  if (m_CurrentStep >= ts) {
      // note m_Series is reset in Init(), so using m_Series->iterations.contains(ts) is only able to check the
      // last written step in m_Series's life time, but not other earlier written steps by other m_Series
      std::string warnMsg = " Warning from openPMD writer: Already written iteration:"+std::to_string(ts);
      std::cout<<warnMsg<<std::endl;
      amrex::Warning(warnMsg);
  }

    m_CurrentStep =  ts;
    Init(openPMD::AccessType::CREATE);

}

void
WarpXOpenPMDPlot::Init(openPMD::AccessType accessType)
{
    // either for the next ts file,
    // or init a single file for all ts
    std::string filename;
    GetFileName(filename);

    // close a previously open series before creating a new one
    // see ADIOS1 limitation: https://github.com/openPMD/openPMD-api/pull/686
    m_Series = nullptr;

    if( amrex::ParallelDescriptor::NProcs() > 1 )
    {
        m_Series = std::make_unique<openPMD::Series>(
            filename, accessType,
            amrex::ParallelDescriptor::Communicator()
        );
        m_MPISize = amrex::ParallelDescriptor::NProcs();
        m_MPIRank = amrex::ParallelDescriptor::MyProc();
    }
    else
    {
        m_Series = std::make_unique<openPMD::Series>(filename, accessType);
        m_MPISize = 1;
        m_MPIRank = 1;
    }

    // input file / simulation setup author
    if( WarpX::authors.size() > 0u )
        m_Series->setAuthor( WarpX::authors );
    // more natural naming for PIC
    m_Series->setMeshesPath( "fields" );
    // conform to ED-PIC extension of openPMD
    uint32_t const openPMD_ED_PIC = 1u;
    m_Series->setOpenPMDextension( openPMD_ED_PIC );
    // meta info
    m_Series->setSoftware( "WarpX", WarpX::Version() );
}


void
WarpXOpenPMDPlot::WriteOpenPMDParticles(const std::unique_ptr<MultiParticleContainer>& mpc)
{
  WARPX_PROFILE("WarpXOpenPMDPlot::WriteOpenPMDParticles()");
  std::vector<std::string> species_names = mpc->GetSpeciesNames();

  for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
    auto& pc  = mpc->GetUniqueContainer(i);
    if (pc->plot_species) {

      // names of amrex::Real and int particle attributes in SoA data
      amrex::Vector<std::string> real_names;
      amrex::Vector<std::string> int_names;
      amrex::Vector<int> int_flags;

      // see openPMD ED-PIC extension for namings
      // note: an underscore separates the record name from its component
      //       for non-scalar records
      real_names.push_back("weighting");

      real_names.push_back("momentum_x");
      real_names.push_back("momentum_y");
      real_names.push_back("momentum_z");

      real_names.push_back("E_x");
      real_names.push_back("E_y");
      real_names.push_back("E_z");

      real_names.push_back("B_x");
      real_names.push_back("B_y");
      real_names.push_back("B_z");

#ifdef WARPX_DIM_RZ
      real_names.push_back("theta");
#endif
      if(pc->do_field_ionization){
         int_names.push_back("ionization_level");
         // int_flags specifies, for each integer attribs, whether it is
         // dumped as particle record in a plotfile. So far, ionization_level is the only
         // integer attribs, and it is automatically dumped as particle record
         // when ionization is on.
         int_flags.resize(1, 1);
      }

      // Convert momentum to SI
      pc->ConvertUnits(ConvertDirection::WarpX_to_SI);
      // real_names contains a list of all real particle attributes.
      // pc->plot_flags is 1 or 0, whether quantity is dumped or not.

      {
        //
        SavePlotFile(pc,
           species_names[i],
           m_CurrentStep,
           pc->plot_flags, // this is protected and accessible by MultiParticleContainer.
           // so kept as is
           int_flags,
           real_names, int_names);
      }

      // Convert momentum back to WarpX units
      pc->ConvertUnits(ConvertDirection::SI_to_WarpX);
    }
  }
}



void
WarpXOpenPMDPlot::SavePlotFile (const std::unique_ptr<WarpXParticleContainer>& pc,
                    const std::string& name,
                    int iteration,
                    const amrex::Vector<int>& write_real_comp,
                    const amrex::Vector<int>& write_int_comp,
                    const amrex::Vector<std::string>& real_comp_names,
                    const amrex::Vector<std::string>&  int_comp_names) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Series != nullptr, "openPMD series must be initialized");

  WarpXParticleCounter counter(pc);

  openPMD::Iteration currIteration = m_Series->iterations[iteration];
  openPMD::ParticleSpecies currSpecies = currIteration.particles[name];
  // meta data for ED-PIC extension
  currSpecies.setAttribute( "particleShape", double( WarpX::noz ) );
  // TODO allow this per direction in the openPMD standard, ED-PIC extension?
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
          case ParticlePusherAlgo::Boris : return "Boris";
          case ParticlePusherAlgo::Vay : return "Vay";
          case ParticlePusherAlgo::HigueraCary : return "HigueraCary";
          default: return "other";
      }
  }() );
  currSpecies.setAttribute( "particleInterpolation", [](){
      switch( WarpX::field_gathering_algo ) {
          case GatheringAlgo::EnergyConserving : return "energyConserving";
          case GatheringAlgo::MomentumConserving : return "momentumConserving";
          default: return "other";
      }
  }() );
  currSpecies.setAttribute( "particleSmoothing", "none" );
  currSpecies.setAttribute( "currentDeposition", [](){
      switch( WarpX::current_deposition_algo ) {
          case CurrentDepositionAlgo::Esirkepov : return "Esirkepov";
          default: return "directMorseNielson";
      }
  }() );

  //
  // define positions & offsets
  //
  SetupPos(pc, currSpecies, counter.GetTotalNumParticles());
  SetupRealProperties(currSpecies, write_real_comp, real_comp_names, counter.GetTotalNumParticles());

  // open files from all processors, in case some will not contribute below
  m_Series->flush();

  for (auto currentLevel = 0; currentLevel <= pc->finestLevel(); currentLevel++)
    {
      uint64_t offset = static_cast<uint64_t>( counter.m_ParticleOffsetAtRank[currentLevel] );

      for (WarpXParIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
         auto const numParticleOnTile = pti.numParticles();
         uint64_t const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

         // get position and particle ID from aos
         // note: this implementation iterates the AoS 4x...
         // if we flush late as we do now, we can also copy out the data in one go
         const auto& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
         {
           // Save positions
           std::vector<std::string> axisNames={"x", "y", "z"};
           for (auto currDim = 0; currDim < AMREX_SPACEDIM; currDim++) {
                std::shared_ptr< amrex::ParticleReal > curr(
                    new amrex::ParticleReal[numParticleOnTile],
                    [](amrex::ParticleReal const *p){ delete[] p; }
                );
                for (auto i=0; i<numParticleOnTile; i++) {
                     curr.get()[i] = aos[i].m_rdata.pos[currDim];
                }
                currSpecies["position"][axisNames[currDim]].storeChunk(curr, {offset}, {numParticleOnTile64});
           }

           // save particle ID after converting it to a globally unique ID
           std::shared_ptr< uint64_t > ids(
               new uint64_t[numParticleOnTile],
               [](uint64_t const *p){ delete[] p; }
           );
           for (auto i=0; i<numParticleOnTile; i++) {
               detail::GlobalID const nextID = { aos[i].m_idata.id, aos[i].m_idata.cpu };
               ids.get()[i] = nextID.global_id;
           }
           auto const scalar = openPMD::RecordComponent::SCALAR;
           currSpecies["id"][scalar].storeChunk(ids, {offset}, {numParticleOnTile64});
        }
         //  save "extra" particle properties in AoS and SoA
         SaveRealProperty(pti,
             currSpecies,
             offset,
             write_real_comp, real_comp_names);

         offset += numParticleOnTile64;
      }
    }
    m_Series->flush();
}

void
WarpXOpenPMDPlot::SetupRealProperties(openPMD::ParticleSpecies& currSpecies,
                      const amrex::Vector<int>& write_real_comp,
                      const amrex::Vector<std::string>& real_comp_names,
                      unsigned long long np) const
{
  auto particlesLineup = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});

  //
  // the beam/input3d showed write_real_comp.size() = 16 while only 10 real comp names
  // so using the min to be safe.
  //
  auto counter = std::min(write_real_comp.size(), real_comp_names.size());
  for (int i = 0; i < counter; ++i)
    if (write_real_comp[i]) {
      // handle scalar and non-scalar records by name
      std::string record_name, component_name;
      std::tie(record_name, component_name) = detail::name2openPMD(real_comp_names[i]);

      auto particleVarComp = currSpecies[record_name][component_name];
      particleVarComp.resetDataset(particlesLineup);
    }

  std::set< std::string > addedRecords; // add meta-data per record only once
  for (auto idx=0; idx<m_NumSoARealAttributes; idx++) {
    auto ii = m_NumAoSRealAttributes + idx;
    if (write_real_comp[ii]) {
      // handle scalar and non-scalar records by name
      std::string record_name, component_name;
      std::tie(record_name, component_name) = detail::name2openPMD(real_comp_names[ii]);
      auto currRecord = currSpecies[record_name];

      // meta data for ED-PIC extension
      bool newRecord = false;
      std::tie(std::ignore, newRecord) = addedRecords.insert(record_name);
      if( newRecord ) {
        currRecord.setUnitDimension( detail::getUnitDimension(record_name) );
        currRecord.setAttribute( "macroWeighted", 0u );
        if( record_name == "momentum" )
            currRecord.setAttribute( "weightingPower", 1.0 );
        else
            currRecord.setAttribute( "weightingPower", 0.0 );
      }
    }
  }
}

void
WarpXOpenPMDPlot::SaveRealProperty(WarpXParIter& pti,
                       openPMD::ParticleSpecies& currSpecies,
                       unsigned long long const offset,
                       amrex::Vector<int> const& write_real_comp,
                       amrex::Vector<std::string> const& real_comp_names) const

{
  int numOutputReal = 0;
  int const totalRealAttrs = m_NumAoSRealAttributes + m_NumSoARealAttributes;

  for( int i = 0; i < totalRealAttrs; ++i )
    if( write_real_comp[i] )
      ++numOutputReal;

  auto const numParticleOnTile = pti.numParticles();
  uint64_t const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );
  auto const& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
  auto const& soa = pti.GetStructOfArrays();

  // properties are saved separately
  {
    for( auto idx=0; idx<m_NumAoSRealAttributes; idx++ ) {
      if( write_real_comp[idx] ) {
          // handle scalar and non-scalar records by name
          std::string record_name, component_name;
          std::tie(record_name, component_name) = detail::name2openPMD(real_comp_names[idx]);
          auto currRecord = currSpecies[record_name];
          auto currRecordComp = currRecord[component_name];

          std::shared_ptr< amrex::ParticleReal > d(
              new amrex::ParticleReal[numParticleOnTile],
              [](amrex::ParticleReal const *p){ delete[] p; }
          );

          for( auto kk=0; kk<numParticleOnTile; kk++ )
               d.get()[kk] = aos[kk].m_rdata.arr[AMREX_SPACEDIM+idx];

          currRecordComp.storeChunk(d,
               {offset}, {numParticleOnTile64});
      }
    }
  }

  {
    for (auto idx=0; idx<m_NumSoARealAttributes; idx++) {
      auto ii = m_NumAoSRealAttributes + idx;
      if (write_real_comp[ii]) {
          // handle scalar and non-scalar records by name
          std::string record_name, component_name;
          std::tie(record_name, component_name) = detail::name2openPMD(real_comp_names[ii]);
          auto& currRecord = currSpecies[record_name];
          auto& currRecordComp = currRecord[component_name];

          currRecordComp.storeChunk(openPMD::shareRaw(soa.GetRealData(idx)),
              {offset}, {numParticleOnTile64});
      }
    }
  }
}



void
WarpXOpenPMDPlot::SetupPos(const std::unique_ptr<WarpXParticleContainer>& pc,
    openPMD::ParticleSpecies& currSpecies,
    const unsigned long long& np) const
{
  auto const realType = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});
  auto const idType = openPMD::Dataset(openPMD::determineDatatype< uint64_t >(), {np});

  for( auto const& comp : {"x", "y", "z"} ) {
      currSpecies["positionOffset"][comp].resetDataset( realType );
      currSpecies["positionOffset"][comp].makeConstant( 0. );
      currSpecies["position"][comp].resetDataset( realType );
  }

  auto const scalar = openPMD::RecordComponent::SCALAR;
  currSpecies["id"][scalar].resetDataset( idType );
  currSpecies["charge"][scalar].resetDataset( realType );
  currSpecies["charge"][scalar].makeConstant( pc->getCharge() );
  currSpecies["mass"][scalar].resetDataset( realType );
  currSpecies["mass"][scalar].makeConstant( pc->getMass() );

  // meta data
  currSpecies["position"].setUnitDimension( detail::getUnitDimension("position") );
  currSpecies["positionOffset"].setUnitDimension( detail::getUnitDimension("positionOffset") );
  currSpecies["charge"].setUnitDimension( detail::getUnitDimension("charge") );
  currSpecies["mass"].setUnitDimension( detail::getUnitDimension("mass") );

  // meta data for ED-PIC extension
  currSpecies["position"].setAttribute( "macroWeighted", 0u );
  currSpecies["position"].setAttribute( "weightingPower", 0.0 );
  currSpecies["positionOffset"].setAttribute( "macroWeighted", 0u );
  currSpecies["positionOffset"].setAttribute( "weightingPower", 0.0 );
  currSpecies["id"].setAttribute( "macroWeighted", 0u );
  currSpecies["id"].setAttribute( "weightingPower", 0.0 );
  currSpecies["charge"].setAttribute( "macroWeighted", 0u );
  currSpecies["charge"].setAttribute( "weightingPower", 1.0 );
  currSpecies["mass"].setAttribute( "macroWeighted", 0u );
  currSpecies["mass"].setAttribute( "weightingPower", 1.0 );
}


//
// this is originally copied from FieldIO.cpp
//
void
WarpXOpenPMDPlot::WriteOpenPMDFields( //const std::string& filename,
                      const std::vector<std::string>& varnames,
                      const amrex::MultiFab& mf,
                      const amrex::Geometry& geom,
                      const int iteration,
                      const double time ) const
{
  //This is AMReX's tiny profiler. Possibly will apply it later
  WARPX_PROFILE("WarpXOpenPMDPlot::WriteOpenPMDFields()");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Series != nullptr, "openPMD series must be initialized");

  int const ncomp = mf.nComp();

  // Create a few vectors that store info on the global domain
  // Swap the indices for each of them, since AMReX data is Fortran order
  // and since the openPMD API assumes contiguous C order
  // - Size of the box, in integer number of cells
  amrex::Box const & global_box = geom.Domain();
  auto const global_size = getReversedVec(global_box.size());
  // - Grid spacing
  std::vector<double> const grid_spacing = getReversedVec(geom.CellSize());
  // - Global offset
  std::vector<double> const global_offset = getReversedVec(geom.ProbLo());
  // - AxisLabels
#if AMREX_SPACEDIM==3
  std::vector<std::string> const axis_labels{"x", "y", "z"};
#else
  std::vector<std::string> const axis_labels{"x", "z"};
#endif

  // Prepare the type of dataset that will be written
  openPMD::Datatype const datatype = openPMD::determineDatatype<amrex::Real>();
  auto const dataset = openPMD::Dataset(datatype, global_size);

  // meta data
  auto series_iteration = m_Series->iterations[iteration];
  series_iteration.setTime( time );

  // meta data for ED-PIC extension
  auto const period = geom.periodicity(); // TODO double-check: is this the proper global bound or of some level?
  std::vector< std::string > fieldBoundary( 6, "reflecting" );
  std::vector< std::string > particleBoundary( 6, "absorbing" );
#if AMREX_SPACEDIM!=3
    fieldBoundary.resize(4);
    particleBoundary.resize(4);
#endif

  for( auto i = 0u; i < fieldBoundary.size() / 2u; ++i )
      if( m_fieldPMLdirections.at( i ) )
          fieldBoundary.at( i ) = "open";

  for( auto i = 0u; i < fieldBoundary.size() / 2u; ++i )
      if( period.isPeriodic( i ) ) {
          fieldBoundary.at(2u*i     ) = "periodic";
          fieldBoundary.at(2u*i + 1u) = "periodic";
          particleBoundary.at(2u*i     ) = "periodic";
          particleBoundary.at(2u*i + 1u) = "periodic";
      }

  auto meshes = series_iteration.meshes;
  meshes.setAttribute( "fieldSolver", [](){
#ifdef WARPX_USE_PSATD
      return "PSATD"; // TODO double-check if WARPX_USE_PSATD_HYBRID is covered
#else
      switch( WarpX::particle_pusher_algo ) {
          case MaxwellSolverAlgo::Yee : return "Yee";
          case MaxwellSolverAlgo::CKC : return "CK";
          default: return "other";
      }
#endif
  }() );
  meshes.setAttribute( "fieldBoundary", fieldBoundary );
  meshes.setAttribute( "particleBoundary", particleBoundary );
  meshes.setAttribute( "currentSmoothing", [](){
      if( WarpX::use_filter ) return "Binomial";
          else return "none";
  }() );
    if( WarpX::use_filter )
        meshes.setAttribute( "currentSmoothingParameters", [](){
            std::stringstream ss;
            ss << "period=1;compensator=false";
            ss << ";numPasses_x=" << WarpX::filter_npass_each_dir[0];
#if (AMREX_SPACEDIM == 3)
            ss << ";numPasses_y=" << WarpX::filter_npass_each_dir[1];
            ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[2];
#else
            ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[1];
#endif
            std::string currentSmoothingParameters = ss.str();
            return std::move(currentSmoothingParameters);
        }() );
  meshes.setAttribute("chargeCorrection", [](){
      if( WarpX::do_dive_cleaning ) return "hyperbolic"; // TODO or "spectral" or something? double-check
      else return "none";
  }() );
  if( WarpX::do_dive_cleaning )
    meshes.setAttribute("chargeCorrectionParameters", "period=1");

  // Loop through the different components, i.e. different fields stored in mf
  for (int icomp=0; icomp<ncomp; icomp++){

    // Check if this field is a vector or a scalar, and extract the field name
    std::string const & varname = varnames[icomp];
    std::string field_name = varname;
    std::string comp_name = openPMD::MeshRecordComponent::SCALAR;
    for( char const* vector_field: {"E", "B", "j"} ) {
        for( char const* comp: {"x", "y", "z"} ) {
            if( varname[0] == *vector_field && varname[1] == *comp ) {
                field_name = varname[0] + varname.substr(2); // Strip component
                comp_name = varname[1];
            }
        }
    }

    // Setup the mesh record accordingly
    auto mesh = meshes[field_name];
    mesh.setDataOrder( openPMD::Mesh::DataOrder::F ); // MultiFab: Fortran order of indices and axes
    mesh.setAxisLabels( axis_labels );
    mesh.setGridSpacing( grid_spacing );
    mesh.setGridGlobalOffset( global_offset );
    mesh.setAttribute( "fieldSmoothing", "none" );
    setOpenPMDUnit( mesh, field_name );

    // Create a new mesh record component, and store the associated metadata
    auto mesh_comp = mesh[comp_name];
    mesh_comp.resetDataset( dataset );
    // Cell-centered data: position is at 0.5 of a cell size.
    mesh_comp.setPosition( std::vector<double>{AMREX_D_DECL(0.5, 0.5, 0.5)} );

    // Loop through the multifab, and store each box as a chunk,
    // in the openPMD file.
    for( amrex::MFIter mfi(mf); mfi.isValid(); ++mfi ) {

      amrex::FArrayBox const& fab = mf[mfi];
      amrex::Box const& local_box = fab.box();

      // Determine the offset and size of this chunk
      amrex::IntVect const box_offset = local_box.smallEnd() - global_box.smallEnd();
      auto const chunk_offset = getReversedVec( box_offset );
      auto const chunk_size = getReversedVec( local_box.size() );

      // Write local data
      amrex::Real const * local_data = fab.dataPtr( icomp );
      mesh_comp.storeChunk( openPMD::shareRaw(local_data),
                            chunk_offset, chunk_size );
    }
  }
  // Flush data to disk after looping over all components
  m_Series->flush();
}
#endif // WARPX_USE_OPENPMD



//
//
//
WarpXParticleCounter::WarpXParticleCounter(const std::unique_ptr<WarpXParticleContainer>& pc)
{
  m_MPISize = amrex::ParallelDescriptor::NProcs();
  m_MPIRank = amrex::ParallelDescriptor::MyProc();

  m_ParticleCounterByLevel.resize(pc->finestLevel()+1);
  m_ParticleOffsetAtRank.resize(pc->finestLevel()+1);
  m_ParticleSizeAtRank.resize(pc->finestLevel()+1);

  for (auto currentLevel = 0; currentLevel <= pc->finestLevel(); currentLevel++)
    {
      long numParticles = 0; // numParticles in this processor

      for (WarpXParIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
    auto numParticleOnTile = pti.numParticles();
    numParticles += numParticleOnTile;
      }

      unsigned long long offset=0; // offset of this level
      unsigned long long sum=0; // numParticles in this level (sum from all processors)

      GetParticleOffsetOfProcessor(numParticles, offset,  sum);

      m_ParticleCounterByLevel[currentLevel] = sum;
      m_ParticleOffsetAtRank[currentLevel] = offset;
      m_ParticleSizeAtRank[currentLevel] = numParticles;

      // adjust offset, it should be numbered after particles from previous levels
      for (auto lv=0; lv<currentLevel; lv++)
    m_ParticleOffsetAtRank[currentLevel] += m_ParticleCounterByLevel[lv];

      m_Total += sum;
    }
}


// get the offset in the overall particle id collection
//
// note: this is a MPI-collective operation
//
// input: num of particles  of from each   processor
//
// output:
//     offset within <all> the particles in the comm
//     sum of all particles in the comm
//
void
WarpXParticleCounter::GetParticleOffsetOfProcessor(const long& numParticles,
                           unsigned long long& offset,
                           unsigned long long& sum)  const


{
      std::vector<long> result(m_MPISize,  0);
      amrex::ParallelGather::Gather (numParticles, result.data(), -1, amrex::ParallelDescriptor::Communicator());

      sum = 0;
      offset = 0;
      for (int i=0;  i<result.size();  i++) {
    sum +=  result[i];
    if (i<m_MPIRank)
      offset +=  result[i];
      }
}
