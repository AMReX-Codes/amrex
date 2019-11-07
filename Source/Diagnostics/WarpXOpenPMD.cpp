#include "WarpXOpenPMD.H"
#include "FieldIO.H"  // for getReversedVec

WarpXOpenPMDParticle::WarpXOpenPMDParticle(bool oneFilePerTS,
                       const std::string& openPMDFileType)
  :m_Series(nullptr),
   m_OneFilePerTS(oneFilePerTS),
   m_OpenPMDFileType(openPMDFileType)
{}

/*
WarpXOpenPMDParticle::WarpXOpenPMDParticle(const std::string& filename,
                       const std::string& openPMDFileType)
  :m_Series(nullptr),
   m_OneFilePerTS(true),
   m_OpenPMDFileType(openPMDFileType)
{
  Init(filename, openPMD::AccessType::CREATE);
}
*/

WarpXOpenPMDParticle::~WarpXOpenPMDParticle()
{
  if (nullptr != m_Series) {
    m_Series->flush();
    delete m_Series;
  }
}


//
//
//
void WarpXOpenPMDParticle::GetFileName(std::string& filename, int ts)
{
  std::string dir = "diags/";
  std::string prefix="/plot";

  // Write openPMD format: only for level 0
  filename = std::string(dir);
  filename.append(m_OpenPMDFileType);

  if (m_OneFilePerTS)
    filename += amrex::Concatenate(prefix, ts);
  else
    filename.append(prefix);

  filename += std::string(".") + m_OpenPMDFileType;
}


void WarpXOpenPMDParticle::SetOpenPMDType(const std::string& type)
{
  m_OpenPMDFileType = type;
}

void WarpXOpenPMDParticle::SetStep(int ts)
{
  if (ts < 0)
    return;

  m_CurrentStep =  ts;

  std::string filename;
  if (m_OneFilePerTS) {
    GetFileName(filename, ts);
    Init(filename, openPMD::AccessType::CREATE);
  } else { // one file
    if (nullptr == m_Series) {
      GetFileName(filename, -1);
      Init(filename, openPMD::AccessType::CREATE);
    }
  }
}

void
WarpXOpenPMDParticle::Init(const std::string& filename,
               openPMD::AccessType accessType)
{
  if (m_Series != nullptr) {
    m_Series->flush();
    delete m_Series;
    m_Series = nullptr;
  }

  if (amrex::ParallelDescriptor::NProcs() > 1) {
    m_Series = new openPMD::Series(filename,
                   accessType,
                   amrex::ParallelDescriptor::Communicator());
    m_MPISize = amrex::ParallelDescriptor::NProcs();
    m_MPIRank = amrex::ParallelDescriptor::MyProc();
  }
  else
    m_Series = new openPMD::Series(filename, accessType);
}


void
WarpXOpenPMDParticle::SaveContainerPlots(const std::unique_ptr<MultiParticleContainer>& mpc)
{
  std::vector<std::string> species_names = mpc->GetSpeciesNames();

  for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
    //auto& pc = allcontainers[i];
    //auto& pc =  mpc->GetParticleContainer(i);
    auto& pc  = mpc->GetUniqueContainer(i);
    if (pc->plot_species) {

      amrex::Vector<std::string> real_names;
      amrex::Vector<std::string> int_names;
      amrex::Vector<int> int_flags;

      real_names.push_back("weight");

      real_names.push_back("momentum_x");
      real_names.push_back("momentum_y");
      real_names.push_back("momentum_z");

      real_names.push_back("Ex");
      real_names.push_back("Ey");
      real_names.push_back("Ez");

      real_names.push_back("Bx");
      real_names.push_back("By");
      real_names.push_back("Bz");

#ifdef WARPX_DIM_RZ
      real_names.push_back("theta");
#endif
      if(pc->do_field_ionization){
    int_names.push_back("ionization_level");
    // int_flags specifies, for each integer attribs, whether it is
    // dumped to plotfiles. So far, ionization_level is the only
    // integer attribs, and it is automatically dumped to plotfiles
    // when ionization is on.
    int_flags.resize(1, 1);
      }

      // Convert momentum to SI
      pc->ConvertUnits(ConvertDirection::WarpX_to_SI);
      // real_names contains a list of all particle attributes.
      // pc->plot_flags is 1 or 0, whether quantity is dumped or not.

      //
      // NOTE: using 0 as iteration number b/c right now
      // because warpx writes one file per iteration in this function
      // so the iteration id in openPMD-api does not matter
      // Needs to add MultipleParticleContainter::WritePlotFile(string, iteration)  and define
      // "openPMDWriter" in the WarpXIO.cpp before calling WritePlotFile(string, iter)
      //
      //std::unique_ptr<WarpXParticleContainer> upc (&pc);
      SavePlotFile(pc,
           species_names[i],
           m_CurrentStep, //-1, // use -1 so levels will be iterations
           pc->plot_flags, // this is protected and accessible by MultiParticleContainer.
           // so kept as is
           int_flags,
           real_names, int_names);

      // Convert momentum back to WarpX units
      pc->ConvertUnits(ConvertDirection::SI_to_WarpX);
    }
  }
}


void
WarpXOpenPMDParticle::SavePlotFile (const std::unique_ptr<WarpXParticleContainer>& pc,
                    const std::string& name,
                    int iteration,
                    const amrex::Vector<int>& write_real_comp,
                    const amrex::Vector<int>& write_int_comp,
                    const amrex::Vector<std::string>& real_comp_names,
                    const amrex::Vector<std::string>&  int_comp_names) const
{
  if ( nullptr == m_Series)
    return;

  if (0 == m_CurrentStep ) // can only define once
    m_Series->setParticlesPath("TestingWarpXParticles_one_level");

  int numLevels = pc->finestLevel();

  for (auto currentLevel = 0; currentLevel <= pc->finestLevel(); currentLevel++)
    {
      //openPMD::Iteration& currIteration = m_Series->iterations[currentLevel];
      std::ostringstream s;       s << name;
      //int iter = iteration;
      s <<":amrLevel="<<currentLevel;


      //openPMD::Iteration& currIteration = m_Series->iterations[iter];
      openPMD::Iteration& currIteration = m_Series->iterations[iteration];
      //openPMD::ParticleSpecies& currSpecies = currIteration.particles[name];

      std::string leveledName(s.str());
      openPMD::ParticleSpecies& currSpecies = currIteration.particles[leveledName];

      long numParticles = 0;

      for (WarpXParIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
    auto numParticleOnTile = pti.numParticles();
    numParticles += numParticleOnTile;
      }

      unsigned long offset=0;
      unsigned long long sum=0;

      GetParticleOffsetOfProcessor(numParticles, offset,  sum);

      //
      // define positions & offsets
      //
      SetupPos(currSpecies, sum);
      SetupRealProperties(currSpecies, write_real_comp, real_comp_names, sum);

      //if return after this, all is fine (although nothing useful is written)

      if (0 == numParticles)
          return;

      // pc->NumIntComp() & NumRealComp() are protected,
      // from WarpXParIter template class definition, we know that
      // soa num real attributes = PIdx::nattribs, and num int in soa is 0

      for (WarpXParIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
    auto numParticleOnTile = pti.numParticles();

    // get position from aos
    const auto& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile

#if (AMREX_SPACEDIM == 3)
    std::vector<amrex::ParticleReal> currX(numParticleOnTile, 0);
    std::vector<amrex::ParticleReal> currY(numParticleOnTile, 0);
    std::vector<amrex::ParticleReal> currZ(numParticleOnTile, 0);

    for (auto i=0; i<numParticleOnTile; i++) {
      currX[i] = aos[i].m_rdata.pos[0];
      currY[i] = aos[i].m_rdata.pos[1];
      currZ[i] = aos[i].m_rdata.pos[2];
    }

    currSpecies["position"]["x"].storeChunk(currX, {offset}, {static_cast<unsigned long long>(numParticleOnTile)});
    currSpecies["position"]["y"].storeChunk(currY, {offset}, {static_cast<unsigned long long>(numParticleOnTile)});
    currSpecies["position"]["z"].storeChunk(currZ, {offset}, {static_cast<unsigned long long>(numParticleOnTile)});

    SaveRealProperty(pti,
             currSpecies,
             offset, numParticles,
             write_real_comp, real_comp_names);
#endif

    offset += numParticleOnTile;
      }
    }
}
void
WarpXOpenPMDParticle::SetupRealProperties(openPMD::ParticleSpecies& currSpecies,
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
      auto& particleVar = currSpecies[real_comp_names[i]];
      auto& particleVarComp = particleVar[openPMD::RecordComponent::SCALAR];
      particleVarComp.resetDataset(particlesLineup);
    }
}

void
WarpXOpenPMDParticle::SaveRealProperty(WarpXParIter& pti,
                       openPMD::ParticleSpecies& currSpecies,
                       unsigned long offset,
                       unsigned long long numParticles,
                       const amrex::Vector<int>& write_real_comp,
                       const amrex::Vector<std::string>& real_comp_names) const

{
  int numOutputReal = 0;
  int totalRealAttrs = m_NumAoSRealAttributes + m_NumSoARealAttributes;

  for (int i = 0; i < totalRealAttrs; ++i)
    if (write_real_comp[i])
      ++numOutputReal;

  auto numParticleOnTile = pti.numParticles();
  const auto& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
  const auto& soa = pti.GetStructOfArrays();

  // properties are saved seperately
  {
    for (auto idx=0; idx<m_NumAoSRealAttributes; idx++) {
      if (write_real_comp[idx]) {
    auto& currVar = currSpecies[real_comp_names[idx]][openPMD::RecordComponent::SCALAR];
    typename amrex::ParticleReal *d =
      static_cast<typename amrex::ParticleReal*> (malloc(sizeof(typename amrex::ParticleReal) *  numParticleOnTile));

    for (auto kk=0; kk<numParticleOnTile; kk++)
      d[kk] = aos[kk].m_rdata.arr[AMREX_SPACEDIM+idx];

    std::shared_ptr <typename amrex::ParticleReal> data(d, free);
    currVar.storeChunk(data,
               {offset}, {static_cast<unsigned long long>(numParticleOnTile)});

      }
    }
  }

  {
    for (auto idx=0; idx<m_NumSoARealAttributes; idx++) {
      auto ii = m_NumAoSRealAttributes + idx;
      if (write_real_comp[ii]) {
    auto& currVar = currSpecies[real_comp_names[ii]][openPMD::RecordComponent::SCALAR];
    currVar.storeChunk(openPMD::shareRaw(soa.GetRealData(idx)),
               {offset}, {static_cast<unsigned long long>(numParticleOnTile)});
      }
    }
  }
}



void
WarpXOpenPMDParticle::SetupPos(openPMD::ParticleSpecies& currSpecies,
                   const unsigned long long& np) const
{
  auto particleLineup = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});

  currSpecies["positionOffset"]["x"].resetDataset(particleLineup);
  currSpecies["positionOffset"]["x"].makeConstant(0);
  currSpecies["positionOffset"]["y"].resetDataset(particleLineup);
  currSpecies["positionOffset"]["y"].makeConstant(0);
  currSpecies["positionOffset"]["z"].resetDataset(particleLineup);
  currSpecies["positionOffset"]["z"].makeConstant(0);

  //auto positions = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});
  currSpecies["position"]["x"].resetDataset(particleLineup);
  currSpecies["position"]["y"].resetDataset(particleLineup);
  currSpecies["position"]["z"].resetDataset(particleLineup);
}


//
// input: num of particles  of from each   processor
//
// output:
//     offset within <all> the particles in the comm
//     sum of all particles in the comm
//
void
WarpXOpenPMDParticle::GetParticleOffsetOfProcessor(const long& numParticles,
                           unsigned long& offset,
                           unsigned long long& sum) const
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
  //std::cout <<"    rank: "<<m_MPIRank<<"  offset:  "<<offset<<", sum="<<sum<<std::endl;
}




//
// this is originally copied from FieldIO.cpp
//
void
WarpXOpenPMDParticle::WriteOpenPMDFields( //const std::string& filename,
                      const std::vector<std::string>& varnames,
                      const amrex::MultiFab& mf,
                      const amrex::Geometry& geom,
                      const int iteration,
                      const double time ) const
{
  BL_PROFILE("WriteOpenPMDFields()");

  if ( nullptr == m_Series)
    return;

  const int ncomp = mf.nComp();

  // Create a few vectors that store info on the global domain
  // Swap the indices for each of them, since AMReX data is Fortran order
  // and since the openPMD API assumes contiguous C order
  // - Size of the box, in integer number of cells
  const amrex::Box& global_box = geom.Domain();
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
  openPMD::Datatype datatype = openPMD::determineDatatype<amrex::Real>();
  auto dataset = openPMD::Dataset(datatype, global_size);

  // Create new file and store the time/iteration info
  /*
  auto series = [filename](){
      if( ParallelDescriptor::NProcs() > 1 )
          return openPMD::Series( filename,
                                  openPMD::AccessType::CREATE,
                                  ParallelDescriptor::Communicator() );
      else
          return openPMD::Series( filename,
                                  openPMD::AccessType::CREATE );
  }();
  */
  auto series_iteration = m_Series->iterations[iteration];
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
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {

      const amrex::FArrayBox& fab = mf[mfi];
      const amrex::Box& local_box = fab.box();

      // Determine the offset and size of this chunk
      amrex::IntVect box_offset = local_box.smallEnd() - global_box.smallEnd();
      auto chunk_offset = getReversedVec(box_offset);
      auto chunk_size = getReversedVec(local_box.size());

      // Write local data
      const double* local_data = fab.dataPtr(icomp);
      mesh_comp.storeChunk(openPMD::shareRaw(local_data),
                           chunk_offset, chunk_size);
    }
  }
  // Flush data to disk after looping over all components
  //std::cout<<" this is optional "<<std::endl;
  m_Series->flush();
}
