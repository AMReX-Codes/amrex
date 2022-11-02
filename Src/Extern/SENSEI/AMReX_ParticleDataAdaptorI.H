#include "Profiler.h"
#include "Error.h"
#include "SVTKUtils.h"
#include "MeshMetadata.h"
// svtk includes
#include <svtkPolyData.h>
#include <svtkDataSetAttributes.h>
#include <svtkCellData.h>
#include <svtkPointData.h>
#include <svtkMultiBlockDataSet.h>



// local includes
#include "AMReX_ParticleDataAdaptor.H"
#include <AMReX_InSituUtils.H>

namespace amrex
{
//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>*
ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::New()
{
  auto result = new ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  result->InitializeObjectBase();
  return result;
}
// senseiNewMacro(ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>);

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetDataSource(
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
    const std::map<std::string, std::vector<int>> & rStructs,
    const std::map<std::string, int> & iStructs,
    const std::map<std::string, std::vector<int>> & rArrays,
    const std::map<std::string, int> & iArrays)
{
  // set the data source on the particles
  this->ReleaseData();
  this->m_particles = particles;

  // set the array names
  int ret = this->SetArrayNames(rStructs, iStructs, rArrays, iArrays);
  if(ret)
  {
    SENSEI_ERROR("problem with array names in SetDataSource");
    return ret;
  }
  return ret;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetArrayNames(
  const std::map<std::string, std::vector<int>> & rStructs,
  const std::map<std::string, int> & iStructs,
  const std::map<std::string, std::vector<int>> & rArrays,
  const std::map<std::string, int> & iArrays)
{
    if(rStructs.size() <= NStructReal)
  {
    // check that no indices from the rStructs map exceed allowable value
    for(auto s : rStructs)
    {
      for(auto i : s.second)
      {
        if(i >= NStructReal)
        {
          SENSEI_ERROR("rStruct index exceeds internal storage size");
          return -1;
        }
      }
    }
    m_realStructs = rStructs;
  }
  else
  {
    SENSEI_ERROR("rStructs array size exceeds internal storage size");
    return -1;
  }

  if(iStructs.size() <= NStructInt)
  {
    // check that no indices from the iStructs map exceed allowable value
    for(auto s : iStructs)
    {
      if(s.second >= NStructInt)
        {
          SENSEI_ERROR("iStructs index exceeds internal storage size");
          return -1;
        }
    }
    m_intStructs = iStructs;
  }
  else
  {
    SENSEI_ERROR("iStructs array size exceeds internal storage size");
    return -1;
  }


  if(rArrays.size() <= NArrayReal)
  {
    // check that no indices from the rArrays map exceed allowable value
    for(auto s : rArrays)
    {
      for(auto i : s.second)
      {
        if(i >= NArrayReal)
        {
          SENSEI_ERROR("rArrays index exceeds internal storage size");
          return -1;
        }
      }
    }
    m_realArrays = rArrays;
  }
  else
  {
    SENSEI_ERROR("rArrays array size exceeds internal storage size");
    return -1;
  }
  if(iArrays.size() <= NArrayInt)
  {
    // check that no indices from the iArrays map exceed allowable value
    for(auto s : iArrays)
    {
      if(s.second >= NArrayInt)
        {
          SENSEI_ERROR("iArray index exceeds internal storage size");
          return -1;
        }
    }
    m_intArrays = iArrays;
  }
  else
  {
    SENSEI_ERROR("iArrays array size exceeds internal storage size");
    return -1;
  }
  return 0;
}
//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetNumberOfMeshes(unsigned int &numMeshes)
{
  numMeshes = 1;
  return 0;
}

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshName(
  unsigned int id, std::string &meshName)
{
  meshName = m_particlesName;
  return 0;
}
#endif

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshHasGhostNodes(
  const std::string &meshName,
  int &nLayers)
{
  nLayers = 0;
  return 0;
}
#endif

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshHasGhostCells(
  const std::string &meshName,
  int &nLayers)
{
  nLayers = 0;
  return 0;
}
#endif

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetNumberOfArrays(
  const std::string &meshName,
  int association,
  unsigned int &numberOfArrays)
{
  numberOfArrays = 0;
  if(association == svtkDataObject::POINT)
  {
    numberOfArrays = m_realStructs.size()
                   + m_intStructs.size()
                   + m_realArrays.size()
                   + m_intArrays.size();
  }
  return 0;
}
#endif
//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetArrayName(
  const std::string &meshName,
  int association,
  unsigned int index,
  std::string &arrayName)
{
  if(association == svtkDataObject::POINT)
  {
    if(index < m_realStructs.size())
    {
      auto a = m_realStructs.begin() + index;
      arrayName = *(a).first;
      return 0;
    }
    if(index < m_intStructs.size())
    {
      int ind = index - m_realStructs.size();
      auto a = m_intStructs.begin() + ind;
      arrayName = *(a).first;
      return 0;
    }
    if(index < m_realArrays.size())
    {
      int ind = index - m_realStructs.size() - m_intStructs.size();
      auto a = m_realArrays.begin() + ind;
      arrayName = *(a).first;
      return 0;
    }
    if(index < m_intArrays.size())
    {
      int ind = index - m_realStructs.size() - m_intStructs.size() - m_realArrays.size();
      auto a = m_intArrays.begin() + ind;
      arrayName = *(a).first;
      return 0;
    }
  }

  return -1;
}
#endif

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMesh(
  const std::string &meshName,
  bool structureOnly,
  svtkDataObject *&mesh)
{
  mesh = nullptr;
  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(this->GetCommunicator(), &nprocs);
  MPI_Comm_rank(this->GetCommunicator(), &rank);

  if (meshName != m_particlesName)
  {
    SENSEI_ERROR("No mesh named \"" << meshName << "\"")
    return -1;
  }
  svtkMultiBlockDataSet* mb = svtkMultiBlockDataSet::New();

  if (structureOnly)
  {
    mesh = mb;
    return 0;
  }

  mb->SetNumberOfBlocks(nprocs);
  svtkPolyData *pd = BuildParticles();
  mb->SetBlock(rank, pd);
  pd->Delete();
  mesh = mb;

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddGhostNodesArray(
  svtkDataObject*,
  const std::string &meshName)
{
  if (meshName != m_particlesName)
  {
    SENSEI_ERROR("no mesh named \"" << meshName << "\"")
    return -1;
  }
  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddGhostCellsArray(
  svtkDataObject*,
  const std::string &meshName)
{
  if (meshName != m_particlesName)
  {
    SENSEI_ERROR("no mesh named \"" << meshName << "\"")
    return -1;
  }
  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddArray(
  svtkDataObject* mesh,
  const std::string &meshName,
  int association,
  const std::string &arrayName)
{
  if (meshName != m_particlesName)
  {
    SENSEI_ERROR("no mesh named \"" << meshName << "\"");
    return -1;
  }

  if (association != svtkDataObject::POINT)
  {
    SENSEI_ERROR("Invalid association " << association);
    return -1;
  }

  if(m_realStructs.find(arrayName) != m_realStructs.end())
  {
    return this->AddParticlesAOSRealArray(arrayName, mesh);
  }

  if(m_intStructs.find(arrayName) != m_intStructs.end())
  {
    return this->AddParticlesAOSIntArray(arrayName, mesh);
  }

  if(m_realArrays.find(arrayName) != m_realArrays.end())
  {
    return this->AddParticlesSOARealArray(arrayName, mesh);
  }

  if(m_intArrays.find(arrayName) != m_intArrays.end())
  {
    return this->AddParticlesSOAIntArray(arrayName, mesh);
  }

  SENSEI_ERROR("Invalid array name " << arrayName);
  return -1;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::ReleaseData()
{
  this->m_particles = nullptr;
  return 0;
}

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR >= 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshMetadata(
  unsigned int id,
  sensei::MeshMetadataPtr &metadata)
{
  sensei::TimeEvent<64> event("AmrMeshDataAdaptor::GetMeshMetadata");

  int nprocs = 1;
  int rank = 0;
  MPI_Comm_size(this->GetCommunicator(), &nprocs);
  MPI_Comm_rank(this->GetCommunicator(), &rank);
  if (id != 0 && id != 1)
  {
    SENSEI_ERROR("invalid mesh id " << id)
    return -1;
  }


  // AMR data is always expected to be a global view
  metadata->GlobalView = false;      // tells if the information describes data
                                     // on this rank or all ranks. Passed into
                                     // Set methods, Get methods generate the
                                     // desired view.

  // name of mesh (all)
  metadata->MeshName = m_particlesName;

  // container mesh type (all)
  metadata->MeshType = SVTK_MULTIBLOCK_DATA_SET;

  // block mesh type (all)
  metadata->BlockType = SVTK_POLY_DATA;

  // global number of blocks (all)
  metadata->NumBlocks = nprocs;

  // number of blocks on each rank (all)
  metadata->NumBlocksLocal = {1};

  // global cell index space extent (Cartesian, AMR, optional)
  // std::array<int,6> Extent;

  // global bounding box (all, optional)
  // std::array<double,6> Bounds;

  // type enum of point data (unstructured, optional)
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  metadata->CoordinateType = SVTK_FLOAT;
#else
  metadata->CoordinateType = SVTK_DOUBLE;
#endif

  // total number of points in all blocks (all, optional)
  // long NumPoints;

  // total number of cells in all blocks (all, optional)
  // long NumCells;

  // total cell array size in all blocks (all, optional)
  // long CellArraySize;

  // number of arrays (all)
  metadata->NumArrays = m_realStructs.size()
                      + m_intStructs.size()
                      + m_realArrays.size()
                      + m_intArrays.size();

  // number of ghost cell layers (all)
  metadata->NumGhostCells = 0;

  // number of ghost node layers (all)
  metadata->NumGhostNodes = 0;

  // number of AMR levels (AMR)
  // metadata->NumLevels;

  // non zero if the mesh does not change in time (all)
  metadata->StaticMesh = 0;

  // name of each data array (all)
  metadata->ArrayName = {};
  for(auto s : m_realStructs)
  {
    metadata->ArrayName.push_back(s.first);
  }
  for(auto s : m_intStructs)
  {
    metadata->ArrayName.push_back(s.first);
  }
  for(auto s : m_realArrays)
  {
    metadata->ArrayName.push_back(s.first);
  }
  for(auto s : m_intArrays)
  {
    metadata->ArrayName.push_back(s.first);
  }

  // centering of each data array (all)
  metadata->ArrayCentering = {};
  for(auto s : m_realStructs)
  {
    metadata->ArrayCentering.push_back(svtkDataObject::POINT);
  }
  for(auto s : m_intStructs)
  {
    metadata->ArrayCentering.push_back(svtkDataObject::POINT);
  }
  for(auto s : m_realArrays)
  {
    metadata->ArrayCentering.push_back(svtkDataObject::POINT);
  }
  for(auto s : m_intArrays)
  {
    metadata->ArrayCentering.push_back(svtkDataObject::POINT);
  }

  // number of components of each array (all)
  metadata->ArrayComponents = {};
  for(auto s : m_realStructs)
  {
    metadata->ArrayComponents.push_back(s.second.size());
  }
  for(auto s : m_intStructs)
  {
    metadata->ArrayComponents.push_back(1);
  }
  for(auto s : m_realArrays)
  {
    metadata->ArrayComponents.push_back(s.second.size());
  }
  for(auto s : m_intArrays)
  {
    metadata->ArrayComponents.push_back(1);
  }

  // type enum of each data array (all)
  metadata->ArrayType = {};
  for(auto s : m_realStructs)
  {
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
    metadata->ArrayType.push_back(SVTK_FLOAT);
#else
    metadata->ArrayType.push_back(SVTK_DOUBLE);
#endif
  }
  for(auto s : m_intStructs)
  {
    metadata->ArrayType.push_back(SVTK_INT);
  }
  for(auto s : m_realArrays)
  {
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
    metadata->ArrayType.push_back(SVTK_FLOAT);
#else
    metadata->ArrayType.push_back(SVTK_DOUBLE);
#endif
  }
  for(auto s : m_intArrays)
  {
    metadata->ArrayType.push_back(SVTK_INT);
  }

  // global min,max of each array (all, optional)
  // metadata->ArrayRange = {};

  // rank where each block resides (all, optional)
  metadata->BlockOwner = {rank};

  // global id of each block (all, optional)
  metadata->BlockIds = metadata->BlockOwner;

  // number of points for each block (all, optional)
  auto nptsOnProc = this->m_particles->TotalNumberOfParticles(true, true);
  metadata->BlockNumPoints = {nptsOnProc};

  // number of cells for each block (all, optional)
  metadata->BlockNumCells = {nptsOnProc};

  // cell array size for each block (unstructured, optional)
  metadata->BlockCellArraySize = {nptsOnProc};

  // note: for AMR BlockExtents and BlockBounds are always global
  // index space extent of each block [i0,i1, j0,j1, k0,k1] (Cartesian, AMR, optional)
  // std::vector<std::array<int,6>> BlockExtents;

  // bounds of each block [x0,x1, y0,y1, z0,z1] (all, optional)
  // std::vector<std::array<double,6>> BlockBounds;
  if (metadata->Flags.BlockBoundsSet())
  {
    std::array<double,6> bounds({
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::lowest(),
#if (AMREX_SPACEDIM == 1)
      0.0, 0.0, 0.0, 0.0
#elif (AMREX_SPACEDIM == 2)
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::lowest(),
      0.0, 0.0
#elif (AMREX_SPACEDIM == 3)
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::lowest()
#endif
      });

    // loop over levels in the particle container
    const auto& particles = this->m_particles->GetParticles();
    for (int lev = 0; lev < particles.size();  lev++)
    {
      // loop over ParticleTiles on level
      auto& pmap = particles[lev];
      for (const auto& kv : pmap)
      {
        // loop over particles in ParticleTile
        const auto& aos = kv.second.GetArrayOfStructs();
        const auto &parts = aos();

        // visit only the "real" particles, skip the "neighbor" particles.
        long long numReal = aos.numRealParticles();
        for (long long i = 0; i < numReal; ++i)
        {
          const auto &part = parts[i];
          if (part.id() > 0)
          {
#if (AMREX_SPACEDIM == 1)
            bounds[0] = bounds[0] > part.pos(0) ? part.pos(0) : bounds[0];
            bounds[1] = bounds[1] < part.pos(0) ? part.pos(0) : bounds[1];
#elif (AMREX_SPACEDIM == 2)
            bounds[0] = bounds[0] > part.pos(0) ? part.pos(0) : bounds[0];
            bounds[1] = bounds[1] < part.pos(0) ? part.pos(0) : bounds[1];
            bounds[2] = bounds[2] > part.pos(1) ? part.pos(1) : bounds[2];
            bounds[3] = bounds[3] < part.pos(1) ? part.pos(1) : bounds[3];
#elif (AMREX_SPACEDIM == 3)
            bounds[0] = bounds[0] > part.pos(0) ? part.pos(0) : bounds[0];
            bounds[1] = bounds[1] < part.pos(0) ? part.pos(0) : bounds[1];
            bounds[2] = bounds[2] > part.pos(1) ? part.pos(1) : bounds[2];
            bounds[3] = bounds[3] < part.pos(1) ? part.pos(1) : bounds[3];
            bounds[4] = bounds[4] > part.pos(2) ? part.pos(2) : bounds[4];
            bounds[5] = bounds[5] < part.pos(2) ? part.pos(2) : bounds[5];
#endif
          }
        }
      }
    }

    metadata->BlockBounds = {bounds};
  }

  // min max of each array on each block.
  // indexed by block then array. (all, optional)
  // std::vector<std::vector<std::array<double,2>>> BlockArrayRange;
  if (metadata->Flags.BlockArrayRangeSet())
  {
    SENSEI_ERROR("BlockArrayRange requested but not yet implemented.")
  }

  // refinement ratio in i,j, and k directions for each level (AMR)
  // std::vector<std::array<int,3>> RefRatio;

  // number of blocks in each level (AMR)
  // std::vector<int> BlocksPerLevel;

  // AMR level of each block (AMR)
  // std::vector<int> BlockLevel;

  // flag indicating presence of a periodic boundary in the i,j,k direction (all)
  // std::array<int,3> PeriodicBoundary;

#if defined(SENSEI_DEBUG)
  metadata->Validate(this->GetCommunicator(), sensei::MeshMetadataFlags());
  metadata->ToStream(std::cerr);
#endif

  return 0;
}
#endif

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
svtkPolyData* ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::BuildParticles()
{
  // return particle data pd
  svtkPolyData* pd  = svtkPolyData::New();

  const auto& particles = this->m_particles->GetParticles();
  long long numParticles = this->m_particles->TotalNumberOfParticles(true, true);

  // allocate vertex storage for particles
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  svtkNew<svtkFloatArray> coords;
#else
  svtkNew<svtkDoubleArray> coords;
#endif
  coords->SetName("coords");
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(numParticles);
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  float *pCoords = coords->GetPointer(0);
#else
  double *pCoords = coords->GetPointer(0);
#endif

  // use this to index into the SVTK array as we copy level by level and tile by
  // tile
  long long ptId = 0;

  // allocate connectivity array for particles
  svtkNew<svtkCellArray> vertex;
  vertex->AllocateExact(numParticles, 1);

  // points->SetNumberOfPoints(numParticles);

  // loop over levels in the particle container
  for (int lev = 0; lev < particles.size();  lev++)
  {
    using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    for (MyParIter pti(*this->m_particles, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto numReal = pti.numParticles();

      for (long long i = 0; i < numReal; ++i)
      {
        const auto &part = aos[i];
        if (part.id() > 0)
        {
          // add a vertex type cell
          vertex->InsertNextCell(1);
          vertex->InsertCellPoint(ptId);
          // copy the partilce coordinates
#if (AMREX_SPACEDIM == 1)
          pCoords[0] = part.pos(0);
          pCoords[1] = amrex_particle_real(0);
          pCoords[2] = amrex_particle_real(0);
#elif (AMREX_SPACEDIM == 2)
          pCoords[0] = part.pos(0);
          pCoords[1] = part.pos(1);
          pCoords[2] = amrex_particle_real(0);
#elif (AMREX_SPACEDIM == 3)
          pCoords[0] = part.pos(0);
          pCoords[1] = part.pos(1);
          pCoords[2] = part.pos(2);
#endif
          pCoords += 3;
          ++ptId;
        }
      }
    }
  }

  // pass the particle coordinates into SVTK's point data structure.
  svtkNew<svtkPoints> points;
  points->SetData(coords);

  // add point and vertex data to output mesh
  pd->SetPoints(points);
  pd->SetVerts(vertex);

  return pd;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesIDArray(
  svtkDataObject* mesh)
{
  auto svtk_particles = dynamic_cast<svtkPolyData*>(mesh);
  const auto& particles = this->m_particles->GetParticles();
  auto nptsOnProc = this->m_particles->TotalNumberOfParticles(true, true);

 // allocate a SVTK array for the data
  svtkNew<svtkIntArray> idArray;
  idArray->SetName("id");
  idArray->SetNumberOfComponents(1);
  idArray->SetNumberOfValues(nptsOnProc);

  // get pointer underlying to data
  int *partIds = idArray->GetPointer(0);

  // loop over particles and append their cpu value to the list
  using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  long ptId = 0;
  for (int level = 0; level < particles.size(); ++level)
  {
    for (MyParIter pti(*this->m_particles, level); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto numReal = pti.numParticles();
      for (long long i = 0; i < numReal; ++i)
      {
        const auto &part = aos[i];
        if (part.id() > 0)
        {
          partIds[i] = part.id();
        }
      }
      partIds += numReal;
    }
  }

  // the association for this array is svtkDataObject::POINT
  svtk_particles->GetPointData()->AddArray(idArray);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesCPUArray(
  svtkDataObject* mesh)
{
  auto svtk_particles = dynamic_cast<svtkPolyData*>(mesh);
  const auto& particles = this->m_particles->GetParticles();
  auto nptsOnProc = this->m_particles->TotalNumberOfParticles(true, true);

  // allocate a SVTK array for the data
  svtkNew<svtkIntArray> cpuArray;
  cpuArray->SetName("cpu");
  cpuArray->SetNumberOfComponents(1);
  cpuArray->SetNumberOfValues(nptsOnProc);

  // get pointer to underlying data
  int* partCpu = cpuArray->GetPointer(0);

  // loop over particles and append their cpu value to the list
  using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  for (int level = 0; level < particles.size(); ++level)
  {
    for (MyParIter pti(*this->m_particles, level); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto numReal = pti.numParticles();
      for (long long i = 0; i < numReal; ++i)
      {
        const auto &part = aos[i];
        if (part.id() > 0)
        {
          partCpu[i] = part.cpu();
        }
      }
      partCpu += numReal;
    }
  }

  // the association for this array is svtkDataObject::POINT
  svtk_particles->GetPointData()->AddArray(cpuArray);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesSOARealArray(
  const std::string &arrayName,
  svtkDataObject* mesh)
{
  const long nParticles = this->m_particles->TotalNumberOfParticles(true, true);

  // check that the name of the arrays is listed in the m_realArrays
  RealDataMapType::iterator ait = m_realArrays.find(arrayName);
  if (ait == m_realArrays.end())
  {
    SENSEI_ERROR("No real SOA named \"" << arrayName << "\"");
    return -1;
  }

  // get the indices of the structs on the particles
  const std::vector<int> &indices = ait->second;
  int nComps = indices.size();

  // check that the indices are within the templated storage spec
  for(auto i : indices)
  {
    if(i >= NArrayReal)
    {
      SENSEI_ERROR("Index out of bounds for real SOA named \"" << arrayName << "\"");
      return -1;
    }
  }

  // allocate the svtkArray
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  svtkNew<svtkFloatArray> data;
#else
  svtkNew<svtkDoubleArray> data;
#endif
  data->SetName(arrayName.c_str());
  data->SetNumberOfComponents(nComps);
  data->SetNumberOfTuples(nParticles);

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  float* pData = data->GetPointer(0);
#else
  double* pData = data->GetPointer(0);
#endif

  using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  for (int level = 0; level < this->m_particles->numLevels(); ++level)
  {
    for (MyParIter pti(*this->m_particles, level); pti.isValid(); ++pti)
    {
      auto& particle_attributes = pti.GetStructOfArrays();
      auto& aos = pti.GetArrayOfStructs();

      auto numReal = pti.numParticles();
      // shuffle from the AMReX component order
      for (int j = 0; j < nComps; ++j)
      {
        int compInd = indices[j];
        // component to copy
        const auto &realData = particle_attributes.GetRealData(compInd);

        for (long long i = 0; i < numReal; ++i)
        {
          const auto &part = aos[i];
          if (part.id() > 0)
          {
            pData[i*nComps + j] = realData[i];
          }
        }
      }
      pData += numReal * nComps;
    }
  }

  // get the block of the mesh
  int rank = 0;
  MPI_Comm_rank(this->GetCommunicator(), &rank);

  auto blocks = dynamic_cast<svtkMultiBlockDataSet*>(mesh);

  auto block = dynamic_cast<svtkPolyData*>(blocks->GetBlock(rank));
  block->GetPointData()->AddArray(data);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesSOAIntArray(
  const std::string &arrayName,
  svtkDataObject* mesh)
{
  // get the particles from the particle container
  auto nptsOnProc = this->m_particles->TotalNumberOfParticles(true, true);

  // check that the name of the arrays is listed in the m_intArrays
  IntDataMapType::iterator ait = m_intArrays.find(arrayName);
  if (ait == m_intArrays.end())
  {
    SENSEI_ERROR("No int SOA named \"" << arrayName << "\"");
    return -1;
  }

  // get the index of the structs on the particles
  int index = ait->second;

  // check that the indices are within the templated storage spec
  if(index >= NArrayInt)
  {
    SENSEI_ERROR("Index out of bounds for int SOA named \"" << arrayName << "\"");
    return -1;
  }

  svtkNew<svtkIntArray> data;
  data->SetName(arrayName.c_str());
  data->SetNumberOfComponents(1);
  data->SetNumberOfValues(nptsOnProc);
  int* pData = data->GetPointer(0);

  // fill array
  using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  for (int level = 0; level< this->m_particles->numLevels(); level++)
  {
    for (MyParIter pti(*this->m_particles, level); pti.isValid(); ++pti)
    {
      auto& particle_attributes = pti.GetStructOfArrays();
      auto& aos = pti.GetArrayOfStructs();

      auto numReal = pti.numParticles();
      // shuffle from the AMReX component order
      // component to copy
      const auto &intData = particle_attributes.GetIntData(index);

      for (long long i = 0; i < numReal; ++i)
      {
        const auto &part = aos[i];
        if (part.id() > 0)
        {
          pData[i] = intData[i];
        }
      }
      pData += numReal;
    }
  }

  // get the block of the mesh
  int rank = 0;
  MPI_Comm_rank(this->GetCommunicator(), &rank);

  auto blocks = dynamic_cast<svtkMultiBlockDataSet*>(mesh);

  auto block = dynamic_cast<svtkPolyData*>(blocks->GetBlock(rank));
  block->GetPointData()->AddArray(data);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesAOSRealArray(
  const std::string &arrayName,
  svtkDataObject* mesh)
{
  // get the particles from the particle container
  const auto& particles = this->m_particles->GetParticles();
  auto nptsOnProc = this->m_particles->TotalNumberOfParticles(true, true);

  // check that the name of the struct data is listed in the m_realStructs
  RealDataMapType::iterator ait = m_realStructs.find(arrayName);
  if (ait == m_realStructs.end())
  {
    SENSEI_ERROR("No real AOS named \"" << arrayName << "\"");
    return -1;
  }

  // get the indices of the struct on the particles
  std::vector<int> indices = ait->second;
  int nComps = indices.size();

  // check that the indices are within the templated storage spec
  for (auto i : indices)
  {
    if (i >= NStructReal)
    {
      SENSEI_ERROR("Index out of bounds for real AOS named \"" << arrayName << "\"");
      return -1;
    }
  }

  // allocate the svtk array
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  svtkNew<svtkFloatArray> data;
#else
  svtkNew<svtkDoubleArray> data;
#endif

  data->SetName(arrayName.c_str());
  data->SetNumberOfComponents(nComps);
  data->SetNumberOfTuples(nptsOnProc);

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  float *pData = data->GetPointer(0);
#else
  double *pData = data->GetPointer(0);
#endif

  // copy the data from each level
  using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  for (int level = 0; level<particles.size();level++)
  {
    for (MyParIter pti(*this->m_particles, level); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();

      auto numReal = pti.numParticles();
      // shuffle from the AMReX component order
      for (int j = 0; j < nComps; ++j)
      {
        for (long long i = 0; i < numReal; ++i)
        {
          const auto &part = aos[i];
          if (part.id() > 0)
          {
            pData[i*nComps + j] = part.rdata(indices[j]);
          }
        }
      }
      pData += numReal * nComps;
    }
  }

  // get the block of the mesh
  int rank = 0;
  MPI_Comm_rank(this->GetCommunicator(), &rank);

  auto blocks = dynamic_cast<svtkMultiBlockDataSet*>(mesh);

  auto block = dynamic_cast<svtkPolyData*>(blocks->GetBlock(rank));
  block->GetPointData()->AddArray(data);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesAOSIntArray(
  const std::string &arrayName,
  svtkDataObject* mesh)
{
  // get the particles from the particle container
  const auto& particles = this->m_particles->GetParticles();

  auto nptsOnProc = this->m_particles->TotalNumberOfParticles(true, true);

  // check that the name of the struct data is listed in the m_intStructs
  IntDataMapType::iterator ait = m_intStructs.find(arrayName);
  if (ait == m_intStructs.end())
  {
    SENSEI_ERROR("No int AOS named \"" << arrayName << "\"");
    return -1;
  }

  // get the index of the struct on the particles
  int index = ait->second;

  // check that the index are within the templated storage spec
  if(index >= NStructInt)
  {
    SENSEI_ERROR("Index out of bounds for int AOS named \"" << arrayName << "\"");
    return -1;
  }

  // allocate svtkArray
  svtkNew<svtkIntArray> data;
  data->SetName(arrayName.c_str());
  data->SetNumberOfComponents(1);
  data->SetNumberOfValues(nptsOnProc);
  int* pData = data->GetPointer(0);

  using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  for (int level = 0; level<particles.size(); ++level)
  {
    for (MyParIter pti(*this->m_particles, level); pti.isValid(); ++pti)
    {
      const auto& aos = pti.GetArrayOfStructs();

      long long numReal = pti.numParticles();
      for (long long i = 0; i < numReal; ++i)
      {
        const auto &part = aos[i];
        if (part.id() > 0)
        {
          pData[i] = part.idata(index);
        }
      }
      pData += numReal;
    }
  }

  // get the block of the mesh
  int rank = 0;
  MPI_Comm_rank(this->GetCommunicator(), &rank);

  auto blocks = dynamic_cast<svtkMultiBlockDataSet*>(mesh);

  auto block = dynamic_cast<svtkPolyData*>(blocks->GetBlock(rank));
  block->GetPointData()->AddArray(data);


  return 0;
}

}
