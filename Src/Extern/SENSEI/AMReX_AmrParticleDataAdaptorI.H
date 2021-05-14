#include "Profiler.h"
#include "Error.h"
#include "VTKUtils.h"
#include "MeshMetadata.h"
// vtk includes
#include <vtkPolyData.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkMultiBlockDataSet.h>



// local includes
#include "AMReX_AmrParticleDataAdaptor.H"
#include <AMReX_InSituUtils.H>

namespace amrex
{
//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>*
AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::New()
{
  auto result = new AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>;
  result->InitializeObjectBase();
  return result;
}
// senseiNewMacro(AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>);

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetDataSource(
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> *tracers,
    std::map<std::string, std::vector<int>> rStructs,
    std::map<std::string, int> iStructs,
    std::map<std::string, std::vector<int>> rArrays,
    std::map<std::string, int> iArrays)
{
  int ret = 0;
  // set the data source on the tracers
  this->ReleaseData();
  this->m_tracers = tracers;

  if(ret)
  {
    SENSEI_ERROR("problem with tracers in SetDataSource");
    return ret;
  }

  // set the array names
  ret = this->SetArrayNames(rStructs, iStructs, rArrays, iArrays);
  if(ret)
  {
    SENSEI_ERROR("problem with array names in SetDataSource");
    return ret;
  }
  return ret;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetArrayNames(
  std::map<std::string, std::vector<int>> rStructs,
  std::map<std::string, int> iStructs,
  std::map<std::string, std::vector<int>> rArrays,
  std::map<std::string, int> iArrays)
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetNumberOfMeshes(unsigned int &numMeshes)
{
  numMeshes = 1;
  return 0;
}

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshName(
  unsigned int id, std::string &meshName)
{
  meshName = m_particlesName;
  return 0;
}
#endif

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR < 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshHasGhostNodes(
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshHasGhostCells(
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetNumberOfArrays(
  const std::string &meshName,
  int association,
  unsigned int &numberOfArrays)
{
  numberOfArrays = 0;
  if(association == vtkDataObject::POINT)
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetArrayName(
  const std::string &meshName,
  int association,
  unsigned int index,
  std::string &arrayName)
{
  if(association == vtkDataObject::POINT)
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMesh(
  const std::string &meshName,
  bool structureOnly,
  vtkDataObject *&mesh)
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
  vtkMultiBlockDataSet *mb = vtkMultiBlockDataSet::New();

  if(structureOnly)
  {
    mesh = mb;
    return 0;
  }

  mb->SetNumberOfBlocks(nprocs);
  vtkPolyData *pd = BuildParticles();
  mb->SetBlock(rank, pd);
  pd->Delete();
  mesh = mb;
  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddGhostNodesArray(
  vtkDataObject*,
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddGhostCellsArray(
  vtkDataObject*,
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddArray(
  vtkDataObject* mesh,
  const std::string &meshName,
  int association,
  const std::string &arrayName)
{
  if (meshName != m_particlesName)
  {
    SENSEI_ERROR("no mesh named \"" << meshName << "\"");
    return -1;
  }

  if (association != vtkDataObject::POINT)
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
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::ReleaseData()
{
  this->m_tracers = nullptr;
  return 0;
}

//-----------------------------------------------------------------------------
#if SENSEI_VERSION_MAJOR >= 3
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshMetadata(
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
  metadata->MeshType = VTK_MULTIBLOCK_DATA_SET;

  // block mesh type (all)
  metadata->BlockType = VTK_POLY_DATA;

  // global number of blocks (all)
  metadata->NumBlocks = nprocs;

  // number of blocks on each rank (all)
  metadata->NumBlocksLocal = {1};

  // global cell index space extent (Cartesian, AMR, optional)
  // std::array<int,6> Extent;

  // global bounding box (all, optional)
  // std::array<double,6> Bounds;

  // type enum of point data (unstructured, optional)
  // int CoordinateType;

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
  metadata->NumLevels;

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
    metadata->ArrayCentering.push_back(vtkDataObject::POINT);
  }
  for(auto s : m_intStructs)
  {
    metadata->ArrayCentering.push_back(vtkDataObject::POINT);
  }
  for(auto s : m_realArrays)
  {
    metadata->ArrayCentering.push_back(vtkDataObject::POINT);
  }
  for(auto s : m_intArrays)
  {
    metadata->ArrayCentering.push_back(vtkDataObject::POINT);
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
    metadata->ArrayType.push_back(VTK_DOUBLE);
  }
  for(auto s : m_intStructs)
  {
    metadata->ArrayType.push_back(VTK_INT);
  }
  for(auto s : m_realArrays)
  {
    metadata->ArrayType.push_back(VTK_DOUBLE);
  }
  for(auto s : m_intArrays)
  {
    metadata->ArrayType.push_back(VTK_INT);
  }

  // global min,max of each array (all, optional)
  // metadata->ArrayRange = {};

  // rank where each block resides (all, optional)
  // metadata->BlockOwner.resize(nprocs);
  metadata->BlockOwner = {rank};
  // std::iota(metadata->BlockOwner.begin(), metadata->BlockOwner.end(), 0);

  // global id of each block (all, optional)
  // TO DO: this is lazy
  metadata->BlockIds = metadata->BlockOwner;

  // note: for AMR BlockNumPoints and BlockNumCells are always global
  // metadata->BlockNumPoints.resize(metadata->NumBlocks);
  // metadata->BlockNumCells.resize(metadata->NumBlocks);
  // metadata->BlockCellArraySize.resize(metadata->NumBlocks);

  // number of points for each block (all, optional)
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);
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

  // min max of each array on each block.
  // indexed by block then array. (all, optional)
  // std::vector<std::vector<std::array<double,2>>> BlockArrayRange;

  // refinement ratio in i,j, and k directions for each level (AMR)
  // std::vector<std::array<int,3>> RefRatio;

  // number of blocks in each level (AMR)
  // std::vector<int> BlocksPerLevel;

  // AMR level of each block (AMR)
  // std::vector<int> BlockLevel;

  // flag indicating presence of a periodic boundary in the i,j,k direction (all)
  // std::array<int,3> PeriodicBoundary;

  return 0;
}
#endif
//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
vtkPolyData* AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::BuildParticles()
{
  // return particle data pd
  vtkPolyData* pd  = vtkPolyData::New();

  const auto& particles = this->m_tracers->GetParticles();
  const int numParticles = this->m_tracers->TotalNumberOfParticles(true, true);

  // allocate vertex storage for particles
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numParticles);

  // allocate connectivity array for particles
  vtkSmartPointer<vtkCellArray> vertex = vtkSmartPointer<vtkCellArray>::New();
  vertex->AllocateExact(numParticles, 1);

  int count = 0;
  // loop over levels in the particle container
  for (int lev = 0; lev < particles.size();  lev++)
  {
    if(count > numParticles)
    {
      SENSEI_ERROR("too many points in BuildParticles!");
      return nullptr;
    }

    // loop over ParticleTiles on level
    auto& pmap = particles[lev];
    for (const auto& kv : pmap)
    {
      // loop over particles in ParticleTile
      const auto& aos = kv.second.GetArrayOfStructs();
      for (auto it = aos.cbegin(); it != aos.cend(); ++it)
      {
        if (it->id() > 0)
        {
          vertex->InsertNextCell(1);
          vertex->InsertCellPoint(count);
          if (AMREX_SPACEDIM == 2)
            points->SetPoint(count, it->pos(0), it->pos(1), 0.0);
          if (AMREX_SPACEDIM == 3)
            points->SetPoint(count, it->pos(0), it->pos(1), it->pos(2));
          count++;
        }
      }
    }
  }

  // add point and vertex data to output mesh
  pd->SetPoints(points);
  pd->SetVerts(vertex);
  return pd;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesIDArray(vtkDataObject* mesh)
{
  auto tracers = dynamic_cast<vtkPolyData*>(mesh);
  const auto& particles = this->m_tracers->GetParticles();
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);

  vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
  intArray->SetName("id");
  intArray->SetNumberOfComponents(1);
  intArray->SetNumberOfValues(nptsOnProc);
  int count = 0;
  for (int level = 0; level < particles.size(); level++)
  {
    auto &pmap = particles[level];
    for (const auto &kv : pmap)
    {
      const auto &aos = kv.second.GetArrayOfStructs();
      for (auto it = aos.cbegin();it != aos.cend();++it)
      {
        if (it->id() > 0)
        {
          intArray->SetValue(count, it->id());
          count++;
        }
      }
    }
  }

  // the association for this array is vtkDataObject::POINT
  tracers->GetPointData()->AddArray(intArray);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesCPUArray(
  vtkDataObject* mesh)
{
  auto tracers = dynamic_cast<vtkPolyData*>(mesh);
  const auto& particles = this->m_tracers->GetParticles();
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);

  vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
  intArray->SetName("cpu");
  intArray->SetNumberOfComponents(1);
  intArray->SetNumberOfValues(nptsOnProc);
  int count = 0;
  for (int level = 0; level < particles.size() ;level++)
  {
    auto &pmap = particles[level];
    for (const auto &kv : pmap)
    {
      const auto &aos = kv.second.GetArrayOfStructs();
      for (auto it = aos.cbegin();it != aos.cend();++it)
      {
        if (it->cpu() > 0)
        {
          intArray->SetValue(count, it->cpu());
          count++;
        }
      }
    }
  }

  // the association for this array is vtkDataObject::POINT
  tracers->GetPointData()->AddArray(intArray);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesSOARealArray(
  const std::string &arrayName,
  vtkDataObject* mesh)
{
  // get the particles from the particle container
  const auto& particles = this->m_tracers->GetParticles();
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);

  // check that the name of the arrays is listed in the m_realArrays
  if (m_realArrays.find(arrayName) == m_realArrays.end())
  {
    SENSEI_ERROR("No real SOA named \"" << arrayName << "\"");
    return -1;
  }

  // get the indices of the structs on the particles
  std::vector<int> indices = this->m_realArrays[arrayName];
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

  int rank = 0;
  MPI_Comm_rank(this->GetCommunicator(), &rank);
  // allocate the vtkArray
  auto blocks = dynamic_cast<vtkMultiBlockDataSet*>(mesh);
  auto tracers = dynamic_cast<vtkPolyData*>(blocks->GetBlock(rank));

  vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
  doubleArray->SetName(arrayName.c_str());
  doubleArray->SetNumberOfComponents(nComps);
  doubleArray->SetNumberOfTuples(nptsOnProc);

  // fill array
  std::vector<double> tuple(nComps);
  int count = 0;
  for(int level = 0; level < particles.size(); level++)
  {
    auto &pmap = particles[level];
    for(const auto &kv : pmap)
    {
      const auto &aos = kv.second.GetArrayOfStructs();
      const auto &soa = kv.second.GetStructOfArrays();
      int element = 0;
      for(auto it = aos.cbegin(); it != aos.cend(); ++it)
      {
        if(it->id() > 0)
        {
          for(int i = 0; i < nComps; ++i)
          {
            tuple[i] = soa.GetRealData(indices[i])[element];
          }

          doubleArray->SetTuple(count, tuple.data());
          element++;
          count++;
        }
      }
    }
  }
  tracers->GetPointData()->AddArray(doubleArray);

  return 0;
}

  //-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesSOAIntArray(
  const std::string &arrayName,
  vtkDataObject* mesh)
{
  // get the particles from the particle container
  const auto& particles = this->m_tracers->GetParticles();
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);

  // check that the name of the arrays is listed in the m_intArrays
  if (m_intArrays.find(arrayName) == m_intArrays.end())
  {
    SENSEI_ERROR("No int SOA named \"" << arrayName << "\"");
    return -1;
  }

  // get the index of the structs on the particles
  int index = this->m_intArrays[arrayName];

  // no multicomponent integer arrays
  int nComps = 1;

  // check that the indices are within the templated storage spec
  if(index >= NArrayInt)
  {
    SENSEI_ERROR("Index out of bounds for int SOA named \"" << arrayName << "\"");
    return -1;
  }

  // allocate vtkArray
  auto tracers = dynamic_cast<vtkPolyData*>(mesh);
  vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
  intArray->SetName(arrayName.c_str());
  intArray->SetNumberOfComponents(nComps);
  intArray->SetNumberOfValues(nptsOnProc);

  // fill array
  std::vector<int> tuple(nComps);
  int count = 0;
  for (int level = 0; level<particles.size();level++)
  {
    auto &pmap = particles[level];
    for (const auto &kv : pmap)
    {
      const auto &aos = kv.second.GetArrayOfStructs();
      const auto &soa = kv.second.GetStructOfArrays();
      int element = 0;
      for (auto it = aos.cbegin();it != aos.cend();++it)
      {
        if (it->id() > 0)
        {
          intArray->SetValue(count, soa.GetIntData(index)[element]);
          element++;
          count++;
        }
      }
    }
  }
  tracers->GetPointData()->AddArray(intArray);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesAOSRealArray(
  const std::string &arrayName,
  vtkDataObject* mesh)
{
  // get the particles from the particle container
  const auto& particles = this->m_tracers->GetParticles();
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);

  // check that the name of the struct data is listed in the m_realStructs
  if (m_realStructs.find(arrayName) == m_realStructs.end())
  {
    SENSEI_ERROR("No real AOS named \"" << arrayName << "\"");
    return -1;
  }

  // get the indices of the struct on the particles
  std::vector<int> indices = this->m_realStructs[arrayName];
  int nComps = indices.size();

  // check that the indices are within the templated storage spec
  for(auto i : indices)
  {
    if(i >= NStructReal)
    {
      SENSEI_ERROR("Index out of bounds for real AOS named \"" << arrayName << "\"");
      return -1;
    }
  }

   int rank = 0;
  MPI_Comm_rank(this->GetCommunicator(), &rank);
  // allocate the vtkArray
  auto blocks = dynamic_cast<vtkMultiBlockDataSet*>(mesh);
  auto tracers = dynamic_cast<vtkPolyData*>(blocks->GetBlock(rank));

  vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
  doubleArray->SetName(arrayName.c_str());
  doubleArray->SetNumberOfComponents(nComps);
  doubleArray->SetNumberOfTuples(nptsOnProc);

  // fill array
  int count = 0;
  std::vector<double> tuple(nComps);
  for (int level = 0; level<particles.size();level++)
  {
    auto &pmap = particles[level];
    for (const auto &kv : pmap)
    {
      const auto &aos = kv.second.GetArrayOfStructs();
      for (auto it = aos.cbegin();it != aos.cend();++it)
      {
        if (it->id() > 0)
        {
          for(int i = 0; i < nComps; ++i)
          {
            tuple[i] = it->rdata(indices[i]);
          }
          doubleArray->SetTuple(count, tuple.data());
          count++;
        }
      }
    }
  }
  tracers->GetPointData()->AddArray(doubleArray);

  return 0;
}

//-----------------------------------------------------------------------------
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddParticlesAOSIntArray(
  const std::string &arrayName,
  vtkDataObject* mesh
)
{
  // get the particles from the particle container
  const auto& particles = this->m_tracers->GetParticles();
  auto nptsOnProc = this->m_tracers->TotalNumberOfParticles(true, true);

  // check that the name of the struct data is listed in the m_intStructs
  if (m_intStructs.find(arrayName) == m_intStructs.end())
  {
    SENSEI_ERROR("No int AOS named \"" << arrayName << "\"");
    return -1;
  }

  // get the index of the struct on the particles
  int index = this->m_intStructs[arrayName];

  // no multicomponent integer arrays
  int nComps = 1;

  // check that the index are within the templated storage spec
  if(index >= NStructInt)
  {
    SENSEI_ERROR("Index out of bounds for int AOS named \"" << arrayName << "\"");
    return -1;
  }

  // allocate vtkArray
  auto tracers = dynamic_cast<vtkPolyData*>(mesh);
  vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
  intArray->SetName(arrayName.c_str());
  intArray->SetNumberOfComponents(nComps);
  intArray->SetNumberOfValues(nptsOnProc);

  int count = 0;
  std::vector<int> tuple(nComps);
  for (int level = 0; level<particles.size();level++)
  {
    auto &pmap = particles[level];
    for (const auto &kv : pmap)
    {
      const auto &aos = kv.second.GetArrayOfStructs();
      for (auto it = aos.cbegin();it != aos.cend();++it)
      {
        if (it->id() > 0)
        {
          intArray->SetValue(count, it->idata(index));
          count++;
        }
      }
    }
  }

  tracers->GetPointData()->AddArray(intArray);

  return 0;
}

}
