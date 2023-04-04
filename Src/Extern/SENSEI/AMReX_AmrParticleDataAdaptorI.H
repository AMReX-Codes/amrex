#include "AMReX_AmrParticleDataAdaptor.H"

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

//-----------------------------------------------------------------------------
template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetDataSource(
  Amr *amr,
  amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
  const std::map<std::string, std::vector<int>> & rStructs,
  const std::map<std::string, int> & iStructs,
  const std::map<std::string, std::vector<int>> & rArrays,
  const std::map<std::string, int> & iArrays)
{
  int retMesh = this->m_meshAdaptor->SetDataSource(amr);
  int retPtls = this->m_particleAdaptor->SetDataSource(particles, rStructs, iStructs, rArrays, iArrays);
  return retMesh + retPtls;
}

// SENSEI API
#if SENSEI_VERSION_MAJOR >= 3
template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshMetadata(
  unsigned int id,
  sensei::MeshMetadataPtr &metadata)
{
  if(id == 0)
  {
    return this->m_meshAdaptor->GetMeshMetadata(0, metadata);
  }
  if(id == 1)
  {
    return this->m_particleAdaptor->GetMeshMetadata(1, metadata);
  }
  SENSEI_ERROR("unknown mesh ID in GetMeshMetadata");
  return -1;
}

#else
template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshName(
  unsigned int id,
  std::string &meshName)
{
  if(id == 0)
  {
    return this->m_meshAdaptor->GetMeshName(0, meshName);
  }
  if(id == 1)
  {
    return this->m_particleAdaptor->GetMeshName(1, meshName);
  }
  SENSEI_ERROR("unknown mesh ID in GetMeshName");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshHasGhostNodes(
  const std::string &meshName, int &nLayers)
{
  if(id == 0)
  {
    return this->m_meshAdaptor->GetMeshHasGhostNodes(meshName, nLayers);
  }
  if(id == 1)
  {
    return this->m_particleAdaptor->GetMeshHasGhostNodes(meshName, nLayers);
  }
  SENSEI_ERROR("unknown mesh ID in GetMeshHasGhostNodes");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMeshHasGhostCells(
  const std::string &meshName, int &nLayers)
{
  if(id == 0)
  {
    return this->m_meshAdaptor->GetMeshHasGhostCells(meshName, nLayers);
  }
  if(id == 1)
  {
    return this->m_particleAdaptor->GetMeshHasGhostCells(meshName, nLayers);
  }
  SENSEI_ERROR("unknown mesh ID in GetMeshHasGhostCells");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetNumberOfArrays(
  const std::string &meshName,
  int association,
  unsigned int &numberOfArrays)
{
  if(id == 0)
  {
    return this->m_meshAdaptor->GetNumberOfArrays(meshName, association, numberOfArrays);
  }
  if(id == 1)
  {
    return this->m_particleAdaptor->GetNumberOfArrays(meshName, association, numberOfArrays);
  }
  SENSEI_ERROR("unknown mesh ID in GetNumberOfArrays");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetArrayName(
  const std::string &meshName,
  int association,
  unsigned int index,
  std::string &arrayName)
{
  if(id == 0)
  {
    return this->m_meshAdaptor->GetArrayName(meshName, association, index, arrayName);
  }
  if(id == 1)
  {
    return this->m_particleAdaptor->GetArrayName(meshName, association, index, arrayName);
  }
  SENSEI_ERROR("unknown mesh ID in GetArrayName");
  return -1;
}
#endif

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetNumberOfMeshes(
  unsigned int &numMeshes)
{
  numMeshes = 2;
  return 0;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetMesh(
  const std::string &meshName,
  bool structureOnly,
  svtkDataObject *&mesh)
{
  if(meshName == m_meshName)
  {
    return this->m_meshAdaptor->GetMesh(meshName, structureOnly, mesh);
  }
  if(meshName == m_particlesName)
  {
    return this->m_particleAdaptor->GetMesh(meshName, structureOnly, mesh);
  }
  SENSEI_ERROR("unknown mesh name in GetMesh");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddGhostNodesArray(
  svtkDataObject* mesh,
  const std::string &meshName)
{
  if(meshName == m_meshName)
  {
    return this->m_meshAdaptor->AddGhostNodesArray(mesh, meshName);
  }
  if(meshName == m_particlesName)
  {
    return this->m_particleAdaptor->AddGhostNodesArray(mesh, meshName);
  }
  SENSEI_ERROR("unknown mesh name in AddGhostNodesArray");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddGhostCellsArray(
  svtkDataObject* mesh,
  const std::string &meshName)
{
  if(meshName == m_meshName)
  {
    return this->m_meshAdaptor->AddGhostCellsArray(mesh, meshName);
  }
  if(meshName == m_particlesName)
  {
    return this->m_particleAdaptor->AddGhostCellsArray(mesh, meshName);
  }
  SENSEI_ERROR("unknown mesh name in AddGhostCellsArray");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::AddArray(
  svtkDataObject* mesh,
  const std::string &meshName,
  int association,
  const std::string &arrayName)
{
  if(meshName == m_meshName)
  {
    return this->m_meshAdaptor->AddArray(mesh, meshName, association, arrayName);
  }
  if(meshName == m_particlesName)
  {
    return this->m_particleAdaptor->AddArray(mesh, meshName, association, arrayName);
  }
  SENSEI_ERROR("unknown mesh name in AddArray");
  return -1;
}

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
int AmrParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::ReleaseData()
{
  this->m_particleAdaptor->ReleaseData();
  this->m_meshAdaptor->ReleaseData();
  return 0;
}

} // end namespace amrex
