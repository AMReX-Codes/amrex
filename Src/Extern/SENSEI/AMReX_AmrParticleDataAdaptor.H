#ifndef AMReX_AmrParticleDataAdaptor_h
#define AMReX_AmrParticleDataAdaptor_h

// AMReX includes
#include <AMReX_Config.H>
#ifdef AMREX_PARTICLES

#include <AMReX_AmrDataAdaptor.H>
#include <AMReX_Particles.H>

#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParticleDataAdaptor.H>
// sensei includes
#include "DataAdaptor.h"

namespace amrex
{
template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
class AmrParticleDataAdaptor : public sensei::DataAdaptor
{
public:
  static AmrParticleDataAdaptor* New();
  senseiTypeMacro(AmrParticleDataAdaptor, sensei::DataAdaptor);

  int SetDataSource(
    Amr *amr,
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
    const std::map<std::string, std::vector<int>> & rStructs = {},
    const std::map<std::string, int> & iStructs = {},
    const std::map<std::string, std::vector<int>> & rArrays = {},
    const std::map<std::string, int> & iArrays = {}
  );

  // SENSEI API
#if SENSEI_VERSION_MAJOR >= 3
  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) override;
#else
  int GetMeshName(unsigned int id, std::string &meshName) override;
  int GetMeshHasGhostNodes(const std::string &meshName, int &nLayers) override;
  int GetMeshHasGhostCells(const std::string &meshName, int &nLayers) override;
  int GetNumberOfArrays(const std::string &meshName, int association, unsigned int &numberOfArrays) override;
  int GetArrayName(const std::string &meshName, int association, unsigned int index, std::string &arrayName) override;
#endif
  int GetNumberOfMeshes(unsigned int &numMeshes) override;
  int GetMesh(const std::string &meshName, bool structureOnly, svtkDataObject *&mesh) override;
  int AddGhostNodesArray(svtkDataObject* mesh, const std::string &meshName) override;
  int AddGhostCellsArray(svtkDataObject* mesh, const std::string &meshName) override;
  int AddArray(svtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;
  int ReleaseData() override;

protected:
  AmrParticleDataAdaptor()
  {
    m_meshAdaptor = AmrDataAdaptor::New();
    m_particleAdaptor = ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>::New();
  }

  ~AmrParticleDataAdaptor()
  {
    m_meshAdaptor->Delete();
    m_particleAdaptor->Delete();
  }

private:
  ParticleDataAdaptor<NStructReal, NStructInt, NArrayReal, NArrayInt>* m_particleAdaptor;
  AmrDataAdaptor* m_meshAdaptor;

  const std::string m_meshName = "mesh";
  const std::string m_particlesName = "particles";
};
}

#include "AMReX_AmrParticleDataAdaptorI.H"
#endif // AMREX_PARTICLES
#endif // AMReX_AmrParticleDataAdaptor_h
