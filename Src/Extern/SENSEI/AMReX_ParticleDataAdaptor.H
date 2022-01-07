#ifndef AMReX_ParticleDataAdaptor_h
#define AMReX_ParticleDataAdaptor_h

#include <AMReX_Config.H>
#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFab.H>

#include <DataAdaptor.h>
class vtkPolyData;

namespace amrex
{

template<int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
class ParticleDataAdaptor : public sensei::DataAdaptor
{
public:
  static ParticleDataAdaptor* New();
  senseiTypeMacro(ParticleDataAdaptor, sensei::DataAdaptor);

  // set the source of the simulation data
  int SetDataSource(
    amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> * particles,
    const std::map<std::string, std::vector<int>> & rStructs = {},
    const std::map<std::string, int> & iStructs = {},
    const std::map<std::string, std::vector<int>> & rArrays = {},
    const std::map<std::string, int> & iArrays = {}
  );

  // set the data array names maps internally
  int SetArrayNames(
    const std::map<std::string, std::vector<int>> & rStructs,
    const std::map<std::string, int> & iStructs,
    const std::map<std::string, std::vector<int>> & rArrays,
    const std::map<std::string, int> & iArrays);

  // force mesh origin to 0,0,0. off by default
  void SetPinMesh(int val);

  // get particle id numbers
  int AddParticlesIDArray(vtkDataObject* mesh);

  // get particle cpu numbers (process each particle was generated on)
  int AddParticlesCPUArray(vtkDataObject* mesh);

  // get particle integer arrays in Structs of Arrays format
  int AddParticlesSOAIntArray(const std::string &arrayName, vtkDataObject* mesh);

  // get particle real arrays in Structs of Arrays format
  int AddParticlesSOARealArray(const std::string &arrayName, vtkDataObject* mesh);

  // get particle integer arrays in Array Of Structs format
  int AddParticlesAOSIntArray(const std::string &arrayName, vtkDataObject* mesh);

  // get particle real arrays in Array Of Structs format
  int AddParticlesAOSRealArray(const std::string &arrayName, vtkDataObject* mesh);

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
  int GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) override;
  int AddGhostNodesArray(vtkDataObject* mesh, const std::string &meshName) override;
  int AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) override;
  int AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;
  int ReleaseData() override;

protected:
  ParticleDataAdaptor() = default;
  ~ParticleDataAdaptor() = default;

private:
  vtkPolyData* BuildParticles();

  const std::string m_particlesName = "particles";

  amrex::ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>* m_particles;

  using RealDataMapType = std::map<std::string, std::vector<int>>;
  using IntDataMapType = std::map<std::string, int>;

  RealDataMapType m_realStructs;
  IntDataMapType m_intStructs;
  RealDataMapType m_realArrays;
  IntDataMapType m_intArrays;
};

}

#include "AMReX_ParticleDataAdaptorI.H"

#endif
#endif
