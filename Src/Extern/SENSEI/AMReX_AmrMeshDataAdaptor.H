#ifndef AMReX_AmrMeshDataAdaptor_h
#define AMReX_AmrMeshDataAdaptor_h

#include "DataAdaptor.h"
#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFab.H>

namespace amrex
{

class AmrMeshDataAdaptor : public sensei::DataAdaptor
{
public:
  static AmrMeshDataAdaptor* New();
  senseiTypeMacro(AmrMeshDataAdaptor, sensei::DataAdaptor);

  // set the source of the simulation data
  int SetDataSource(amrex::AmrMesh *mesh,
    const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
    const std::vector<std::vector<std::string>> &arrayNames);

  // force mesh origin to 0,0,0. off by default
  void SetPinMesh(int val);

  // SENSEI API
  int GetNumberOfMeshes(unsigned int &numMeshes) override;
#if SENSEI_VERSION_MAJOR >= 3
  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) override;
#else
  int GetMeshName(unsigned int id, std::string &meshName) override;
  int GetMeshHasGhostNodes(const std::string &meshName, int &nLayers) override;
  int GetMeshHasGhostCells(const std::string &meshName, int &nLayers) override;
  int GetNumberOfArrays(const std::string &meshName, int association, unsigned int &numberOfArrays) override;
  int GetArrayName(const std::string &meshName, int association, unsigned int index, std::string &arrayName) override;
#endif
  int GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) override;
  int AddGhostNodesArray(vtkDataObject* mesh, const std::string &meshName) override;
  int AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) override;
  int AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;
  int ReleaseData() override;

protected:
  AmrMeshDataAdaptor();
  ~AmrMeshDataAdaptor();

private:
  struct InternalsType;
  InternalsType *Internals;
};

}

#endif
