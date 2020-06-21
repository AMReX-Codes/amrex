#ifndef AMReX_AmrDataAdaptor_h
#define AMReX_AmrDataAdaptor_h

#include "DataAdaptor.h"
#include <AMReX_Amr.H>

namespace amrex
{
class AmrDataAdaptor : public sensei::DataAdaptor
{
public:
  static AmrDataAdaptor* New();
  senseiTypeMacro(AmrDataAdaptor, sensei::DataAdaptor);

  // set the source of the simulation data
  int SetDataSource(amrex::Amr *amr);

  // force mesh origin to 0,0,0. Off by default
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
  AmrDataAdaptor();
  ~AmrDataAdaptor();

private:
  struct InternalsType;
  InternalsType *Internals;
};

}

#endif
