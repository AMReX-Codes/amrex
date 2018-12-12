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
  int GetMeshName(unsigned int id, std::string &meshName) override;
  int GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) override;
  int GetMeshHasGhostNodes(const std::string &meshName, int &nLayers) override;
  int AddGhostNodesArray(vtkDataObject* mesh, const std::string &meshName) override;
  int GetMeshHasGhostCells(const std::string &meshName, int &nLayers) override;
  int AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) override;
  int AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;
  int GetNumberOfArrays(const std::string &meshName, int association, unsigned int &numberOfArrays) override;
  int GetArrayName(const std::string &meshName, int association, unsigned int index, std::string &arrayName) override;
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
