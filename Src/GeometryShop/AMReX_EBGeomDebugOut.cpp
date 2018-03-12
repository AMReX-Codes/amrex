
#include "AMReX_EBGeomDebugOut.H"

namespace amrex
{
  void dumpIFData1(const IFData<1>* a_ifData)
  {
    a_ifData->print(pout());
  }

  void dumpIFData2(const IFData<2>* a_ifData)
  {
    a_ifData->print(pout());
  }

  void dumpIFData3(const IFData<3>* a_ifData)
  {
    a_ifData->print(pout());
  }

  void dumpCCM1(const CutCellMoments<1>* a_ccm)
  {
    a_ccm->print(pout());
  }

  void dumpCCM2(const CutCellMoments<2>* a_ccm)
  {
    a_ccm->print(pout());
  }

  void dumpCCM3(const CutCellMoments<3>* a_ccm)
  {
    a_ccm->print(pout());
  }
}

