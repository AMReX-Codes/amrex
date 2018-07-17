#include "AMReX_AllRegularService.H"

namespace amrex
{
  /*******************/
  /*******************/
  AllRegularService::AllRegularService()
  {
  }
  /*******************/
  /*******************/
  AllRegularService::~AllRegularService()
  {
  }
           
  /*******************/
  /*******************/
  bool
  AllRegularService::
  isRegular(const Box& a_region,
            const Box& a_domain,
            const RealVect& a_origin,
            const Real& a_dx) const
  {
    return true;
  }
           
  /*******************/
  /*******************/
  bool
  AllRegularService::isCovered(const Box& a_region,
                               const Box& a_domain,
                               const RealVect& a_origin,
                               const Real& a_dx) const
  {
    return false;
  }
           
  /*******************/
  /*******************/
  void
  AllRegularService::fillGraph(BaseFab<int>&        a_regIrregCovered,
                               Vector<IrregNode>&   a_nodes,
                               NodeMap&             a_intersections,
                               const Box&           a_validRegion,
                               const Box&           a_ghostRegion,
                               const Box& a_domain,
                               const RealVect&      a_origin,
                               const Real&          a_dx) const
  {
    a_nodes.resize(0);
    a_regIrregCovered.resize(a_ghostRegion, 1);
    //set all cells to regular
    a_regIrregCovered.setVal(1);
    a_intersections.clear();
  }
}
