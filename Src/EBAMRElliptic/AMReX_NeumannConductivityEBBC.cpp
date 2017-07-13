 
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */
#include "AMReX_NeumannConductivityEBBC.H"

namespace amrex
{
  void NeumannConductivityEBBC::applyEBFlux(EBCellFAB              &       a_lphi,
                                            const EBCellFAB        &       a_phi,
                                            const vector<VolIndex> &       a_vofsToChange,
                                            const MFIter           &       a_mfi,
                                            const Real             &       a_factor,
                                            const bool             &       a_useHomogeneous)
  {
    CH_TIME("NeumannConductivityEBBC::applyEBFlux");
    CH_assert(a_lphi.nComp() == 1 );
    CH_assert(a_phi.nComp()  == 1);

    Real flux = 0.0;

    const EBISBox&   ebisBox = a_phi.getEBISBox();
    for(int ivof = 0; ivof < a_vofsToChange.size(); ivof++)
    {
      const VolIndex& vof = a_vofsToChange[ivof];
      const RealVect& centroid = ebisBox.bndryCentroid(vof);
      centroid *= m_dx;
      centroid += m_probLo;
      Real point = EBArith::getVoFLocation(vof, m_dx, centroid);
      Real value = bcvaluefunc(point);
      flux = value;

      const Real& areaFrac = ebisBox.bndryArea(vof);
      flux *= areaFrac;
      Real bcoef = (*m_bcoe)[a_mfi].getEBFlux()(vof,0);
      flux *= m_beta*bcoef;
      a_lphi(vof,0) += flux * a_factor;
    }
  }
}

