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

#include "AMReX_DirichletConductivityEBBC.H"
namespace amrex
{
  void
  DirichletConductivityEBBC::
  defineStencils()
  {
    m_fluxStencil.define(m_eblg.getDBL());
    m_fluxWeight .define(m_eblg.getDBL());
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {

      const Box     & grid = m_eblg.getDBL()  [mfi];
      const EBISBox & ebis = m_eblg.getEBISL()[mfi];
      IntVectSet ivsIrreg  = ebis.getIrregIVS(grid);

      m_fluxStencil[mfi].define(ivsIrreg, ebis.getEBGraph(), 1);
      m_fluxWeight [mfi].define(ivsIrreg, ebis.getEBGraph(), 1);

      const IntVectSet& cfivs = (*m_eblg.getCFIVS())[mfi];
      for (VoFIterator vofit(ivs, ebis.getEBGraph()); vofit.ok(); ++vofit)
      {
        VoFStencil gradStencil;
        Real curWeight;
        if (m_order == 1)
        {
          getFirstOrderStencil( gradStencil,curWeight,vof,ebis,m_dx);
        }
        else 
        {
          getSecondOrderStencil(gradStencil,curWeight,vof,ebis,m_dx,cfivs);
        }
        
        VoFStencil fluxStencilPt = gradStencil;
        m_fluxStencil[dit()](vofit(), 0) = poissSten[dit()](vofit(), 0);
        Real factor = (*m_bcoe)[dit()](vofit(), 0);
        factor *= m_beta;
        fluxStencilPt *= factor;

        m_fluxStencil[dit()](vofit(), 0) = fluxStencilPt;
        m_fluxWeight[mfi](vofit(), 0)    = curWeight;
      }
    }
  }

  /*****************/
  void
  DirichletConductivityEBBC::
  applyEBFlux(EBCellFAB              &       a_lphi,
              const EBCellFAB        &       a_phi,
              const vector<VolIndex> &       a_vofsToChange,
              const MFIter           &       a_mfi,
              const Real             &       a_factor,
              const bool             &       a_useHomogeneous)
  {
    BL_PROFILE("DirichletConductivityEBBC::applyEBFlux");
    BL_ASSERT(a_lphi.nComp() == 1 );
    BL_ASSERT( a_phi.nComp() == 1);

    const BaseIVFAB<Real>& poissWeight = m_fluxWeight[a_mfi];
    Real value = 0.0;

    const EBISBox&   ebisBox = a_phi.getEBISBox();
    for(int ivof = 0; ivof < a_vofsToChange.size(); ivof++)
    {
      const VolIndex& vof = a_vofsToChange[ivof];


      const RealVect& centroid = ebisBox.bndryCentroid(vof);
      centroid *= m_dx;
      centroid += m_probLo;
      Real point = EBArith::getVoFLocation(vof, m_dx, centroid);
      Real value = bcvaluefunc(point);

      Real poissWeightPt = poissWeight(vof, 0);
      const Real& areaFrac = ebisBox.bndryArea(vof);
      const Real& bcoef    = (*m_bcoe)[a_dit](vof,0);
      Real flux = poissWeightPt*value*areaFrac;
      Real compFactor = a_factor*bcoef*m_beta;
      a_lphi(vof,0) += flux * compFactor;
    }
  }


  bool
  DirichletConductivityEBBC::
  getSecondOrderStencil(VoFStencil&          a_stencil,
                        Real&                a_weight,
                        vector<VoFStencil>&  a_pointStencils,
                        vector<Real>&        a_distanceAlongLine,
                        const VolIndex&      a_vof,
                        const EBISBox&       a_ebisBox,
                        const RealVect&      a_dx,
                        const IntVectSet&    a_cfivs)
  {
    BL_PROFILE("DirichletPoissonEBBC::getSecondOrderStencil1");

    a_stencil.clear();
    bool dropOrder = false;
    EBArith::johanStencil(dropOrder, a_pointStencils, a_distanceAlongLine,
                          a_vof, a_ebisBox, a_dx, a_cfivs);
    if (dropOrder)
    {
      return true;
    }

    //if we got this far, sizes should be at least 2
    CH_assert(a_distanceAlongLine.size() >= 2);
    CH_assert(a_pointStencils.size() >= 2);
    Real x1 = a_distanceAlongLine[0];
    Real x2 = a_distanceAlongLine[1];
    //fit quadratic function to line and find the gradient at the origin
    //grad = (x2*x2*(phi1-phi0)-x1*x1(phi2-phi0))/(x2*x2*x1 - x1*x1*x2);
    Real denom = x2*x2*x1 - x1*x1*x2;
    //not done by reference because i want point stencils checkable externally.
    VoFStencil phi1Sten = a_pointStencils[0];
    VoFStencil phi2Sten = a_pointStencils[1];
    phi1Sten *=-x2*x2/denom;
    phi2Sten *= x1*x1/denom;
    //weight is the multiplier of the inhomogeneous value (phi0)
    a_weight =-x1*x1/denom + x2*x2/denom;
    a_stencil += phi1Sten;
    a_stencil += phi2Sten;
    //if we got this far, we have a second order stencil;
    return false;
  }

  void 
  DirichletConductivityEBBC::
  getFirstOrderStencil(VoFStencil&     a_stencil,
                       Real&           a_weight,
                       const VolIndex& a_vof,
                       const EBISBox&  a_ebisBox,
                       const RealVect& a_dx)
  {
    BL_PROFILE("DirichletPoissonEBBC::getFirstOrderStencil1");
    RealVect centroid = a_ebisBox.centroid(a_vof);
    EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisBox, a_dx, m_domain, 0);
  }

  void 
  DirichletConductivityEBBC::
  getSecondOrderStencil(VoFStencil&       a_stencil,
                        Real&             a_weight,
                        const VolIndex&   a_vof,
                        const EBISBox&    a_ebisBox,
                        const RealVect&   a_dx,
                        const IntVectSet& a_cfivs)
  {
    BL_PROFILE("DirichletPoissonEBBC::getSecondOrderStencil2");
    vector<VoFStencil>  pointStencils;
    vector<Real>        distanceAlongLine;
    bool needToDropOrder = getSecondOrderStencil(a_stencil, a_weight,
                                                 pointStencils, distanceAlongLine,
                                                 a_vof, a_ebisBox, a_dx, a_cfivs);
    if (needToDropOrder)
    {
      getFirstOrderStencil(a_stencil,a_weight,a_vof,a_ebisBox,a_dx);
    }
  }
}
