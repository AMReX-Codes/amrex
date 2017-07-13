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
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const IntVectSet& ivs = poissSten[dit()].getIVS();
      const EBGraph&    ebg = poissSten[dit()].getEBGraph();
      m_fluxStencil[dit()].define(ivs, ebg, 1);
      for (VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
      {
        m_fluxStencil[dit()](vofit(), 0) = poissSten[dit()](vofit(), 0);
        Real factor = (*m_bcoe)[dit()](vofit(), 0);
        factor *= m_beta;
        m_fluxStencil[dit()](vofit(), 0) *= factor;
      }
    }
  }

  /*****************/
  void
  DirichletConductivityEBBC::
  applyEBFlux(EBCellFAB&                    a_lphi,
              const EBCellFAB&              a_phi,
              VoFIterator&                  a_vofit,
              const LayoutData<IntVectSet>& a_cfivs,
              const DataIndex&              a_dit,
              const RealVect&               a_probLo,
              const RealVect&               a_dx,
              const Real&                   a_factor,
              const bool&                   a_useHomogeneous,
              const Real&                   a_time)
  {
    CH_TIME("DirichletConductivityEBBC::applyEBFlux");
    CH_assert(a_lphi.nComp() == 1 );
    CH_assert(a_phi.nComp()  == 1);

    const BaseIVFAB<Real>& poissWeight = (m_bc.getFluxWeight())[a_dit];
    Real value = 0.0;

    const EBISBox&   ebisBox = a_phi.getEBISBox();
    for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();
      if (m_bc.m_isFunction)
      {
        const RealVect& centroid = ebisBox.bndryCentroid(vof);
        const RealVect&   normal = ebisBox.normal(vof);

        value = m_bc.m_func->value(vof,centroid,normal,a_dx,a_probLo,a_dit,a_time,0);
      }
      else
      {
        if (m_bc.m_onlyHomogeneous)
        {
          MayDay::Error("DirichletConductivityEBBC::getFaceFlux called with undefined inhomogeneous BC");
        }

        value = m_bc.m_value;
      }
      Real poissWeightPt = poissWeight(vof, 0);
      const Real& areaFrac = ebisBox.bndryArea(vof);
      const Real& bcoef    = (*m_bcoe)[a_dit](vof,0);
      Real flux = poissWeightPt*value*areaFrac;
      Real compFactor = a_factor*bcoef*m_beta;
      a_lphi(vof,0) += flux * compFactor;
    }
  }
}
