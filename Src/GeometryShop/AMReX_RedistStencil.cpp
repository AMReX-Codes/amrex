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

#include "AMReX_RedistStencil.H"
#include "AMReX_EBArith.H"

namespace amrex
{
  RedistStencil::RedistStencil()
  {
    m_isDefined = false;
  }
    
  RedistStencil::~RedistStencil()
  {
  }
    
  RedistStencil::RedistStencil(const EBLevelGrid & a_eblg,
                               const int         & a_redistRadius)
  {
    define(m_eblg, a_redistRadius);
  }
    
  void RedistStencil::define(const EBLevelGrid & a_eblg,
                             const int         & a_redistRadius)
  {
    BL_PROFILE("RedistStencil::define");
    m_isDefined = true;
    m_hasDefaultWeights = true;
    m_eblg = a_eblg;
    m_redistRadius = a_redistRadius;
    m_stencil.define(m_eblg.getDBL(), m_eblg.getDM());
    m_volsten.define(m_eblg.getDBL(), m_eblg.getDM());
    for (MFIter  mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box region = m_eblg.getDBL()[mfi];
      region.grow(m_redistRadius);
      region &= m_eblg.getDomain();
      const EBISBox& ebisBox = m_eblg.getEBISL()[mfi];
      IntVectSet irregIVS = ebisBox.getIrregIVS(region);
      BaseIVFAB<VoFStencil >& stenFAB =    m_stencil[mfi];
      BaseIVFAB<VoFStencil >& volstenFAB = m_volsten[mfi];
      stenFAB.define(   irregIVS, ebisBox.getEBGraph(), 1);
      volstenFAB.define(irregIVS, ebisBox.getEBGraph(), 1);
      for (VoFIterator vofit(irregIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        VoFStencil thisSten;
        computePointStencil(thisSten, vof,  mfi);
        stenFAB(vof, 0) = thisSten;
        volstenFAB(vof, 0) = thisSten;
      }
    }
  }
    
  bool RedistStencil::isDefined() const
  {
    return m_isDefined;
  }
    
  int RedistStencil::getRedistRadius() const
  {
    return m_redistRadius;
  }
    
  void RedistStencil::resetWeights(const FabArray<EBCellFAB> & a_modifier,
                                   const int                 & a_ivar)
  {
    BL_ASSERT(isDefined());
    m_hasDefaultWeights = false;
    for (MFIter mfi(a_modifier); mfi.isValid(); ++mfi)
    {
      const EBISBox& ebisBox = m_eblg.getEBISL()[mfi];
      //initiate with the volume weighted stencil
      BaseIVFAB<VoFStencil >& stenFAB    =  m_stencil[mfi];
      BaseIVFAB<VoFStencil >& volstenFAB =  m_volsten[mfi];
      const EBCellFAB& modFAB            = a_modifier[mfi];
      const IntVectSet& irregIVS = volstenFAB.getIVS();
      for (VoFIterator vofit(irregIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        VoFStencil oldSten = volstenFAB(vof, 0);
        VoFStencil newSten;
        Real sum = 0.0;
        for (int isten = 0; isten < oldSten.size(); isten++)
        {
          const VolIndex& thatVoF = oldSten.vof(isten);
    
          Real weight  = modFAB(thatVoF, a_ivar);
          Real volfrac = ebisBox.volFrac(thatVoF);
          //it is weight*volfrac that is normalized
          sum += weight*volfrac;
          newSten.add(thatVoF, weight);
        }
        Real eps = 1.0e-12;
        if (std::abs(sum) > eps)
        {
          Real scaling = 1.0/sum;
          newSten *= scaling;
        }
        else
        {
          //if there is nowhere to put the mass
          newSten *= 0.;
        }
        stenFAB(vof, 0) = newSten;
      }
    }
  }
    
  const BaseIVFAB<VoFStencil>& RedistStencil::operator[](const MFIter& a_mfi) const
  {
    return m_stencil[a_mfi];
  }
    
  void RedistStencil::computePointStencil(VoFStencil     & a_stencil,
                                          const VolIndex & a_srcVoF,
                                          const MFIter   & a_mfi)
  {
    BL_PROFILE("RedistStencil::computePointStencil");
    
    const EBISBox& ebisBox = m_eblg.getEBISL()[a_mfi];
    
    //now set the weights according to the volumefrac/sum(volfrac)
    //you can reset the weights later if you like
    a_stencil.clear();
    Real sum = 0.0;
    
    //get the vofs.  these are the intvects
    //it must be called with the first time
    IntVect timesMoved = IntVect::TheZeroVector();
    IntVect pathSign   = IntVect::TheZeroVector();
    Array<VolIndex> vofsStencil;
    EBArith::getAllVoFsInMonotonePath(vofsStencil, timesMoved,
                                      pathSign, a_srcVoF, ebisBox,
                                      m_redistRadius);
    
    for (int isten = 0; isten < vofsStencil.size(); isten++)
    {
      const VolIndex& vof = vofsStencil[isten];
      Real weight = ebisBox.volFrac(vof);
      //w*volfrac is normalized
      //so the default weight is 1/sum
      sum += weight;
      //add with the sum of
      a_stencil.add(vof, 1.0);
    }
    
    //normalize the stencil
    if (std::abs(sum) > 1.0e-10)
    {
      Real scale = 1.0/sum;
      a_stencil *= scale;
    }
    //if there is not enough volume to which to redistribute
    //set the stencil to zero.   "Enough volume" shall be defined
    //as at least the volume of one full cell
    //
    if(std::abs(sum) < 1.0)
    {
      a_stencil *= 0.0;
    }
  }
}

