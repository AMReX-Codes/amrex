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

#include "AMReX_VCAggStencil.H"

namespace amrex
{
/**************/
  VCAggStencil::
  VCAggStencil(const vector<shared_ptr<BaseIndex>   >    & a_dstVoFs,
               const vector<shared_ptr<BaseStencil> >    & a_vofStencil,
               const EBCellFAB                           & a_phiData,
               const EBCellFAB                           & a_rhsData,
               const EBCellFAB                           & a_relCoef,
               const BaseIVFAB<Real>                     & a_alphaWt,
               const int                                 & a_ncomp)
    :AggStencil<EBCellFAB, EBCellFAB>(a_dstVoFs, a_vofStencil, a_phiData, a_rhsData)
  {
    BL_PROFILE("VCAggSten.constructor");

    m_phiAccess.resize(a_dstVoFs.size());
    m_relAccess.resize(a_dstVoFs.size());
    m_alpAccess.resize(a_dstVoFs.size());
    m_iv.resize(a_dstVoFs.size());
    for (int idst = 0; idst < a_dstVoFs.size(); idst++)
    {
      const BaseIndex& dstVoF = *a_dstVoFs[idst];
      m_phiAccess[idst].dataID = a_phiData.dataType(dstVoF);
      m_phiAccess[idst].offset = a_phiData.offset(dstVoF, 0);
      m_relAccess[idst].dataID = a_relCoef.dataType(dstVoF);
      m_relAccess[idst].offset = a_relCoef.offset(dstVoF, 0);
      m_alpAccess[idst].dataID = a_alphaWt.dataType(dstVoF);
      m_alpAccess[idst].offset = a_alphaWt.offset(dstVoF, 0);

      VolIndex* vofPtr = dynamic_cast<VolIndex*>(&(*a_dstVoFs[idst]));
      if(vofPtr == NULL)
      {
        amrex::Error("dynamic cast error--VCAggStencil just handles VolIndicies");
      }
      m_iv[idst] = vofPtr->gridIndex();
    }
    m_cachePhi.resize( m_ebstencil.size(), vector<Real>(a_ncomp, 0.));
  }
/************/
  void
  VCAggStencil::
  cachePhi(const EBCellFAB& a_phi) const
  {
    BL_PROFILE("VCAggStencil::cachePhi");
    vector<const Real*> dataPtrsPhi(a_phi.numDataTypes());
    for (int ivar = 0; ivar < a_phi.nComp(); ivar++)
    {
      for (int ivec = 0; ivec < dataPtrsPhi.size(); ivec++)
      {
        dataPtrsPhi[ivec] = a_phi.dataPtr(ivec, ivar);
      }

      for (int idst = 0; idst < m_ebstencil.size(); idst++)
      {
        const Real* phiPtr =  dataPtrsPhi[m_phiAccess[idst].dataID] + m_phiAccess[idst].offset;
        m_cachePhi[idst][ivar] = *phiPtr;
      }
    }
  }
/**************/
  void
  VCAggStencil::
  uncachePhi(EBCellFAB& a_phi) const
  {
    BL_PROFILE("VCAggSten::uncache");
    vector<Real*> dataPtrsPhi(a_phi.numDataTypes());
    for (int ivar = 0; ivar < a_phi.nComp(); ivar++)
    {
      for (int ivec = 0; ivec < dataPtrsPhi.size(); ivec++)
      {
        dataPtrsPhi[ivec] = a_phi.dataPtr(ivec, ivar);
      }

      for (int idst = 0; idst < m_ebstencil.size(); idst++)
      {
        Real* phiPtr =  dataPtrsPhi[m_phiAccess[idst].dataID] + m_phiAccess[idst].offset;
        *phiPtr = m_cachePhi[idst][ivar];
      }
    }
  }
/**************/
  void 
  VCAggStencil::
  relax(EBCellFAB              & a_phi,
        const EBCellFAB        & a_rhs,
        const EBCellFAB        & a_relCoef,
        const BaseIVFAB<Real>  & a_alphaWt,
        const Real             & a_alpha,
        const Real             & a_beta,
        const int              & a_varDest,
        const IntVect          & a_color)
  {
    BL_PROFILE("VCAggSten::relax");
    const int numtyperhs = a_rhs.numDataTypes();
    const int numtypephi = a_phi.numDataTypes();
    const int numtyperel = a_relCoef.numDataTypes();
    const int numtypealp = a_alphaWt.numDataTypes();
  
    vector<const Real*> dataPtrsRhs(numtyperhs);
    vector<const Real*> dataPtrsRel(numtyperel);
    vector<const Real*> dataPtrsAlp(numtypealp);
    //phi is the source (what the stencil gets applied to)
    //and the destination (where the answer goes).   For the 
    //source, we need  the variable to be zero because the stencil
    //variable is taken into account in aggstencil
    vector<const Real*> dataPtrsSrc(numtypephi);
    vector<Real*>       dataPtrsDst(numtypephi);
    int varDst = a_varDest;
    int varSrc = 0;  //stencil variable taken into account in aggstencil
    for (int ivec = 0; ivec < numtyperhs; ivec++)
    {
      //rhs has the same variable number as destination (in this case phi)
      dataPtrsRhs[ivec] = a_rhs.dataPtr(ivec, varDst);
    }
    for (int ivec = 0; ivec < numtypealp; ivec++)
    {
      //alphaweight has the same variable number as destination (in this case phi)
      dataPtrsAlp[ivec] = a_alphaWt.dataPtr(ivec, varDst);
    }
    for (int ivec = 0; ivec < numtypephi; ivec++)
    {
      dataPtrsSrc[ivec] = a_phi.dataPtr(ivec, varSrc);
    }
    for (int ivec = 0; ivec < numtypephi; ivec++)
    {
      dataPtrsDst[ivec] = a_phi.dataPtr(ivec, varDst);
    }
    for (int ivec = 0; ivec < numtypephi; ivec++)
    {
      //relaxation coeff has the same variable number as as destination (in this case phi)
      dataPtrsRel[ivec] = a_relCoef.dataPtr(ivec, varDst);
    }

    for (int idst = 0; idst < m_ebstencil.size(); idst++)
    {
      const IntVect& iv = m_iv[idst];
      bool doThisVoF = true;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (iv[idir] % 2 != a_color[idir])
        {
          doThisVoF = false;
          break;
        }
      }
      if(doThisVoF)
      {
//begin debug
        int idebug = 0;
        if((iv[0]==14) && (iv[1]==12))
        {
          idebug = 1;
        }
//end debug

        const Real* rhsiPtr =  dataPtrsRhs[m_dstAccess[idst].dataID] + m_dstAccess[idst].offset;
        const Real& rhsi = *rhsiPtr;
        const Real* relcoPtr =  dataPtrsRel[m_relAccess[idst].dataID] + m_relAccess[idst].offset;
        const Real& relco = *relcoPtr;
        const Real* alpWtPtr =  dataPtrsAlp[m_alpAccess[idst].dataID] + m_alpAccess[idst].offset;
        const Real& alphaWeight = *alpWtPtr;


        const stencil_t& ebstencil = m_ebstencil[idst];
        Real lphi = 0;
        for (int isten = 0; isten < ebstencil.size(); isten++)
        {
          const Real& weight = ebstencil[isten].second;
          const long& offset = ebstencil[isten].first.offset;
          const int & dataID = ebstencil[isten].first.dataID;
          const Real& phiVal = *(dataPtrsSrc[dataID] + offset);
          lphi += phiVal*weight;
        }
        Real* phiiPtr =  dataPtrsDst[m_phiAccess[idst].dataID] + m_phiAccess[idst].offset;
        Real& phii = *phiiPtr;
        //multiply by beta and add in identity term
        lphi = a_beta*lphi + a_alpha*alphaWeight*phii;

        phii = phii + relco*(rhsi - lphi);

//begin debug
        idebug = 0;
//end debug
      }
    }
  }



/**************/
  void 
  VCAggStencil::
  apply(EBCellFAB              & a_lph,
        const EBCellFAB        & a_phi,
        const BaseIVFAB<Real>  & a_alp, //alphaDiagWeight
        const Real             & a_alpha,
        const Real             & a_beta,
        const int              & a_varDest,
        const bool             & a_incrementOnly)
  {
    BL_PROFILE("VCAggSten::apply");
    //this makes lphi = divf
    AggStencil<EBCellFAB,EBCellFAB>::apply(a_lph, a_phi, a_varDest, a_incrementOnly);

    const int numtypelph = a_lph.numDataTypes();
    const int numtypephi = a_phi.numDataTypes();
    const int numtypealp = a_alp.numDataTypes();
  
    vector<      Real*> dataPtrsLph(numtypelph);
    vector<const Real*> dataPtrsPhi(numtypephi);
    vector<const Real*> dataPtrsAlp(numtypealp);
    //no stencils here so everything is on the same var
    for (int ivec = 0; ivec < numtypelph; ivec++)
    {
      dataPtrsLph[ivec] = a_lph.dataPtr(ivec, a_varDest);
    }
    for (int ivec = 0; ivec < numtypealp; ivec++)
    {
      dataPtrsAlp[ivec] = a_alp.dataPtr(ivec, a_varDest);
    }
    for (int ivec = 0; ivec < numtypephi; ivec++)
    {
      dataPtrsPhi[ivec] = a_phi.dataPtr(ivec, a_varDest);
    }

    for (int idst = 0; idst < m_ebstencil.size(); idst++)
    {
      Real*        lphiPtr    =  dataPtrsLph[m_dstAccess[idst].dataID] + m_dstAccess[idst].offset;
      Real&           lphi    = *lphiPtr;
      const Real* alpWtPtr    =  dataPtrsAlp[m_alpAccess[idst].dataID] + m_alpAccess[idst].offset;
      const Real& alphaWeight = *alpWtPtr;
      const Real* phiiPtr     =  dataPtrsPhi[m_phiAccess[idst].dataID] + m_phiAccess[idst].offset;
      const Real& phii        = *phiiPtr;
      //multiply by beta and add in identity term
      lphi = a_beta*lphi + a_alpha*alphaWeight*phii;
    }
  }
}

