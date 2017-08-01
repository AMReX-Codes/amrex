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

#include "AMReX_EBConductivityOp.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBConductivityOpFactory.H"
#include "AMReX_EBCoarseAverage.H"

namespace amrex
{

  int EBConductivityOpFactory::s_testRef = 2;
  int EBConductivityOpFactory::s_maxBoxSize = 32;

  //-----------------------------------------------------------------------
  void
  nwoebcoCoarsenStuff(LevelData<EBCellFAB>               & a_acoefCoar,
                      LevelData<EBFluxFAB>               & a_bcoefCoar,
                      const EBLevelGrid                  & a_eblgFine,
                      const EBLevelGrid                  & a_eblgCoar,
                      const LevelData<EBCellFAB>         & a_acoefFine,
                      const LevelData<EBFluxFAB>         & a_bcoefFine,
                      const int                          & a_refToDepth)
  {
    BL_ASSERT(a_refToDepth > 0);
    BL_ASSERT(a_acoefFine.nGrow() == a_acoefCoar.nGrow());
    BL_ASSERT(a_acoefFine.nGrow() == a_bcoefCoar.nGrow());
    BL_ASSERT(a_acoefFine.nGrow() == a_bcoefFine.nGrow());

    if (a_refToDepth == 1)
    {
      a_acoefCoar.copy(a_acoefFine 0, 0, 1, 0, 0);
      a_bcoefCoar.copy(a_bcoefCoar 0, 0, 1, 0, 0);
    }
    else
    {
      bool useKappaWeightingInStencil = true;
      bool enableFaceAveraging = true;
      EBCoarseAverage averageOp(a_eblgFine, a_eblgCoar, a_refToDepth, 
                                a_acoefFine.nGrow(),
                                useKappaWeightingInStencil, enableFaceAveraging);


      //amrex::Warning("might want to figure out what harmonic averaging is in this context");
      averageOp.average( a_acoefCoar     ,  a_acoefFine     , 0, 0, 1);
      averageOp.average( a_bcoefCoar     ,  a_bcoefFine     , 0, 0, 1);
      averageOp.average( a_bcoefCoarIrreg,  a_bcoefFineIrreg, 0, 0, 1);
    }
  }
  //-----------------------------------------------------------------------
  EBConductivityOpFactory::~EBConductivityOpFactory()
  {
  }
  //-----------------------------------------------------------------------
  int
  EBConductivityOpFactory::
  refToFiner(const ProblemDomain& a_domain) const
  {
    int retval = -1;
    bool found = false;
    for (int ilev = 0; ilev < m_eblgs.size(); ilev++)
    {
      if (m_eblgs[ilev].getDomain() == a_domain)
      {
        retval = m_refRatio[ilev];
        found = true;
      }
    }
    if (!found)
    {
      amrex::Error("Domain not found in AMR hierarchy");
    }
    return retval;
  }
  //-----------------------------------------------------------------------
  EBConductivityOpFactory::
  EBConductivityOpFactory(onst vector<EBLevelGrid>&                                   a_eblgs,
                          const Real&                                                 a_alpha,
                          const Real&                                                 a_beta,
                          const vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                          const vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                          const Real&                                                 a_dxCoarse,
                          const vector<int>&                                          a_refRatio,
                          const RefCountedPtr<ConductivityBaseDomainBCFactory>&       a_domainBCFactory,
                          const RefCountedPtr<ConductivityBaseEBBCFactory>    &       a_ebBCFactory,
                          const int     &                                             a_ghostCellsPhi,
                          const int     &                                             a_ghostCellsRhs,
                          int a_numLevels)
  {
    BL_ASSERT(a_eblgs.size() <= a_refRatio.size());
    m_dataBased = false;
    m_ghostCellsRhs = a_ghostCellsRhs;
    m_ghostCellsPhi = a_ghostCellsPhi;
    m_acoef         = a_acoef;
    m_bcoef         = a_bcoef;
    m_alpha         = a_alpha;
    m_beta          = a_beta;
    if (a_numLevels > 0)
    {
      m_numLevels = a_numLevels;
    }
    else
    {
      m_numLevels = a_eblgs.size();
    }

    m_domainBCFactory = a_domainBCFactory;
    m_ebBCFactory     = a_ebBCFactory;

    m_eblgs           = a_eblgs;
    m_refRatio        = a_refRatio;
    m_dx.resize(m_numLevels);

    m_dx[0] = a_dxCoarse;
    for (int ilev = 1; ilev < m_numLevels; ilev++)
    {
      m_dx[ilev] = m_dx[ilev-1];
      m_dx[ilev] /= m_refRatio[ilev-1];
    }
    m_alpha = a_alpha;
    m_beta = a_beta;
  }
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  EBConductivityOp*
  EBConductivityOpFactory::
  MGnewOp(const ProblemDomain& a_domainFine,
          int                  a_depth,
          bool                 a_homoOnly)
  {
    ParmParse pp;
    bool turn_off_mg = false;
    pp.query("turn_off_multigrid_for_nwoc", turn_off_mg);
    if(turn_off_mg)
    {
      pout() << "turn off EBMultiGrid for EBConductivityOp because turn_off_multigrid_for_nwoc = true " << endl;
      return NULL;
    }

    //find out if there is a real starting point here.
    int ref=-1;
    bool found = false;

    RefCountedPtr<LevelData<EBCellFAB> >               acoef;
    RefCountedPtr<LevelData<EBFluxFAB> >               bcoef;
    RefCountedPtr<LevelData<BaseIVFAB<Real> > >   bcoefIrreg;

    for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_eblgs[ilev].getDomain())
      {
        found = true;
        ref = ilev ;
      }
    }
    if (!found)
    {
      amrex::Error("No corresponding AMRLevel to starting point of MGnewOp");
    }
    

    EBConductivityOp* newOp = createOperator(ref, a_depth, false);
    
    return newOp;

  }
  //-----------------------------------------------------------------------
  EBConductivityOp*
  EBConductivityOpFactory::
  createOperator(const int&  a_amrlevel,
                 const int&  a_depth,
                 const bool  a_amrLevelOp)
  {
    EBLevelGrid                           eblg;
    EBLevelGrid                           eblgFine;
    EBLevelGrid                           eblgCoar;
    int                                   refToFine=-1;
    int                                   refToCoar=-1;
    bool hasFiner   = false;
    bool hasCoarser = false;
    int refRatMG = 1;
    for(int ipow = 0; ipow < a_depth; ipow++)
    {
      refRatMG *= 2;
    }

    if(a_amrLevelOp)
    {
      eblg = m_eblgs[a_amrlevel];
      if(a_amrLevel > 0)
      {
        eblgCoar = m_eblgs[a_amrlevel-1];
        refToCoar = m_refRatio[a_amrlevel-1]
        hasCoarser = true;
      }
      if(a_amrLevel < (m_numLevels - 1))
      {
        eblgFine = m_eblgs[a_amrlevel+1];
        refToFine = m_refRatio[a_amrlevel]
        hasFiner = true;
      }
    }
    else
    {
      //check to see if we are able to create an EBLG at this depth.   If not, return NULL
      if(a_depth == 0)
      {
        eblg = m_eblgs[a_amrlevel];
      }
      else
      {
        bool coarsenable = EBArith::createCoarserEBLG(eblg,  m_eblgs[a_amrlevel], refRatMG,  s_testRef, s_maxBoxSize);
        if(!coarsenable)
        {
          return NULL;
        }
      }
    }

    //now see if there is a coarser MG level
    bool hasCoarMG = EBArith::createCoarserEBLG(eblgCoarMG , eblg, 2, s_testRef, s_maxBoxSize);
    Real          dxMGLevel = m_dx[a_amrlevel]*refRatMG;
    Real dxCoar = -1.0;
    if (a_amrlevel > 0)
    {
      dxCoar = m_dx[a_amrlevel-1];
    }

    shared_ptr<ConductivityBaseDomainBC>  domainBC(m_domainBCFactory->create());
    shared_ptr<ConductivityBaseEBBC>      ebBC(m_ebbcFactory->create());
    shared_ptr<FabArray<EBCellFAB> >      acoef;
    shared_ptr<FabArray<EBFluxFAB> >      bcoef;

    //get coefficients
    nwoebcoCoarsenStuff(acoef,
                        bcoef,
                        m_eblgs[a_amrlevel],
                        eblg,
                        m_acoef[a_amrlevel],
                        m_bcoef[a_amrlevel],
                        refRatMG);


    EBConductivityOp* newOp = 
      new EBConductivityOp(eblgFine,
                           eblg,
                           eblgCoar,
                           eblgCoarMG,
                           domainBC,
                           ebBC,
                           dxMGLevel,
                           dxCoar,
                           refToFine,
                           refToCoar,
                           hasFiner,
                           hasCoarser,
                           hasCoarMG,
                           m_alpha,
                           m_beta,
                           acoef,
                           bcoef,
                           m_ghostCellsPhi,
                           m_ghostCellsRHS);
      }    

  //-----------------------------------------------------------------------
  EBConductivityOp*
  EBConductivityOpFactory::
  AMRnewOp(const ProblemDomain& a_domainFine)
  {
    //figure out which level we are at.
    int ref=-1;
    bool found = false;
    EBLevelGrid eblgFine, eblgCoar;
    for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_eblgs[ilev].getDomain())
      {
        found = true;
        ref = ilev ;
      }
    }
    if (!found)
    {
      amrex::Error("No corresponding AMRLevel to starting point of AMRnewOp");
    }

    EBConductivityOp* newOp = createOperator(ref, 0, true);

    return newOp;
  }
}
