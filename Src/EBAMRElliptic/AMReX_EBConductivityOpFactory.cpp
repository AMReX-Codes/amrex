#include "AMReX_EBConductivityOp.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBConductivityOpFactory.H"
#include "AMReX_EBCoarseAverage.H"
#include "AMReX_ParmParse.H"

namespace amrex
{

  int EBConductivityOpFactory::s_testRef = 2;
  int EBConductivityOpFactory::s_maxBoxSize = 32;

  //-----------------------------------------------------------------------
  void
  EBConductivityOpFactory::
  nwoebcoCoarsenStuff(FabArray<EBCellFAB>               & a_acoefCoar,
                      FabArray<EBFluxFAB>               & a_bcoefCoar,
                      const EBLevelGrid                 & a_eblgFine,
                      const EBLevelGrid                 & a_eblgCoar,
                      const FabArray<EBCellFAB>         & a_acoefFine,
                      const FabArray<EBFluxFAB>         & a_bcoefFine,
                      const int                         & a_refToDepth)
  {
    BL_ASSERT(a_refToDepth > 0);
    BL_ASSERT(a_acoefFine.nGrow() == m_ghost);
    BL_ASSERT(a_bcoefFine.nGrow() == m_ghost);

    if (a_refToDepth == 1)
    {
      a_acoefCoar.copy(a_acoefFine, 0, 0, 1, 0, 0);
      a_bcoefCoar.copy(a_bcoefFine, 0, 0, 1, 0, 0);
    }
    else
    {
      bool useKappaWeightingInStencil = true;
      bool enableFaceAveraging = true;
      EBCoarseAverage averageOp(a_eblgFine, a_eblgCoar, a_refToDepth, 
                                m_ghost,
                                useKappaWeightingInStencil, enableFaceAveraging);


      //amrex::Warning("might want to figure out what harmonic averaging is in this context");
      averageOp.average( a_acoefCoar     ,  a_acoefFine     , 0, 0, 1);
      averageOp.average( a_bcoefCoar     ,  a_bcoefFine     , 0, 0, 1);
    }
    a_acoefCoar.FillBoundary();
    a_bcoefCoar.FillBoundary();
  }
  //-----------------------------------------------------------------------
  EBConductivityOpFactory::~EBConductivityOpFactory()
  {
  }
  //-----------------------------------------------------------------------
  EBConductivityOpFactory::
  EBConductivityOpFactory(const vector<EBLevelGrid>&                              a_eblgs,
                          const Real&                                             a_alpha,
                          const Real&                                             a_beta,
                          const vector<shared_ptr<FabArray<EBCellFAB> > >&        a_acoef,
                          const vector<shared_ptr<FabArray<EBFluxFAB> > >&        a_bcoef,
                          const Real&                                             a_dxCoarse,
                          const vector<int>&                                      a_refRatio,
                          const shared_ptr<ConductivityBaseDomainBCFactory>&      a_domainBCFactory,
                          const shared_ptr<ConductivityBaseEBBCFactory>    &      a_ebBCFactory,
                          const int     &                                         a_ghost,
                          int a_numLevels)
  {
    BL_ASSERT(a_eblgs.size() <= a_refRatio.size());

    m_ghost         = a_ghost;
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
  MGnewOp(const Box& a_domainFine,
          int        a_depth,
          bool       a_homoOnly)
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

    shared_ptr<FabArray<EBCellFAB> >               acoef;
    shared_ptr<FabArray<EBFluxFAB> >               bcoef;
    shared_ptr<FabArray<BaseIVFAB<Real> > >   bcoefIrreg;

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
  createOperator(const int&  a_amrLevel,
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
      eblg = m_eblgs[a_amrLevel];
      if(a_amrLevel > 0)
      {
        eblgCoar = m_eblgs[a_amrLevel-1];
        refToCoar = m_refRatio[a_amrLevel-1];
        hasCoarser = true;
      }
      if(a_amrLevel < (m_numLevels - 1))
      {
        eblgFine = m_eblgs[a_amrLevel+1];
        refToFine = m_refRatio[a_amrLevel];
        hasFiner = true;
      }
    }
    else
    {
      //check to see if we are able to create an EBLG at this depth.   If not, return NULL
      if(a_depth == 0)
      {
        eblg = m_eblgs[a_amrLevel];
      }
      else
      {
        bool coarsenable = EBArith::createCoarserEBLG(eblg,  m_eblgs[a_amrLevel], refRatMG,  s_testRef, s_maxBoxSize);
        if(!coarsenable)
        {
          return NULL;
        }
      }
    }

    //now see if there is a coarser MG level
    EBLevelGrid eblgCoarMG;
    bool hasCoarMG = EBArith::createCoarserEBLG(eblgCoarMG , eblg, 2, s_testRef, s_maxBoxSize);
    Real          dxMGLevel = m_dx[a_amrLevel]*refRatMG;
    Real dxCoar = -1.0;
    if (a_amrLevel > 0)
    {
      dxCoar = m_dx[a_amrLevel-1];
    }

    shared_ptr<ConductivityBaseDomainBC>  domainBC(m_domainBCFactory->create());
    shared_ptr<ConductivityBaseEBBC>      ebBC(        m_ebBCFactory->create());

    //these have to be the same
    EBCellFactory cellfact(eblg.getEBISL());
    EBFluxFactory fluxfact(eblg.getEBISL());
    shared_ptr<FabArray<EBCellFAB> > acoef(new FabArray<EBCellFAB>(eblg.getDBL(), eblg.getDM(), 1, m_ghost, MFInfo(), cellfact));
    shared_ptr<FabArray<EBFluxFAB> > bcoef(new FabArray<EBFluxFAB>(eblg.getDBL(), eblg.getDM(), 1, m_ghost, MFInfo(), fluxfact));
    
    //get coefficients
    nwoebcoCoarsenStuff(*acoef,
                        *bcoef,
                        m_eblgs[a_amrLevel],
                        eblg,
                        *m_acoef[a_amrLevel],
                        *m_bcoef[a_amrLevel],
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
                           m_ghost);

    return newOp;
  }    
  
  //-----------------------------------------------------------------------
  EBConductivityOp*
  EBConductivityOpFactory::
  AMRnewOp(const Box& a_domainFine)
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
