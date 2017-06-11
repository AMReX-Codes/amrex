#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoadBalance.H"
#include "EBArith.H"

#include "NWOEBConductivityOp.H"
#include "EBArith.H"
#include "ParmParse.H"
#include "CH_Timer.H"
#include "NWOEBConductivityOpFactory.H"
#include "EBCoarseAverage.H"
#include "NamespaceHeader.H"

int NWOEBConductivityOpFactory::s_testRef = 2;
int NWOEBConductivityOpFactory::s_maxBoxSize = 32;

//-----------------------------------------------------------------------
void
nwoebcoCoarsenStuff(LevelData<EBCellFAB>               & a_acoefCoar,
                    LevelData<EBFluxFAB>               & a_bcoefCoar,
                    LevelData<BaseIVFAB<Real> >        & a_bcoefCoarIrreg,
                    const EBLevelGrid                  & a_eblgFine,
                    const EBLevelGrid                  & a_eblgCoar,
                    const LevelData<EBCellFAB>         & a_acoefFine,
                    const LevelData<EBFluxFAB>         & a_bcoefFine,
                    const LevelData<BaseIVFAB<Real> >  & a_bcoefFineIrreg,
                    const int                          & a_refToDepth)
{
  CH_assert(a_refToDepth > 0);

  Interval interv(0, 0);
  if (a_refToDepth == 1)
    {
      a_acoefFine.     copyTo(interv,  a_acoefCoar,     interv);
      a_bcoefFine.     copyTo(interv,  a_bcoefCoar,     interv);
      a_bcoefFineIrreg.copyTo(interv,  a_bcoefCoarIrreg,interv);
    }
  else
    {
      EBCoarseAverage averageOp(a_eblgFine.getDBL(),    a_eblgCoar.getDBL(),
                                a_eblgFine.getEBISL(),  a_eblgCoar.getEBISL(),
                                a_eblgCoar.getDomain(), a_refToDepth, 1,
                                a_eblgCoar.getEBIS());

      //MayDay::Warning("might want to figure out what harmonic averaging is in this context");
      averageOp.average( a_acoefCoar     ,  a_acoefFine     , interv);
      averageOp.average( a_bcoefCoar     ,  a_bcoefFine     , interv);
      averageOp.average( a_bcoefCoarIrreg,  a_bcoefFineIrreg, interv);
    }
  a_acoefCoar.exchange(Interval(0,0));
  a_bcoefCoar.exchange(Interval(0,0));
  a_bcoefCoarIrreg.exchange(Interval(0,0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
NWOEBConductivityOpFactory::~NWOEBConductivityOpFactory()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
NWOEBConductivityOpFactory::
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
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
NWOEBConductivityOpFactory::
NWOEBConductivityOpFactory(const Vector<EBLevelGrid>&                               a_eblgs,
                        const Vector<RefCountedPtr<NWOEBQuadCFInterp> >&            a_quadCFI,
                        const Real&                                                 a_alpha,
                        const Real&                                                 a_beta,
                        const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                        const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                        const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg,
                        const Real&                                                 a_dxCoarse,
                        const Vector<int>&                                          a_refRatio,
                        const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                        const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBCFactory,
                        const IntVect&                                              a_ghostCellsPhi,
                        const IntVect&                                              a_ghostCellsRhs,
                        const int &                                                 a_relaxType,
                        int a_numLevels)
{
  CH_assert(a_eblgs.size() <= a_refRatio.size());
  m_dataBased = false;
  m_relaxType = a_relaxType;
  m_quadCFI = a_quadCFI;
  m_ghostCellsRhs = a_ghostCellsRhs;
  m_ghostCellsPhi = a_ghostCellsPhi;
  m_acoef         = a_acoef;
  m_bcoef         = a_bcoef;
  m_bcoefIrreg    = a_bcoefIrreg;
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

  m_eblgsMG.resize(m_numLevels);
  m_acoefMG.resize(m_numLevels);
  m_bcoefMG.resize(m_numLevels);
  m_bcoefIrregMG.resize(m_numLevels);
  m_hasMGObjects.resize(m_numLevels);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
  {
    if ((ilev==0) || (m_refRatio[ilev] > 2))
    {
      m_hasMGObjects[ilev] = (s_testRef<s_maxBoxSize);

      int mgRef = 2;
      m_eblgsMG[ilev]     .resize(0);
      m_acoefMG[ilev]     .resize(0);
      m_bcoefMG[ilev]     .resize(0);
      m_bcoefIrregMG[ilev].resize(0);
      m_eblgsMG[ilev]     .push_back(m_eblgs[ilev]);
      m_acoefMG[ilev]     .push_back(m_acoef[ilev]);
      m_bcoefMG[ilev]     .push_back(m_bcoef[ilev]);
      m_bcoefIrregMG[ilev].push_back(m_bcoefIrreg[ilev]);

      bool hasCoarser = true;
      hasCoarser = (s_testRef<s_maxBoxSize);
      while (hasCoarser)
      {
        int imgsize = m_eblgsMG[ilev].size();
        const EBLevelGrid& eblgFine=  m_eblgsMG[ilev][imgsize-1];
        DisjointBoxLayout dblCoarMG;
        ProblemDomain  domainCoarMG;
        bool dumbool;
        hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarMG,
                                                       domainCoarMG,
                                                       eblgFine.getDBL(),
                                                       eblgFine.getEBISL(),
                                                       eblgFine.getDomain(),
                                                       mgRef,
                                                       eblgFine.getEBIS(),
                                                       s_maxBoxSize,
                                                       dumbool,
                                                       s_testRef);

        if (hasCoarser)
        {
          m_eblgsMG[ilev].push_back(EBLevelGrid(dblCoarMG, domainCoarMG, 4, eblgFine.getEBIS()));
          int img = m_eblgsMG[ilev].size() - 1;
          const EBLevelGrid& eblgCoar = m_eblgsMG[ilev][img  ];
          const EBLevelGrid& eblgFine = m_eblgsMG[ilev][img-1];

          int nghost = 1;
          LayoutData<IntVectSet> irregSets(eblgCoar.getDBL());
          for (DataIterator dit = eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            Box grownBox = grow(eblgCoar.getDBL().get(dit()), nghost);
            grownBox &= domainCoarMG;
            irregSets[dit()] = eblgCoar.getEBISL()[dit()].getIrregIVS(grownBox);
          }
          EBFluxFactory       ebfluxfact(eblgCoar.getEBISL());
          EBCellFactory       ebcellfact(eblgCoar.getEBISL());
          BaseIVFactory<Real> baseivfact(eblgCoar.getEBISL(), irregSets);

          RefCountedPtr<LevelData<BaseIVFAB<Real> > > bcoefIrregCoar( new LevelData<BaseIVFAB<Real> >(eblgCoar.getDBL(), 1, nghost*IntVect::Unit, baseivfact) );
          RefCountedPtr<LevelData<EBCellFAB> >             acoefCoar( new LevelData<EBCellFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebcellfact) );
          RefCountedPtr<LevelData<EBFluxFAB> >             bcoefCoar( new LevelData<EBFluxFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebfluxfact) );

          nwoebcoCoarsenStuff(*acoefCoar, *bcoefCoar, *bcoefIrregCoar, eblgFine, eblgCoar,
                              *m_acoefMG[ilev][img-1], *m_bcoefMG[ilev][img-1], *m_bcoefIrregMG[ilev][img-1], mgRef);

          m_acoefMG[ilev].push_back(acoefCoar);
          m_bcoefMG[ilev].push_back(bcoefCoar);
          m_bcoefIrregMG[ilev].push_back(bcoefIrregCoar);
        }
      }

    }
    else
    {
      m_hasMGObjects[ilev] = false;
    }
  }
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
NWOEBConductivityOp*
NWOEBConductivityOpFactory::
MGnewOp(const ProblemDomain& a_domainFine,
        int                  a_depth,
        bool                 a_homoOnly)
{
  ParmParse pp;
  bool turn_off_mg = false;
  pp.query("turn_off_multigrid_for_nwoc", turn_off_mg);
  if(turn_off_mg)
    {
      pout() << "turn off multigrid for NWOEBConductivityOp because turn_off_multigrid_for_nwoc = true " << endl;
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
    MayDay::Error("No corresponding AMRLevel to starting point of MGnewOp");
  }

  //multigrid operator.  coarn by two from depth  no
  //coarse or finer to worry about.
  Real          dxMGLevel;
  EBLevelGrid eblgMGLevel;
  EBLevelGrid eblgCoarMG;
  RefCountedPtr<NWOEBQuadCFInterp> quadCFI; //only defined if on an amr level
  Real dxCoar = 1.0;
  dxCoar *= -1.0;
  //  int refToDepth = 1;
  if (ref > 0)
  {
    dxCoar = m_dx[ref-1];
  }
  bool hasCoarMGObjects = false;
  if (a_depth == 0)
  {
    eblgMGLevel    = m_eblgs[ref];

    acoef = m_acoef[ref];

    bcoef          = m_bcoef[ref];
    bcoefIrreg     = m_bcoefIrreg[ref];
    dxMGLevel      = m_dx[ref];
    quadCFI        = m_quadCFI[ref];

    hasCoarMGObjects = m_hasMGObjects[ref];
    if (hasCoarMGObjects)
    {
      eblgCoarMG = m_eblgsMG[ref][1];
    }
  }
  else
  {
    int icoar = 1;
    for (int idep = 0; idep < a_depth; idep++)
    {
      icoar *= 2;
    }
    //    refToDepth = icoar;
    const ProblemDomain domainFine = m_eblgs[ref].getDomain();
    ProblemDomain domainBoxMGLevel = coarsen(domainFine, icoar);
    bool foundMGLevel = false;
    int numMGLevels = m_eblgsMG[ref].size();
    for (int img = 0; img < numMGLevels; img++)
    {
      if (m_eblgsMG[ref][img].getDomain() == domainBoxMGLevel)
      {
        eblgMGLevel = m_eblgsMG[ref][img];
        acoef = m_acoefMG[ref][img];
        bcoef = m_bcoefMG[ref][img];
        bcoefIrreg  = m_bcoefIrregMG[ref][img];
        foundMGLevel = true;

        hasCoarMGObjects = ((img+1) < (numMGLevels));
        if (hasCoarMGObjects)
        {
          eblgCoarMG = m_eblgsMG[ref][img+1];
        }
        break;
      }
    }
    bool coarsenable = foundMGLevel;

    dxMGLevel = m_dx[ref];
    dxMGLevel *= Real(icoar);

    if (!coarsenable)
    {
      //not coarsenable.
      //return null
      return NULL;
    }
  }

  ConductivityBaseEBBC*     viscEBBC  = (ConductivityBaseEBBC*)         m_ebBCFactory->create(eblgMGLevel.getDomain(), eblgMGLevel.getEBISL(),
      dxMGLevel*RealVect::Unit,
      &m_ghostCellsPhi, &m_ghostCellsRhs);
  ConductivityBaseDomainBC* viscDomBC = (ConductivityBaseDomainBC*) m_domainBCFactory->create(eblgMGLevel.getDomain(), eblgMGLevel.getEBISL(),
      dxMGLevel*RealVect::Unit);
  RefCountedPtr<ConductivityBaseEBBC>      ebbc(viscEBBC);
  RefCountedPtr<ConductivityBaseDomainBC> dombc(viscDomBC);
  //no need for finer or coarser amr levels here
  bool hasFine   = false;
  bool hasCoar   = false;
  int bogRef = 2;
  //hook for optimization.
  //need to put this in  ala EBAMRPoissonOp to get it.
  bool layoutChanged = true;
  //
  NWOEBConductivityOp* newOp = NULL;
  // Time-independent A coefficient.
  newOp = new NWOEBConductivityOp(EBLevelGrid(), eblgMGLevel, EBLevelGrid(), eblgCoarMG, quadCFI,
                               dombc, ebbc, dxMGLevel, dxCoar, bogRef, bogRef, hasFine, hasCoar,
                               hasCoarMGObjects, layoutChanged, m_alpha, m_beta,
                               acoef, bcoef, bcoefIrreg, m_ghostCellsPhi, m_ghostCellsRhs, m_relaxType);

  return newOp;

}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
NWOEBConductivityOp*
NWOEBConductivityOpFactory::
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
      MayDay::Error("No corresponding AMRLevel to starting point of AMRnewOp");
    }

  int refToFiner   = 2;
  int refToCoarser = 2;
  Real dxCoar = -1.0;
  if (ref > 0)
    {
      eblgCoar = m_eblgs[ref-1];
      dxCoar = m_dx[ref-1];
      refToCoarser= m_refRatio[ref-1];
    }
  if (ref < m_numLevels-1)
    {
      eblgFine = m_eblgs[ref+1];
      refToFiner = m_refRatio[ref];
    }
  //creates coarse and finer info and bcs and all that
  EBLevelGrid      eblgMGLevel = m_eblgs[ref];
  Real               dxMGLevel =    m_dx[ref];

  bool hasCoarMGObjects = m_hasMGObjects[ref];
  EBLevelGrid eblgCoarMG;
  if (hasCoarMGObjects)
    {
      eblgCoarMG = m_eblgsMG[ref][1];
    }
  ConductivityBaseEBBC*     viscEBBC  = (ConductivityBaseEBBC*)     m_ebBCFactory->create(    m_eblgs[ref].getDomain(),
                                                                                              m_eblgs[ref].getEBISL(), dxMGLevel*RealVect::Unit,
                                                                                              &m_ghostCellsPhi, &m_ghostCellsRhs);
  if (m_dataBased)
    {
      viscEBBC->setData(m_data[ref]);
    }

  ConductivityBaseDomainBC* viscDomBC = (ConductivityBaseDomainBC*) m_domainBCFactory->create(m_eblgs[ref].getDomain(),
                                                                                              m_eblgs[ref].getEBISL(), dxMGLevel*RealVect::Unit);
  RefCountedPtr<ConductivityBaseEBBC>      ebbc(viscEBBC);
  RefCountedPtr<ConductivityBaseDomainBC> dombc(viscDomBC);

  bool hasFine = (ref < (m_eblgs.size()-1));
  bool hasCoar = (ref > 0);
  //optimization hook.  need to store the result out of EBArith::getCoarserLayoutss
  bool layoutChanged = true;

  NWOEBConductivityOp* newOp = NULL;

  newOp = new NWOEBConductivityOp(eblgFine, eblgMGLevel, eblgCoar, eblgCoarMG, m_quadCFI[ref],
                               dombc, ebbc,  dxMGLevel,dxCoar, refToFiner, refToCoarser,
                               hasFine, hasCoar, hasCoarMGObjects,  layoutChanged,
                               m_alpha, m_beta, m_acoef[ref], m_bcoef[ref], m_bcoefIrreg[ref],
                               m_ghostCellsPhi, m_ghostCellsRhs, m_relaxType);

  return newOp;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
NWOEBConductivityOpFactory::
reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim)
{
  delete a_reclaim;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
NWOEBConductivityOpFactory::
AMRreclaim(NWOEBConductivityOp* a_reclaim)
{
  delete a_reclaim;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
