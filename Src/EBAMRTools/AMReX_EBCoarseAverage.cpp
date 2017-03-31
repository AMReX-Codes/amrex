#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri aug 17, 2001
//major rework 4/2009 dtg

#include "EBCoarseAverage.H"
#include "EBAverageF_F.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "EBLoadBalance.H"
#include "NamespaceHeader.H"

/************************************/
void
EBCoarseAverage::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBCoarseAverage::EBCoarseAverage()
{
  setDefaultValues();
}
/************************************/
EBCoarseAverage::~EBCoarseAverage()
{
}
/************************************/
EBCoarseAverage::EBCoarseAverage(const DisjointBoxLayout& a_dblFine,
                                 const DisjointBoxLayout& a_dblCoar,
                                 const EBISLayout& a_ebislFine,
                                 const EBISLayout& a_ebislCoar,
                                 const ProblemDomain& a_domainCoar,
                                 const int& a_nref,
                                 const int& a_nvar,
                                 const EBIndexSpace* ebisPtr)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, ebisPtr);
}
/************************************/
EBCoarseAverage::EBCoarseAverage(const EBLevelGrid& a_eblgFine,
                                 const EBLevelGrid& a_eblgCoar,
                                 const EBLevelGrid& a_eblgCoFi,
                                 const int& a_nref,
                                 const int& a_nvar)
{
  setDefaultValues();

  define(a_eblgFine, a_eblgCoar, a_eblgCoFi, a_nref, a_nvar);
}
/************************************/
void
EBCoarseAverage::define(const DisjointBoxLayout& a_dblFine,
                        const DisjointBoxLayout& a_dblCoar,
                        const EBISLayout& a_ebislFine,
                        const EBISLayout& a_ebislCoar,
                        const ProblemDomain& a_domainCoar,
                        const int& a_nref,
                        const int& a_nvar,
                        const EBIndexSpace* ebisPtr)
{
  CH_TIME("EBCoarseAverage::define");
  CH_assert(ebisPtr->isDefined());

  ProblemDomain domainFine = a_domainCoar;
  domainFine.refine(a_nref);
  EBLevelGrid eblgFine;
  EBLevelGrid eblgCoar = EBLevelGrid(a_dblCoar, a_ebislCoar, a_domainCoar);
  EBLevelGrid eblgCoFi;

  //check to see if the input layout is coarsenable.
  //if so, proceed with ordinary drill
  //otherwise, see if the layout covers the domain.
  //if it does, we can use domainsplit
  if (a_dblFine.coarsenable(a_nref))
    {
      eblgFine = EBLevelGrid(a_dblFine, a_ebislFine,   domainFine);
      m_useFineBuffer = false;
    }
  else
    {
      Box fineDomBox = refine(a_domainCoar.domainBox(), a_nref);
      int numPtsDom = fineDomBox.numPts();
      //no need for gathers here because the meta data is global
      int numPtsLayout = 0;
      for (LayoutIterator lit = a_dblFine.layoutIterator(); lit.ok(); ++lit)
        {
          numPtsLayout += a_dblFine.get(lit()).numPts();
        }
      bool coveringDomain = (numPtsDom == numPtsLayout);
      if (coveringDomain)
        {
          m_useFineBuffer = true;
          int maxBoxSize = 4*a_nref;
          Vector<Box> boxes;
          Vector<int> procs;
          domainSplit(fineDomBox, boxes, maxBoxSize);
          mortonOrdering(boxes);
          LoadBalance(procs, boxes);
          DisjointBoxLayout dblBufFine(boxes, procs);

          eblgFine = EBLevelGrid(dblBufFine, domainFine, 2, eblgCoar.getEBIS());
        }
      else
        {
          pout() << "EBCoarseAverage::input layout is not coarsenable and does not cover the domain--bailing out" << endl;
          MayDay::Error();
        }
    }

  coarsen(eblgCoFi, eblgFine, a_nref);
  define(eblgFine, eblgCoar, eblgCoFi, a_nref, a_nvar);
}
/************************************/
void
EBCoarseAverage::define(const EBLevelGrid& a_eblgFine,
                        const EBLevelGrid& a_eblgCoar,
                        const EBLevelGrid& a_eblgCoFi,
                        const int& a_nref,
                        const int& a_nvar)
{
  CH_TIME("EBCoarseAverage::define(EBLG)");
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);
  CH_assert(a_eblgFine.getEBISL().getGhost() >= 2);
  CH_assert(a_eblgCoar.getEBISL().getGhost() >= 2);


  m_isDefined = true;
  m_refRat    = a_nref;
  m_nComp     = a_nvar;
  m_eblgCoar  = a_eblgCoar;

  m_eblgFine  = a_eblgFine;
  m_eblgCoFi  = a_eblgCoFi;
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_nref, m_eblgCoar.getEBIS());

  m_irregSetsCoFi.define(m_eblgCoFi.getDBL());
  for (DataIterator dit = m_eblgCoFi.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_irregSetsCoFi[dit()] = m_eblgCoFi.getEBISL()[dit()].getIrregIVS(m_eblgCoFi.getDBL().get(dit()));
    }

  if (m_useFineBuffer)
    {
      m_irregSetsFine.define(m_eblgFine.getDBL());
      for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          m_irregSetsFine[dit()] = m_eblgFine.getEBISL()[dit()].getIrregIVS(m_eblgFine.getDBL().get(dit()));
        }
    }
}
/************************************/
bool
EBCoarseAverage::isDefined() const
{
  return m_isDefined;
}
/************************************/
void
EBCoarseAverage::average(LevelData<EBCellFAB>& a_coarData,
                         const LevelData<EBCellFAB>& a_fineData,
                         const Interval& a_variables)
{
  CH_TIME("EBCoarseAverage::average(LD<EBCellFAB>)");
  Interval fineInterv, coarInterv;
  LevelData<EBCellFAB> coarFiData;
  LevelData<EBCellFAB> fineBuffer;
  CH_assert(isDefined());
  {
    CH_TIME("buffer allocation");
    EBCellFactory  factCoFi(m_eblgCoFi.getEBISL());
    coarFiData.define(  m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, factCoFi);
    if (m_useFineBuffer)
      {
        EBCellFactory    factFine(m_eblgFine.getEBISL());
        fineBuffer.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, factFine);
      }
  }

  if (m_useFineBuffer)
    {
      CH_TIME("copy_fine");
      a_fineData.copyTo(a_variables, fineBuffer, a_variables);
    }

  {
    CH_TIME("averaging");
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        const EBCellFAB* fineFAB = NULL;
        if (m_useFineBuffer)
          {
            fineInterv = Interval(0, a_variables.size()-1);
            fineFAB = &fineBuffer[dit()];
          }
        else
          {
            fineInterv = a_variables;
            fineFAB = &a_fineData[dit()];
          }

        coarInterv = Interval(0, a_variables.size()-1);
        averageFAB(coarFiData[dit()],
                   *fineFAB,
                   dit(),
                   fineInterv,
                   coarInterv);
      }
  }
  {
    CH_TIME("copy_coar");
    coarFiData.copyTo(coarInterv, a_coarData, fineInterv);
  }
}
/************************************/
void
EBCoarseAverage::average(LevelData<EBFluxFAB>&       a_coarData,
                         const LevelData<EBFluxFAB>& a_fineData,
                         const Interval& a_variables)
{
  CH_TIME("EBCoarseAverage::average(LD<EBFluxFAB>)");
  Interval fineInterv, coarInterv;
  LevelData<EBFluxFAB> coarFiData;
  LevelData<EBFluxFAB> fineBuffer;
  CH_assert(isDefined());
  {
    CH_TIME("buffer allocation");
    EBFluxFactory factCoFi(m_eblgCoFi.getEBISL());
    coarFiData.define(  m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, factCoFi);
    if (m_useFineBuffer)
      {
        EBFluxFactory    factFine(m_eblgFine.getEBISL());
        fineBuffer.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, factFine);
      }
  }

  if (m_useFineBuffer)
    {
      CH_TIME("fine_copy");
      a_fineData.copyTo(a_variables, fineBuffer, a_variables);
    }

  {
    CH_TIME("averaging");
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        const EBFluxFAB* fineFABPtr = NULL;
        CH_assert(a_variables.size() <= m_nComp);
        if (m_useFineBuffer)
          {
            fineInterv = Interval(0, a_variables.size()-1);
            fineFABPtr = &fineBuffer[dit()];
          }
        else
          {
            fineInterv = a_variables;
            fineFABPtr = &a_fineData[dit()];
          }
        EBFluxFAB&       cofiFAB = coarFiData[dit()];
        const EBFluxFAB& fineFAB = *fineFABPtr;
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            coarInterv = Interval(0, a_variables.size()-1);
            averageFAB(cofiFAB[idir],
                       fineFAB[idir],
                       dit(),
                       fineInterv,
                       coarInterv,
                       idir);

          }
      }
  }
  {
    CH_TIME("copy_coar");
    coarFiData.copyTo(coarInterv, a_coarData, fineInterv);
  }
}
/************************************/
void
EBCoarseAverage::average(LevelData<BaseIVFAB<Real> >&        a_coarData,
                         const LevelData<BaseIVFAB<Real> >&  a_fineData,
                         const Interval&                     a_variables)
{
  CH_TIME("EBCoarseAverage::average(LD<BaseIVFAB>)");
  LevelData<BaseIVFAB<Real> > coarFiData;
  LevelData<BaseIVFAB<Real> > fineBuffer;
  CH_assert(a_variables.begin() == 0);
  CH_assert(a_variables.size() <= m_nComp);
  CH_assert(isDefined());
  {
    CH_TIME("buffer allocation");
    BaseIVFactory<Real> factCoFi(m_eblgCoFi.getEBISL(), m_irregSetsCoFi);
    coarFiData.define(m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, factCoFi);
    if (m_useFineBuffer)
      {
        BaseIVFactory<Real> factFine(m_eblgFine.getEBISL(), m_irregSetsFine);
        coarFiData.define(m_eblgFine.getDBL(), m_nComp, IntVect::Zero, factFine);
      }
  }

  if (m_useFineBuffer)
  {
    CH_TIME("fine_copy");
    a_fineData.copyTo(a_variables, fineBuffer, a_variables);
  }

  {
    CH_TIME("averaging");
    int ifnerg = 0;
    for (DataIterator dit = m_eblgFine.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        const BaseIVFAB<Real>* fineFABPtr = NULL;
        if (m_useFineBuffer)
          {
            fineFABPtr = &fineBuffer[dit()];
          }
        else
          {
            fineFABPtr = &a_fineData[dit()];
          }
        BaseIVFAB<Real>&       cofiFAB = coarFiData[dit()];
        const BaseIVFAB<Real>& fineFAB = *fineFABPtr;
        averageFAB(cofiFAB,
                   fineFAB,
                   dit(),
                   a_variables);

        ifnerg++;
      }
  }
  {
    CH_TIME("copy_coar");
    coarFiData.copyTo(a_variables, a_coarData, a_variables);
  }
}
/************************************/
void
EBCoarseAverage::averageFAB(EBCellFAB&       a_coar,
                            const EBCellFAB& a_fine,
                            const DataIndex& a_datInd,
                            const Interval&  a_fineInterv,
                            const Interval&  a_coarInterv) const
{
  CH_TIME("EBCoarseAverage::averageFAB(EBCellFAB)");
  CH_assert(isDefined());
  //do all cells as if they were regular
  BaseFab<Real>& coarRegFAB =  a_coar.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();
  Box refbox(IntVect::Zero, (m_refRat-1)*IntVect::Unit);

  const Box& coarBox = m_eblgCoFi.getDBL().get(a_datInd);

#ifndef NDEBUG
  Box fineBox = refine(coarBox, m_refRat);
  CH_assert(coarRegFAB.box().contains(coarBox));
  CH_assert(fineRegFAB.box().contains(fineBox));
#endif

  for (int ioff = 0;
       ioff < a_fineInterv.size(); ioff++)
    {
      int ivarf = a_fineInterv.begin() + ioff;
      int ivarc = a_coarInterv.begin() + ioff;
      FORT_EBAVERAGE(CHF_FRA1(coarRegFAB,ivarc),
                     CHF_CONST_FRA1(fineRegFAB,ivarf),
                     CHF_BOX(coarBox),
                     CHF_CONST_INT(m_refRat),
                     CHF_BOX(refbox));
    }

  //overwrite irregular cells

  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];
  IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarBox);

  //fine cell volume is normalized to one.
  //compute coarse volume
  Real dxCoar = Real(m_refRat);
  Real cellVolCoar = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    cellVolCoar *= dxCoar;

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph());
      vofitCoar.ok(); ++vofitCoar)
    {
      const VolIndex& coarVoF = vofitCoar();
      Vector<VolIndex> fineVoFs =
        m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

      Real volCoar = cellVolCoar*ebisBoxCoar.volFrac(coarVoF);
      for (int ioff = 0;
           ioff < a_fineInterv.size(); ioff++)
        {
          int ivarf = a_fineInterv.begin() + ioff;
          int ivarc = a_coarInterv.begin() + ioff;
          Real dataVal = 0.0;
          if (volCoar > 0.)
            {
              for (int ifine = 0; ifine < fineVoFs.size(); ifine++)
                {
                  const VolIndex& fineVoF = fineVoFs[ifine];
                  //fine cell volume is normalized to one...
                  Real volFine  = ebisBoxFine.volFrac(fineVoF);
                  Real fineVal =  a_fine(fineVoF, ivarf);
                  if (volFine > 0.)
                    {
                      dataVal += fineVal*volFine;
                    }
                }
              dataVal /= volCoar;
            }
          else
            {
              //if there is no real volume, just take the ave
              //of fine values
              for (int ifine = 0; ifine < fineVoFs.size(); ifine++)
                {
                  const VolIndex& fineVoF = fineVoFs[ifine];
                  Real fineVal =  a_fine(fineVoF, ivarf);
                  dataVal += fineVal;
                }
              if (fineVoFs.size() > 0)
                {
                  dataVal /= Real(fineVoFs.size());
                }
            }
          a_coar(coarVoF, ivarc) = dataVal;
        }
    }
}
/***/
void
EBCoarseAverage::averageFAB(BaseIVFAB<Real>&       a_coar,
                            const BaseIVFAB<Real>& a_fine,
                            const DataIndex&       a_datInd,
                            const Interval&        a_variables) const
{
  CH_assert(isDefined());
  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];
  const IntVectSet& coarIrregIVS = a_coar.getIVS();
  const IntVectSet& fineIrregIVS = a_fine.getIVS();

  for (VoFIterator vofitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph());
      vofitCoar.ok(); ++vofitCoar)
    {
      const VolIndex& coarVoF = vofitCoar();
      Vector<VolIndex> fineVoFs =
        m_eblgCoFi.getEBISL().refine(coarVoF, m_refRat, a_datInd);

      for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
        {
          int  numVoFs = 0;
          Real areaTot = 0;
          Real dataVal = 0;
          for (int ifine = 0; ifine < fineVoFs.size(); ifine++)
            {
              const VolIndex& fineVoF = fineVoFs[ifine];
              if (fineIrregIVS.contains(fineVoF.gridIndex()))
                {
                  Real bndryArea = ebisBoxFine.bndryArea(fineVoF);
                  if (bndryArea > 0)
                    {
                      areaTot += bndryArea;
                      numVoFs++;
                      dataVal += bndryArea*a_fine(fineVoF, ivar);
                    }
                }
            }
          if (areaTot > 0)
            {
              dataVal /= areaTot;
            }
          a_coar(coarVoF, ivar) = dataVal;
        }
    }
}
/************************************/
void
EBCoarseAverage::averageFAB(EBFaceFAB&       a_coar,
                            const EBFaceFAB& a_fine,
                            const DataIndex& a_datInd,
                            const Interval&  a_fineInterv,
                            const Interval&  a_coarInterv,
                            const int&       a_dir) const
{
  CH_TIME("EBCoarseAverage::averageFAB(EBFaceFAB)");
  CH_assert(isDefined());

  //do all cells as if they were regular
  BaseFab<Real>& coarRegFAB =  a_coar.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();
  Box refbox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);
  refbox.surroundingNodes(a_dir);
  const Box& coarDBLBox = m_eblgCoFi.getDBL().get(a_datInd);
  Box coarFaceBox = coarDBLBox;
  coarFaceBox.surroundingNodes(a_dir);

#ifndef NDEBUG
  Box fineBox = refine(coarFaceBox, m_refRat);
  CH_assert(coarRegFAB.box().contains(coarFaceBox));
  CH_assert(fineRegFAB.box().contains(fineBox));
#endif

  CH_assert(a_fineInterv.size() == a_coarInterv.size());
  for (int ioff = 0;
       ioff < a_fineInterv.size(); ioff++)
    {
      int ivarf = a_fineInterv.begin() + ioff;
      int ivarc = a_coarInterv.begin() + ioff;
      FORT_EBAVERAGEFACE(CHF_FRA1(coarRegFAB,ivarc),
                         CHF_CONST_FRA1(fineRegFAB,ivarf),
                         CHF_BOX(coarFaceBox),
                         CHF_CONST_INT(m_refRat),
                         CHF_BOX(refbox),
                         CHF_CONST_INT(a_dir));
    }

  //overwrite irregular cells

  //recall that datInd is from the fine layout.
  const EBISBox& ebisBoxCoar = m_eblgCoFi.getEBISL()[a_datInd];
  const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[a_datInd];
  IntVectSet coarIrregIVS = ebisBoxCoar.getIrregIVS(coarDBLBox);

  //fine face area is normalized to one.
  //compute coarse volume
  Real dxCoar = Real(m_refRat);
  Real faceAreaCoar = 1.0;
  for (int idir = 0; idir < (SpaceDim - 1); idir++)
    faceAreaCoar *= dxCoar;

  for (FaceIterator faceitCoar(coarIrregIVS, ebisBoxCoar.getEBGraph(),a_dir,
                              FaceStop::SurroundingWithBoundary);
      faceitCoar.ok(); ++faceitCoar)
    {
      const FaceIndex& coarFace = faceitCoar();
      Vector<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(coarFace,m_refRat,a_datInd);

      Real areaCoar = faceAreaCoar*ebisBoxCoar.areaFrac(coarFace);
      for (int ioff = 0;
           ioff < a_fineInterv.size(); ioff++)
        {
          int ivarf = a_fineInterv.begin() + ioff;
          int ivarc = a_coarInterv.begin() + ioff;

          Real dataVal = 0.0;
          if (areaCoar > 0.)
            {
              for (int ifine = 0; ifine < fineFaces.size(); ifine++)
                {
                  const FaceIndex& fineFace = fineFaces[ifine];
                  //fine face area is normalized to one...
                  Real areaFine  = ebisBoxFine.areaFrac(fineFace);
                  Real fineVal =  a_fine(fineFace, ivarf);
                  if (areaFine > 0.)
                    {
                      dataVal += fineVal*areaFine;
                    }
                }
              dataVal /= areaCoar;
            }
          else
            {
              //if there is no real area, just take the ave
              //of fine values
              for (int ifine = 0; ifine < fineFaces.size(); ifine++)
                {
                  const FaceIndex& fineFace = fineFaces[ifine];
                  Real fineVal =  a_fine(fineFace, ivarf);
                  dataVal += fineVal;
                }
              if (fineFaces.size() > 0)
                {
                  dataVal /= Real(fineFaces.size());
                }
            }

          a_coar(coarFace, ivarc) = dataVal;
        }
    }
}


/************************************/
#include "NamespaceFooter.H"
