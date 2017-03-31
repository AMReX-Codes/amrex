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

#include "EBPWLFineInterp.H"
#include "EBInterpolateF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

/************************************/
/************************************/
void
EBPWLFineInterp::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
/************************************/
EBPWLFineInterp::EBPWLFineInterp()
{
  setDefaultValues();
}
/************************************/
/************************************/
EBPWLFineInterp::~EBPWLFineInterp()
{
}
/************************************/
/************************************/
EBPWLFineInterp::EBPWLFineInterp(const DisjointBoxLayout& a_dblFine,
                                 const DisjointBoxLayout& a_dblCoar,
                                 const EBISLayout& a_ebislFine,
                                 const EBISLayout& a_ebislCoar,
                                 const ProblemDomain& a_domainCoar,
                                 const int& a_nref,
                                 const int& a_nvar,
                                 const EBIndexSpace* const a_ebisPtr)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, a_ebisPtr);
}
/************************************/
/************************************/
void
EBPWLFineInterp::define(const DisjointBoxLayout& a_dblFine,
                        const DisjointBoxLayout& a_dblCoar,
                        const EBISLayout& a_ebislFine,
                        const EBISLayout& a_ebislCoar,
                        const ProblemDomain& a_domainCoar,
                        const int& a_nref,
                        const int& a_nvar,
                        const EBIndexSpace* const a_ebisPtr)
{
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);

  m_isDefined = true;
  m_refRat = a_nref;
  m_nComp = a_nvar;

  m_coarDomain = a_domainCoar;
  //need one ghost cell because stencil is one wide
  IntVect ghostiv = IntVect::Unit;
  int nghost = 1;

  CH_assert(a_ebisPtr->isDefined());
  m_coarsenedFineGrids = DisjointBoxLayout();
  coarsen(m_coarsenedFineGrids, a_dblFine, m_refRat);
  a_ebisPtr->fillEBISLayout(m_coarsenedFineEBISL,
                          m_coarsenedFineGrids,
                          a_domainCoar, nghost);
  m_coarsenedFineEBISL.setMaxRefinementRatio(m_refRat, a_ebisPtr);

  EBCellFactory ebcellfact(m_coarsenedFineEBISL);
  m_coarsenedFineData.define(m_coarsenedFineGrids, m_nComp,
                             ghostiv, ebcellfact);
  makeDerivStencils();
}
/************************************/
/************************************/
void
EBPWLFineInterp::makeDerivStencils()
{
  //define irregular regions to be irregular
  //cells plus all regular neighborsof
  //multivalued cells
  m_irregRegions.define(m_coarsenedFineGrids);
  for (DataIterator dit = m_coarsenedFineGrids.dataIterator();
      dit.ok(); ++dit)
    {
      const Box& localBox = m_coarsenedFineGrids.get(dit());
      const EBISBox& ebisBox  = m_coarsenedFineEBISL[dit()];
      IntVectSet&    localIrreg   = m_irregRegions[dit()];

      localIrreg = ebisBox.getMultiCells(localBox);
      //if there are any multivalued cells, grow the ivs and
      //intersect it with the local box and subtrac off covered boxes
      if (!localIrreg.isEmpty())
        {
          localIrreg.grow(1);
          localIrreg &= localBox;
          IntVectSet lo(localIrreg); // need to copy for iteration.  bvs
          IVSIterator ivsit(lo);
          for (ivsit.reset(); ivsit.ok(); ++ivsit)
            {
              const IntVect& iv = ivsit();
              if (ebisBox.isCovered(iv))
                {
                  localIrreg -= iv;
                }
            }
        }
      localIrreg |= ebisBox.getIrregIVS(localBox);
    }

  //now build the actual stencils
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_hiStencils[idir].define(m_coarsenedFineGrids);
      m_loStencils[idir].define(m_coarsenedFineGrids);

      LayoutData<BaseIVFAB<VoFStencil > >&
        lostenLDF= m_loStencils[idir];
      LayoutData<BaseIVFAB<VoFStencil > >&
        histenLDF= m_hiStencils[idir];

      for (DataIterator dit = m_coarsenedFineGrids.dataIterator();
          dit.ok(); ++dit)
        {
          const EBISBox& ebisBox  =  m_coarsenedFineEBISL[dit()];
          const IntVectSet& localIrreg  =  m_irregRegions[dit()];
          BaseIVFAB<VoFStencil >& histenFAB = histenLDF[dit()];
          BaseIVFAB<VoFStencil >& lostenFAB = lostenLDF[dit()];

          histenFAB.define(localIrreg, ebisBox.getEBGraph(), 1);
          lostenFAB.define(localIrreg, ebisBox.getEBGraph(), 1);

          for (VoFIterator vofit(localIrreg, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Vector<FaceIndex> loFaces=
                ebisBox.getFaces(vof, idir, Side::Lo);
              Vector<FaceIndex> hiFaces=
                ebisBox.getFaces(vof, idir, Side::Hi);
              Vector<VolIndex> loVoFs;
              Vector<VolIndex> hiVoFs;
              for (int iface = 0; iface < loFaces.size(); iface++)
                {
                  //use one-sided diff at edge of domain boundary
                  if (!loFaces[iface].isBoundary())
                    loVoFs.push_back(loFaces[iface].getVoF(Side::Lo));
                }
              for (int iface = 0; iface < hiFaces.size(); iface++)
                {
                  //use one-sided diff at edge of domain boundary
                  if (!hiFaces[iface].isBoundary())
                    hiVoFs.push_back(hiFaces[iface].getVoF(Side::Hi));
                }
              if (loVoFs.size() > 0)
                {
                  // have vofs on the low side.
                  //one sided diff with low side
                  Vector<VolIndex> stenVoFs;
                  Vector<Real>     stenWeig;

                  Real rfacesLo = Real(loVoFs.size());
                  Vector<Real> loWeig(loVoFs.size(), -1.0/rfacesLo);
                  stenVoFs = loVoFs;
                  stenWeig = loWeig;
                  stenVoFs.push_back(vof);
                  stenWeig.push_back(1.0);

                  //finally put the stencil into its container
                  VoFStencil& vofsten = lostenFAB(vof, 0);
                  vofsten.clear();
                  CH_assert(stenVoFs.size() == stenWeig.size());
                  for (int isten = 0; isten < stenVoFs.size(); isten++)
                    {
                      vofsten.add(stenVoFs[isten], stenWeig[isten]);
                    }
                }
              if (hiVoFs.size() > 0)
                {
                  // have vofs on the high side.
                  //one sided diff with high side
                  Vector<VolIndex> stenVoFs;
                  Vector<Real>     stenWeig;

                  Real rfacesHi = Real(hiVoFs.size());
                  Vector<Real> hiWeig(hiVoFs.size(),  1.0/rfacesHi);
                  stenVoFs = hiVoFs;
                  stenWeig = hiWeig;
                  stenVoFs.push_back(vof);
                  stenWeig.push_back(-1.0);

                  //finally put the stencil into its container
                  VoFStencil& vofsten = histenFAB(vof, 0);
                  vofsten.clear();
                  CH_assert(stenVoFs.size() == stenWeig.size());
                  for (int isten = 0; isten < stenVoFs.size(); isten++)
                    {
                      vofsten.add(stenVoFs[isten], stenWeig[isten]);
                    }
                }
            }
        }
    }
}
/************************************/
/************************************/
bool
EBPWLFineInterp::isDefined() const
{
  return m_isDefined;
}
/************************************/
/************************************/
void
EBPWLFineInterp::interpolate(LevelData<EBCellFAB>& a_fineData,
                             const LevelData<EBCellFAB>& a_coarData,
                             const Interval& a_variables)
{
  CH_assert(isDefined());
  a_coarData.copyTo(a_variables, m_coarsenedFineData, a_variables);

  for (DataIterator fineit = m_coarsenedFineGrids.dataIterator();
      fineit.ok(); ++fineit)
    {
      // interpolateFAB interpolates from an entire coarse grid onto an
      // entire fine grid.
      interpolateFAB(a_fineData[fineit()],
                     m_coarsenedFineData[fineit()],
                     fineit(),
                     a_variables);
    }

}
/************************************/
/************************************/
void
EBPWLFineInterp::interpolateFAB(EBCellFAB& a_fine,
                                const EBCellFAB& a_coar,
                                const DataIndex& a_datInd,
                                const Interval& a_variables) const
{
  const Box& bCoar = m_coarsenedFineGrids[a_datInd];
  Box refBox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);

  const BaseFab<Real>& coarRegFAB = a_coar.getSingleValuedFAB();
  BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();
  for (int ivar = a_variables.begin();
      ivar <= a_variables.end(); ivar++)
    {
      // fill fine data with piecewise constant coarse data
      FORT_EBINTERPCONSTANT(CHF_FRA1(fineRegFAB, ivar),
                            CHF_CONST_FRA1(coarRegFAB,ivar),
                            CHF_BOX(bCoar),
                            CHF_CONST_INT(m_refRat),
                            CHF_BOX(refBox));

      //add in minmod slope in each direction
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          BaseFab<Real> mmSlope(bCoar, 1);
          BaseFab<Real> loSlope;
          BaseFab<Real> hiSlope;
          mmSlope.setVal(0.0);

          Box interiorDomain = m_coarDomain.domainBox();
          Box loDiffDomain = m_coarDomain.domainBox();
          Box hiDiffDomain = m_coarDomain.domainBox();

          Box bInt = bCoar & interiorDomain;

          Box bLo = bCoar & loDiffDomain;
          Box bHi = bCoar & hiDiffDomain;

          if (!bLo.isEmpty())
            {
              loSlope.resize(bLo, 1);
              loSlope.setVal(0.0);
              FORT_EBLOSIDESLOPE(CHF_FRA1(loSlope, 0),
                                 CHF_CONST_FRA1(coarRegFAB, ivar),
                                 CHF_BOX(bLo),
                                 CHF_CONST_INT(idir));
              //copy over intersection.  this will get theslopes on
              //the edge of domain correct
              mmSlope.copy(loSlope);
            }
          if (!bHi.isEmpty())
            {
              hiSlope.resize(bHi, 1);
              hiSlope.setVal(0.0);
              FORT_EBHISIDESLOPE(CHF_FRA1(hiSlope, 0),
                                 CHF_CONST_FRA1(coarRegFAB, ivar),
                                 CHF_BOX(bHi),
                                 CHF_CONST_INT(idir));
              //copy over intersection.  this will get theslopes on
              //the edge of domain correct
              mmSlope.copy(hiSlope);
            }

          //do maxminmod thing over interior of domain
          if (!bInt.isEmpty())
            {
              //figure out minmod slope in interior
              FORT_EBMAXMINMOD(CHF_FRA1(mmSlope, 0),
                               CHF_CONST_FRA1(loSlope, 0),
                               CHF_CONST_FRA1(hiSlope, 0),
                               CHF_BOX(bInt));
            }

          FORT_EBINTERPLINEAR(CHF_FRA1(fineRegFAB, ivar),
                              CHF_CONST_FRA1(mmSlope, 0),
                              CHF_BOX(bCoar),
                              CHF_CONST_INT(idir),
                              CHF_CONST_INT(m_refRat),
                              CHF_BOX(refBox));

        } //end loop over directions

      //now do irregular cells.
      const EBISBox& coarEBISBox = m_coarsenedFineEBISL[a_datInd];
      for (VoFIterator vofit(m_irregRegions[a_datInd], coarEBISBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& coarVoF = vofit();
          const IntVect& coarIV = coarVoF.gridIndex();
          Vector<VolIndex> fineVoFs = m_coarsenedFineEBISL.refine(coarVoF, m_refRat, a_datInd);
          for (int ifinevof = 0; ifinevof < fineVoFs.size(); ifinevof++)
            {
              const VolIndex& fineVoF = fineVoFs[ifinevof];
              const IntVect& fineIV = fineVoF.gridIndex();

              //first set the fine stuff to the coar stuff
              Real fineVal = a_coar(coarVoF, ivar);

              //now add first derivative terms
              for (int idir = 0; idir  < SpaceDim; idir++)
                {

                  const LayoutData<BaseIVFAB<VoFStencil > >&
                    lodirSten = m_loStencils[idir];
                  const LayoutData<BaseIVFAB<VoFStencil > >&
                    hidirSten = m_hiStencils[idir];
                  const BaseIVFAB<VoFStencil >& lostenBF = lodirSten[a_datInd];
                  const BaseIVFAB<VoFStencil >& histenBF = hidirSten[a_datInd];

                  //compute delta idir direction (normalized s.t. dxcoar == 1)
                  Real deltahi = 0.0;
                  Real deltalo = 0.0;
                  bool hasHiSten;
                  bool hasLoSten;
                  {
                    const VoFStencil& hivofsten = histenBF(coarVoF,0);
                    for (int isten = 0; isten < hivofsten.size(); isten++)
                      {
                        const VolIndex&  stenvof = hivofsten.vof(isten);
                        const Real& stenwgt = hivofsten.weight(isten);
                        deltahi += stenwgt*a_coar(stenvof, ivar);
                      }
                    hasHiSten = (hivofsten.size() > 0);
                  }
                  {
                    const VoFStencil& lovofsten = lostenBF(coarVoF,0);
                    for (int isten = 0; isten < lovofsten.size(); isten++)
                      {
                        const VolIndex&  stenvof = lovofsten.vof(isten);
                        const Real& stenwgt = lovofsten.weight(isten);
                        deltalo += stenwgt*a_coar(stenvof, ivar);
                      }
                    hasLoSten = (lovofsten.size() > 0);
                  }

                  Real deltaminmod;
                  if (hasHiSten && hasLoSten)
                    {
                      //have deltas in both directions, do minmod thing.
                      Real mono = deltahi*deltalo;
                      if (mono > 0.0)
                        {
                          Real rsign = 1.0;
                          if ((deltahi + deltalo) < 0.0)
                            rsign = -1.0;
                          deltaminmod = rsign*Min(Abs(deltalo), Abs(deltahi));
                        }
                      else
                        {
                          deltaminmod = 0.0;
                        }
                    }
                  else if (hasHiSten)
                    {
                      deltaminmod = deltahi;
                    }
                  else if (hasLoSten)
                    {
                      deltaminmod = deltalo;
                    }
                  else
                    {
                      //no derivs exist.
                      deltaminmod = 0.0;
                    }
                  //these are the locations in space of the data
                  Real coarLoc = Real(coarIV[idir]) + 0.5;
                  //the factor of 1/refrat is because we normalized
                  //for dxcoar == 1
                  Real fineLoc  = (Real(fineIV[idir]) + 0.5)/Real(m_refRat);;
                  //add in contribution of deriv*dist
                  fineVal += deltaminmod*(fineLoc-coarLoc);
                }
              a_fine(fineVoF, ivar) = fineVal;
            } //end loop over fine vofs = refine(coarvof)
        } //end loop over coarse vofs

    }//end loop over variables
}
/************************************/
/************************************/
#include "NamespaceFooter.H"
