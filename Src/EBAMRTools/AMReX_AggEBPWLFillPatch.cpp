#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves wednesday, june 1, 2011

#include "AggEBPWLFillPatch.H"
#include "EBInterpolateF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "LayoutIterator.H"
#include "EBIndexSpace.H"
#include "EBAlias.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

/************************************/
/************************************/
void
AggEBPWLFillPatch::
setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
  m_radius = -1;
  m_patcher = NULL;
}
/************************************/
/************************************/
AggEBPWLFillPatch::
AggEBPWLFillPatch()
{
  setDefaultValues();
}
/************************************/
/************************************/
AggEBPWLFillPatch::
~AggEBPWLFillPatch()
{
  if (m_patcher != NULL)
    {
      delete m_patcher;
    }
}
/************************************/
/************************************/
AggEBPWLFillPatch::
AggEBPWLFillPatch(const DisjointBoxLayout& a_dblFine,
                  const DisjointBoxLayout& a_dblCoar,
                  const EBISLayout& a_ebislFine,
                  const EBISLayout& a_ebislCoar,
                  const ProblemDomain& a_domainCoar,
                  const int& a_nref,
                  const int& a_nvar,
                  const int& a_radius,
                  const IntVect& a_ghost,
                  const bool& a_forceNoEBCF,
                  const EBIndexSpace* const a_ebisPtr)
{
  CH_TIME("AggEBPWLFillPatch::AggEBPWLFillPatch");
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, a_radius,  a_ghost, 
         a_forceNoEBCF, a_ebisPtr);
}
/************************************/
/************************************/
void
AggEBPWLFillPatch::
define(const DisjointBoxLayout& a_dblFine,
       const DisjointBoxLayout& a_dblCoar,
       const EBISLayout& a_ebislFine,
       const EBISLayout& a_ebislCoar,
       const ProblemDomain& a_domainCoar,
       const int& a_nref,
       const int& a_nvar,
       const int& a_radius,
       const IntVect& a_ghost,
       const bool& a_forceNoEBCF,
       const EBIndexSpace* const a_ebisPtr)
{
  CH_TIME("AggEBPWLFillPatch::define");
  CH_assert(a_nref > 0);                               //
  CH_assert(a_nvar > 0);
  //first cell is interpolated, rest are extrapolated.
  //if this is 0, there is not even interpolation.
  CH_assert(a_radius > 0);
  m_forceNoEBCF = a_forceNoEBCF;
  m_radius = a_radius;
  m_isDefined = true;
  m_refRat = a_nref;
  m_nComp = a_nvar;
  m_ghost = a_ghost;
  m_coarDomain = a_domainCoar;

  //setup the non-EB PWLFillPatch
  definePieceWiseLinearFillPatch(a_dblFine,
                                 a_dblCoar);

  if(!m_forceNoEBCF)
    {
      m_fineEBISL = a_ebislFine;
      if (m_fineEBISL.getMaxCoarseningRatio() < m_refRat)
        {
          m_fineEBISL.setMaxCoarseningRatio(m_refRat,a_ebisPtr);
        }
      m_fineGrids = a_dblFine;
      m_coarGrids = a_dblCoar;
      //needs to be extra because the stencil of
      //the derivs has an extra coarse cell
      const int coarse_slope_radius =
        (m_radius + m_refRat - 1) / m_refRat;
      const int coarse_ghost_radius = coarse_slope_radius + 1;

      m_coarGhostRad = coarse_ghost_radius;
      IntVect ghostiv = m_coarGhostRad*IntVect::Unit;

      CH_assert(a_ebisPtr->isDefined());
      m_coarsenedFineGrids = DisjointBoxLayout();
      coarsen(m_coarsenedFineGrids, a_dblFine, m_refRat);
      a_ebisPtr->fillEBISLayout(m_coarsenedFineEBISL,
                                m_coarsenedFineGrids,
                                a_domainCoar, m_coarGhostRad);
      //  m_coarsenedFineEBISL.setMaxCoarseningRatio(m_refRat);
      EBCellFactory ebcellfact(m_coarsenedFineEBISL);
      m_coarOnFDataOld.define(m_coarsenedFineGrids, m_nComp,
                              ghostiv, ebcellfact);
      m_coarOnFDataNew.define(m_coarsenedFineGrids, m_nComp,
                              ghostiv, ebcellfact);
      makeStencils();
    }
}
/************************************/
/************************************/
void
AggEBPWLFillPatch::
makeStencils()
{
  CH_TIME("AggEBPWLFillPatch::makeStencil");
  CH_assert(!m_forceNoEBCF);
  LayoutData< IntVectSet > irregRegionsFine;
  LayoutData< IntVectSet > irregRegionsCoFi;
  LayoutData< Vector<VolIndex> > srcVoFs;
  LayoutData< IntVectSet >  coarCeInterp[SpaceDim];
  LayoutData< IntVectSet >  coarLoInterp[SpaceDim];
  LayoutData< IntVectSet >  coarHiInterp[SpaceDim];
  LayoutData< Vector <VoFStencil> >  hiStencils[SpaceDim];
  LayoutData< Vector <VoFStencil> >  loStencils[SpaceDim];
  irregRegionsFine.define(m_fineGrids);
  irregRegionsCoFi.define(m_fineGrids);
  srcVoFs.define(m_fineGrids);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      coarCeInterp[idir].define(m_coarsenedFineGrids);
      coarLoInterp[idir].define(m_coarsenedFineGrids);
      coarHiInterp[idir].define(m_coarsenedFineGrids);
    }

  getIVS(irregRegionsFine, irregRegionsCoFi, srcVoFs);

  getLoHiCenIVS(coarLoInterp,  coarHiInterp,  coarCeInterp);

  getSten(loStencils, hiStencils,
          coarLoInterp,  coarHiInterp, coarCeInterp,
          srcVoFs);


  defineSlopeHolders(irregRegionsCoFi);

  getOffsets(srcVoFs, irregRegionsFine, loStencils, hiStencils,
             coarLoInterp,  coarHiInterp, coarCeInterp);

  defineAggStencils(loStencils, hiStencils, srcVoFs);
}
/*********/
void
AggEBPWLFillPatch::
defineAggStencils(LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
                  LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
                  const LayoutData< Vector<VolIndex> >&   a_srcVoFs)

{
  CH_TIME("AggEBPWLFillPatch::defineAggStencils");
  CH_assert(!m_forceNoEBCF);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_stenLo[idir].define(m_coarsenedFineGrids);
      m_stenHi[idir].define(m_coarsenedFineGrids);

      for (DataIterator dit= m_coarsenedFineGrids.dataIterator(); dit.ok(); ++dit)
        {
          Vector<RefCountedPtr<BaseIndex   > > baseindice(a_srcVoFs[dit()].size());
          Vector<RefCountedPtr<BaseStencil > > basestenlo(a_srcVoFs[dit()].size());
          Vector<RefCountedPtr<BaseStencil > > basestenhi(a_srcVoFs[dit()].size());

          for (int ivof= 0; ivof < a_srcVoFs[dit()].size(); ivof++)
            {
              baseindice[ivof] =   RefCountedPtr<BaseIndex>(new   VolIndex(         a_srcVoFs[dit()][ivof]));
              basestenlo[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(a_loStencils[idir][dit()][ivof]));
              basestenhi[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(a_hiStencils[idir][dit()][ivof]));
            }
          //all the slopes are the same size so we can use any of them really
          m_stenLo[idir][dit()] = RefCountedPtr< AggStencil <EBCellFAB, BaseIVFAB<Real> > >
            (new AggStencil <EBCellFAB, BaseIVFAB<Real> >(baseindice, basestenlo, m_coarOnFDataOld[dit()], m_slopeLoOld[idir][dit()]));
          m_stenHi[idir][dit()] = RefCountedPtr< AggStencil <EBCellFAB, BaseIVFAB<Real> > >
            (new AggStencil <EBCellFAB, BaseIVFAB<Real> >(baseindice, basestenhi, m_coarOnFDataOld[dit()], m_slopeHiOld[idir][dit()]));
        }
    }
}
/*********/
void
AggEBPWLFillPatch::
defineSlopeHolders(const LayoutData<IntVectSet>& a_irregRegionsCoFi)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_slopeLoOld[idir].define(m_coarsenedFineGrids);
      m_slopeHiOld[idir].define(m_coarsenedFineGrids);
      m_slopeCeOld[idir].define(m_coarsenedFineGrids);
      m_slopeLoNew[idir].define(m_coarsenedFineGrids);
      m_slopeHiNew[idir].define(m_coarsenedFineGrids);
      m_slopeCeNew[idir].define(m_coarsenedFineGrids);

      for (DataIterator dit = m_coarsenedFineGrids.dataIterator(); dit.ok(); ++dit)
        {
          const EBGraph& ebgraph = m_coarsenedFineEBISL[dit()].getEBGraph();
          const IntVectSet& ivs = a_irregRegionsCoFi[dit()];
          m_slopeLoOld[idir][dit()].define(ivs, ebgraph, m_nComp);
          m_slopeHiOld[idir][dit()].define(ivs, ebgraph, m_nComp);
          m_slopeCeOld[idir][dit()].define(ivs, ebgraph, m_nComp);
          m_slopeLoNew[idir][dit()].define(ivs, ebgraph, m_nComp);
          m_slopeHiNew[idir][dit()].define(ivs, ebgraph, m_nComp);
          m_slopeCeNew[idir][dit()].define(ivs, ebgraph, m_nComp);
        }
    }
}
/************************************/
void
AggEBPWLFillPatch::
getOffsets(const LayoutData< Vector<VolIndex> >&  a_srcVoFsCoar,
           const LayoutData<IntVectSet>        &  a_irregRegionsFine,
           const LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
           const LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
           const LayoutData<IntVectSet>           a_coarLoInte[SpaceDim],
           const LayoutData<IntVectSet>           a_coarHiInte[SpaceDim],
           const LayoutData<IntVectSet>           a_coarCeInte[SpaceDim])

{
  CH_assert(!m_forceNoEBCF);
  //first get the coarse offsets for the slopes and so on
  m_coarOffsets.define(m_coarsenedFineGrids);
  int ibox = 0;
  for (DataIterator dit = m_coarsenedFineGrids.dataIterator(); dit.ok(); ++dit, ++ibox)
    {
      m_coarOffsets[dit()].resize(a_srcVoFsCoar[dit()].size());
      for (int ivof = 0; ivof < a_srcVoFsCoar[dit()].size(); ivof++)
        {
          //all the slope containers are the same size so I can
          //use any of them
          const VolIndex& vof = a_srcVoFsCoar[dit()][ivof];
          m_coarOffsets[dit()][ivof].slop_access.offset =     m_slopeLoOld[0][dit()].offset(  vof, 0);
          m_coarOffsets[dit()][ivof].slop_access.dataID =     m_slopeLoOld[0][dit()].dataType(vof);
          m_coarOffsets[dit()][ivof].data_access.offset =    m_coarOnFDataOld[dit()].offset(  vof, 0);
          m_coarOffsets[dit()][ivof].data_access.dataID =    m_coarOnFDataOld[dit()].dataType(vof);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              m_coarOffsets[dit()][ivof].has_lo[idir] = (a_loStencils[idir][dit()][ivof].size() != 0);
              m_coarOffsets[dit()][ivof].has_hi[idir] = (a_hiStencils[idir][dit()][ivof].size() != 0);

              m_coarOffsets[dit()][ivof].ivs_lo[idir] = (a_coarLoInte[idir][dit()].contains(vof.gridIndex()));
              m_coarOffsets[dit()][ivof].ivs_hi[idir] = (a_coarHiInte[idir][dit()].contains(vof.gridIndex()));
              m_coarOffsets[dit()][ivof].ivs_ce[idir] = (a_coarCeInte[idir][dit()].contains(vof.gridIndex()));
            }
        }
    }

  //now get the offsets for the fine data so we can put the answers in
  //the proper place.  the dummy data holder is so we can get the
  //offsets
  m_fineOffsets.define(    m_fineGrids);
  EBCellFactory ebcellfact(m_fineEBISL);
  LevelData<EBCellFAB> dummyfine(m_fineGrids, 1, m_ghost, ebcellfact);
  for (DataIterator dit = m_fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      const IntVectSet& ivs = a_irregRegionsFine[dit()];
      VoFIterator vofitfine(ivs, m_fineEBISL[dit()].getEBGraph());
      const Vector<VolIndex>& vofVec = vofitfine.getVector();
      m_fineOffsets[dit()].resize(vofVec.size());
      for (int ivof = 0; ivof < vofVec.size(); ivof++)
        {
          const VolIndex& vof = vofVec[ivof];
          m_fineOffsets[dit()][ivof].dest_access.offset = dummyfine[dit()].offset(  vof, 0);
          m_fineOffsets[dit()][ivof].dest_access.dataID = dummyfine[dit()].dataType(vof);
          //now find which coarse VoF this find vof coarsens to
          //(so we can use the slopes there)
          int coarind = -1;
          bool found = false;
          VolIndex coarVoF =  m_fineEBISL.coarsen(vof, m_refRat, dit());
          for (int icoar = 0; icoar < a_srcVoFsCoar[dit()].size(); icoar++)
            {
              const VolIndex& listVoF = a_srcVoFsCoar[dit()][icoar];
              if (listVoF == coarVoF)
                {
                  found = true;
                  coarind = icoar;
                  break;
                }
            }
          if (!found)
            {
              MayDay::Error("coarse-fine inconsistency in AggEBPWLFP");
            }

          m_fineOffsets[dit()][ivof].slop_index = coarind;
          m_fineOffsets[dit()][ivof].coariv = coarVoF.gridIndex();
          m_fineOffsets[dit()][ivof].fineiv =     vof.gridIndex();

        }
    }

}
/************************************/
void
AggEBPWLFillPatch::
getIVS(LayoutData<IntVectSet>& a_irregRegionsFine,
       LayoutData<IntVectSet>& a_irregRegionsCoar,
       LayoutData< Vector<VolIndex> >&   a_srcVoFs)

{
  CH_TIME("AggEBPWLFillPatch::getIVS");
  ProblemDomain domFine = refine(m_coarDomain, m_refRat);
  CH_assert(!m_forceNoEBCF);
  for (DataIterator dit = m_coarsenedFineGrids.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet&    localIrregF = a_irregRegionsFine[dit()];
      const EBISBox&   ebisBoxCF = m_coarsenedFineEBISL[dit()];
      const Box&        regionCF = ebisBoxCF.getRegion();
      //NOTE::this needs to be done for ALL boxes that are "near" the EB interface.
      //If one uses an isRegular check (as we did before) for this, you end up with inconsistent (yet accurate)
      //interpolations. Inconsistent (non-EB vs EB) interpolations cause problems for some codes that
      //rely on the fact that ghost cells are filled in exactly the same way.
      //This was a problem for the mac-projections because the advection algorithms' face velocities depend on
      //fillpatched ghost cells that must contain the same data where ghost cells overlap from one box to the next.

      //make localIrreg for all c/f ghost ivs within m_radius
      Box bigBox= m_fineGrids.get(dit());
      bigBox.grow(m_radius);
      bigBox &= domFine;
      localIrregF = IntVectSet(bigBox);
      for (LayoutIterator lit = m_fineGrids.layoutIterator();
          lit.ok(); ++lit)
        {
          localIrregF -= m_fineGrids.get(lit());
        }

      //keep only fine ivs that are within 1 coarse radius of irregular coarse ivs
      //the remainder of fine ivs are done with the non-EB FillPatch
      IntVectSet  irregCF = ebisBoxCF.getIrregIVS(regionCF);
      irregCF.grow(1);

      IntVectSet irregRefinedCF = irregCF;
      irregRefinedCF.refine(m_refRat);
      localIrregF &= irregRefinedCF;

      IntVectSet irregCoar = localIrregF;
      irregCoar.coarsen(m_refRat);

      a_irregRegionsCoar[dit()] = irregCoar;
      VoFIterator vofit(irregCoar, ebisBoxCF.getEBGraph());
      a_srcVoFs[dit()] = vofit.getVector();

    }
}

/************************************/
/************************************/
void
AggEBPWLFillPatch::
getLoHiCenIVS(LayoutData<IntVectSet>  a_coarLoInterp[SpaceDim],
              LayoutData<IntVectSet>  a_coarHiInterp[SpaceDim],
              LayoutData<IntVectSet>  a_coarCeInterp[SpaceDim])
{
  CH_TIME("AggEBPWLFillPatch::getLoHiCenIVS");
  //now create the coarse lo hi and center intvectsets--
  //these sets determine where you have to do one sided differences
  //because of the configuration of the coarse layout.
  CH_assert(!m_forceNoEBCF);
  int ibox = 0;
  for (DataIterator dit = m_coarsenedFineGrids.dataIterator();
      dit.ok(); ++dit)
    {
      const Box& fineBox = m_fineGrids.get(dit());
      Box coarsenedFineBox = coarsen(grow(fineBox, m_radius), m_refRat);
      coarsenedFineBox &= m_coarDomain;

      IntVectSet coarsenedFineInterp(coarsenedFineBox);

      // Iterate over boxes in coarsened fine domain, and subtract off
      // from the set of coarse cells from which the fine ghost cells
      // will be interpolated.

      for (LayoutIterator litCF = m_coarsenedFineGrids.layoutIterator();
          litCF.ok(); ++litCF)
        {
          const Box& otherCoarsenedBox = m_coarsenedFineGrids.get(litCF());
          coarsenedFineInterp -= otherCoarsenedBox;
        }

      // Now that we have the coarsened cells required for interpolation,
      // construct IntvectSets specifying the one-sided and centered
      // stencil locations.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet& coarseCeInterp = a_coarCeInterp[idir][dit()];
          IntVectSet& coarseLoInterp = a_coarLoInterp[idir][dit()];
          IntVectSet& coarseHiInterp = a_coarHiInterp[idir][dit()];

          coarseCeInterp = coarsenedFineInterp;

          coarseLoInterp = coarseCeInterp;
          coarseLoInterp.shift(BASISV(idir));

          coarseHiInterp = coarseCeInterp;
          coarseHiInterp.shift(-BASISV(idir));

          // We iterate over the coarse grids and subtract them off of the
          // one-sided stencils.
          for (LayoutIterator litC = m_coarGrids.layoutIterator(); litC.ok(); ++litC)
            {
              const Box& bx = m_coarGrids.get(litC());
              coarseLoInterp -= bx;
              coarseHiInterp -= bx;
            }

          coarseLoInterp.shift(-BASISV(idir));
          coarseHiInterp.shift(BASISV(idir));

          coarseCeInterp -= coarseLoInterp;
          coarseCeInterp -= coarseHiInterp;

        }
      ibox++;
    }
}

/************************************/
/************************************/
void
AggEBPWLFillPatch::
getSten(LayoutData<Vector<VoFStencil> >  a_loStencils[SpaceDim],
        LayoutData<Vector<VoFStencil> >  a_hiStencils[SpaceDim],
        LayoutData<IntVectSet>           a_coarLoInterp[SpaceDim],
        LayoutData<IntVectSet>           a_coarHiInterp[SpaceDim],
        LayoutData<IntVectSet>           a_coarCeInterp[SpaceDim],
        const LayoutData< Vector<VolIndex> >&   a_srcVoFs)
{
  CH_TIME("AggEBPWLFillPatch::getSten");               //
  CH_assert(!m_forceNoEBCF);

  /////////////////////////
  for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      a_hiStencils[derivDir].define(m_coarsenedFineGrids);
      a_loStencils[derivDir].define(m_coarsenedFineGrids);

      int ibox=0;
      for (DataIterator dit = m_coarsenedFineGrids.dataIterator();
          dit.ok(); ++dit, ++ibox)
        {
          const EBISBox& ebisBox  =  m_coarsenedFineEBISL[dit()];
          a_loStencils[derivDir][dit()].resize(a_srcVoFs[dit()].size());
          a_hiStencils[derivDir][dit()].resize(a_srcVoFs[dit()].size());

          for (int ivof = 0; ivof < a_srcVoFs[dit()].size(); ivof++)
            {
              const VolIndex& vof = a_srcVoFs[dit()][ivof];
              Vector<FaceIndex> loFaces=
                ebisBox.getFaces(vof, derivDir, Side::Lo);
              Vector<FaceIndex> hiFaces=
                ebisBox.getFaces(vof, derivDir, Side::Hi);
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
                  a_loStencils[derivDir][dit()][ivof].clear();
                  CH_assert(stenVoFs.size() == stenWeig.size());
                  for (int isten = 0; isten < stenVoFs.size(); isten++)
                    {
                      a_loStencils[derivDir][dit()][ivof].add(stenVoFs[isten], stenWeig[isten]);
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
                  a_hiStencils[derivDir][dit()][ivof].clear();
                  CH_assert(stenVoFs.size() == stenWeig.size());
                  for (int isten = 0; isten < stenVoFs.size(); isten++)
                    {
                      a_hiStencils[derivDir][dit()][ivof].add(stenVoFs[isten], stenWeig[isten]);
                    }
                }
            } //end loop over vofs
        } //end loop over grids
    } //end loop over derivative directions
}
/************************************/
void
AggEBPWLFillPatch::
interpolate(LevelData<EBCellFAB>& a_fineData,
            const LevelData<EBCellFAB>& a_coarDataOld,
            const LevelData<EBCellFAB>& a_coarDataNew,
            const Real& a_coarTimeOld,
            const Real& a_coarTimeNew,
            const Real& a_fineTime,
            const Interval& a_variables) const
{
  CH_TIME("AggEBPWLFillPatch::interpolate");
  CH_assert(isDefined());
  CH_assert(a_coarTimeNew >= a_coarTimeOld);
  CH_assert(a_coarTimeNew >= a_fineTime);
  CH_assert(a_fineTime >= a_coarTimeOld);
  CH_assert(a_fineData.ghostVect() == m_ghost);

  //do non-EB PWLFillPatch
  {
    LevelData<FArrayBox> fineDataLDFAB, coarOldDataLDFAB,coarNewDataLDFAB;
    aliasEB(fineDataLDFAB, a_fineData);
    aliasEB(coarOldDataLDFAB, (LevelData<EBCellFAB>&)a_coarDataOld);
    aliasEB(coarNewDataLDFAB, (LevelData<EBCellFAB>&)a_coarDataNew);
    Real timeCoef = 0.0;
    if (a_fineTime>a_coarTimeOld)
      {
        timeCoef = (a_fineTime - a_coarTimeOld)/(a_coarTimeNew - a_coarTimeOld);
      }

    m_patcher->fillInterp(fineDataLDFAB,
                          coarOldDataLDFAB,
                          coarNewDataLDFAB,
                          timeCoef,
                          a_variables.begin(),
                          a_variables.begin(),
                          a_variables.size());
  }

  if(!m_forceNoEBCF)
    {

      //level data copy fills coarse ghost cells as well as interior
      a_coarDataOld.copyTo(a_variables, m_coarOnFDataOld, a_variables);
      a_coarDataNew.copyTo(a_variables, m_coarOnFDataNew, a_variables);

      m_coarOnFDataOld.exchange(a_variables);
      m_coarOnFDataNew.exchange(a_variables);
      int ibox = 0;
      for (DataIterator fineit = m_coarsenedFineGrids.dataIterator();
           fineit.ok(); ++fineit)
        {
          Box fineBox = m_fineGrids.get(fineit());
          Box coarsenedFineBox = m_coarsenedFineGrids.get(fineit());
          EBCellFAB&  fineData = a_fineData[fineit()];
          const EBCellFAB&  coarDataOld = m_coarOnFDataOld[fineit()];
          const EBCellFAB&  coarDataNew = m_coarOnFDataNew[fineit()];
          // interpolateFAB interpolates from an entire coarse grid onto an
          // entire fine grids coarse-fine regions.
          interpolateFAB(fineData,
                         coarDataOld,
                         coarDataNew,
                         a_coarTimeOld,
                         a_coarTimeNew,
                         a_fineTime, fineit(),
                         a_variables);
          ibox++;
        }
    }
}
/************************************/
void
AggEBPWLFillPatch::
deltaMinMod(Real& deltaminmod, Real& deltalo, Real& deltahi) const
{
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
/************************************/
void
AggEBPWLFillPatch::
getSlopes(const EBCellFAB& a_coarDataOld,
          const EBCellFAB& a_coarDataNew,
          const DataIndex& a_dit,
          const Interval& a_variables) const
{
  CH_TIME("AggEBPWLFillPatch::getSlopes");
  CH_assert(!m_forceNoEBCF);
  int ibeg = a_variables.begin();  int isiz = a_variables.size();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //get low and high versions of the slopes by applying the stencils
      m_stenLo[idir][a_dit]->apply(m_slopeLoOld[idir][a_dit], a_coarDataOld, ibeg, ibeg, isiz, false);
      m_stenLo[idir][a_dit]->apply(m_slopeLoNew[idir][a_dit], a_coarDataNew, ibeg, ibeg, isiz, false);
      m_stenHi[idir][a_dit]->apply(m_slopeHiOld[idir][a_dit], a_coarDataOld, ibeg, ibeg, isiz, false);
      m_stenHi[idir][a_dit]->apply(m_slopeHiNew[idir][a_dit], a_coarDataNew, ibeg, ibeg, isiz, false);

      //now do delta min mod logic to get centered slope
      for (int ivar = ibeg; ivar < isiz; ivar++)
        {
          //baseivfabs only have one data type
          CH_assert(m_slopeLoOld[idir][a_dit].numDataTypes() == 1);
          CH_assert(m_slopeHiOld[idir][a_dit].numDataTypes() == 1);
          CH_assert(m_slopeCeOld[idir][a_dit].numDataTypes() == 1);
          CH_assert(m_slopeLoNew[idir][a_dit].numDataTypes() == 1);
          CH_assert(m_slopeHiNew[idir][a_dit].numDataTypes() == 1);
          CH_assert(m_slopeCeNew[idir][a_dit].numDataTypes() == 1);


          Real* ptrLoOld = m_slopeLoOld[idir][a_dit].dataPtr(ivar);
          Real* ptrHiOld = m_slopeHiOld[idir][a_dit].dataPtr(ivar);
          Real* ptrCeOld = m_slopeCeOld[idir][a_dit].dataPtr(ivar);
          Real* ptrLoNew = m_slopeLoNew[idir][a_dit].dataPtr(ivar);
          Real* ptrHiNew = m_slopeHiNew[idir][a_dit].dataPtr(ivar);
          Real* ptrCeNew = m_slopeCeNew[idir][a_dit].dataPtr(ivar);
          for (int ivof = 0; ivof < m_coarOffsets[a_dit].size(); ivof++)
            {
              const size_t& offset = m_coarOffsets[a_dit][ivof].slop_access.offset;
              Real& sloLoOld =  *(ptrLoOld + offset);
              Real& sloHiOld =  *(ptrHiOld + offset);
              Real& sloCeOld =  *(ptrCeOld + offset);
              Real& sloLoNew =  *(ptrLoNew + offset);
              Real& sloHiNew =  *(ptrHiNew + offset);
              Real& sloCeNew =  *(ptrCeNew + offset);
              const bool& hasLo = m_coarOffsets[a_dit][ivof].has_lo[idir];
              const bool& hasHi = m_coarOffsets[a_dit][ivof].has_hi[idir];

              const bool& ivsLo = m_coarOffsets[a_dit][ivof].ivs_lo[idir];
              const bool& ivsHi = m_coarOffsets[a_dit][ivof].ivs_hi[idir];
              const bool& ivsCe = m_coarOffsets[a_dit][ivof].ivs_ce[idir];
              if (ivsCe)
                {
                  if (hasLo && hasHi)
                    {
                      deltaMinMod(sloCeOld, sloLoOld, sloHiOld);
                      deltaMinMod(sloCeNew, sloLoNew, sloHiNew);
                    }
                  else if (hasLo)
                    {
                      //we only have slopes on the low side
                      sloCeOld = sloLoOld;
                      sloCeNew = sloLoNew;
                    }
                  else if (hasHi)
                    {
                      //we only have slopes on the high side
                      sloCeOld = sloHiOld;
                      sloCeNew = sloHiNew;
                    }
                  else
                    {
                      //we have no slopes
                      sloCeOld = 0;
                      sloCeNew = 0;
                    }
                }
              else if (ivsHi && hasHi)
                {
                  //we only have slopes on the low side
                  sloCeOld = sloHiOld;
                  sloCeNew = sloHiNew;
                }
              else if (ivsLo && hasLo)
                {
                  //we only have slopes on the low side
                  sloCeOld = sloLoOld;
                  sloCeNew = sloLoNew;
                }
              else
                {
                  //we have no slopes
                  sloCeOld = 0;
                  sloCeNew = 0;
                }
            }

        }
    }
}
/************************************/
void
AggEBPWLFillPatch::
interpolateFAB(EBCellFAB& a_fine,
               const EBCellFAB& a_coarDataOld,
               const EBCellFAB& a_coarDataNew,
               const Real& a_coarTimeOld,
               const Real& a_coarTimeNew,
               const Real& a_fineTime,
               const DataIndex& a_dit,
               const Interval& a_variables) const
{
  CH_TIME("AggEBPWLFillPatch::interpolateFAB");
  getSlopes(a_coarDataOld, a_coarDataNew, a_dit, a_variables);
  //EBCellFAB has two data types, baseivfab has one
  CH_assert(a_coarDataOld.numDataTypes() == 2);
  CH_assert(a_coarDataNew.numDataTypes() == 2);
  CH_assert(a_fine.numDataTypes() == 2);
  //interpolation factor
  Real factor = 0.0;
  if ((a_coarTimeNew - a_coarTimeOld) > 1.0e-8)
    factor = (a_fineTime - a_coarTimeOld)/(a_coarTimeNew - a_coarTimeOld);

  for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
    {
      //BaseIVFAB has only one ptr
      const Real* slopPtrsDirOld[SpaceDim];
      const Real* slopPtrsDirNew[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          slopPtrsDirOld[idir] = m_slopeCeOld[idir][a_dit].dataPtr(icomp);
          slopPtrsDirNew[idir] = m_slopeCeNew[idir][a_dit].dataPtr(icomp);
        }
      const Real* coarPtrsOld[2];
      const Real* coarPtrsNew[2];
      Real*       finePtrsMid[2];
      for (int itype = 0; itype < 2; itype++)
        {
          coarPtrsOld[itype] = a_coarDataOld.dataPtr(itype, icomp);
          coarPtrsNew[itype] = a_coarDataNew.dataPtr(itype, icomp);
          finePtrsMid[itype] =        a_fine.dataPtr(itype, icomp);
        }
      for (int ifine = 0; ifine < m_fineOffsets[a_dit].size(); ifine++)
        {
          const int&    icoar   = m_fineOffsets[a_dit][ifine].slop_index;
          const size_t& fineOff = m_fineOffsets[a_dit][ifine].dest_access.offset;
          const int   & fineDat = m_fineOffsets[a_dit][ifine].dest_access.dataID;
          const size_t& coarOff = m_coarOffsets[a_dit][icoar].data_access.offset;
          const int   & coarDat = m_coarOffsets[a_dit][icoar].data_access.dataID;
          //only one data id for slopes
          const size_t& slopOff = m_coarOffsets[a_dit][icoar].slop_access.offset;

          const Real& coarOld = *(coarPtrsOld[coarDat] + coarOff);
          const Real& coarNew = *(coarPtrsNew[coarDat] + coarOff);
          Real&     interpVal = *(finePtrsMid[fineDat] + fineOff);

          //first set the fine stuff to the coar stuff
          Real fineValOld = coarOld;
          Real fineValNew = coarNew;

          //now add first derivative terms
          const IntVect& coariv = m_fineOffsets[a_dit][ifine].coariv;
          const IntVect& fineiv = m_fineOffsets[a_dit][ifine].fineiv;
          for (int derivDir = 0; derivDir  < SpaceDim; derivDir++)
            {
              const Real& slopOld = *(slopPtrsDirOld[derivDir] + slopOff);
              const Real& slopNew = *(slopPtrsDirNew[derivDir] + slopOff);
              Real coarLoc = Real(coariv[derivDir]) + 0.5;
              Real fineLoc  = (Real(fineiv[derivDir]) + 0.5)/Real(m_refRat);
              fineValOld += slopOld*(fineLoc - coarLoc);
              fineValNew += slopNew*(fineLoc - coarLoc);
            }

          //finally interpolate in time
         interpVal =  fineValOld + factor*(fineValNew-fineValOld);
        }
    }
}
/************************************/
void
AggEBPWLFillPatch::
definePieceWiseLinearFillPatch(const DisjointBoxLayout& a_dblFine,
                               const DisjointBoxLayout& a_dblCoar)
{
  if (m_patcher != NULL)
    {
      delete m_patcher;
    }

  m_patcher = new PiecewiseLinearFillPatch();
  m_patcher->define(a_dblFine,
                    a_dblCoar,
                    m_nComp,
                    m_coarDomain,
                    m_refRat,
                    m_radius);
}
/************************************/
#include "NamespaceFooter.H"
