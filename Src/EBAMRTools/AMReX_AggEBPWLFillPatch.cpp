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


#include "AMReX_AggEBPWLFillPatch.H"
#include "AMReX_EBCellFactory.H"

namespace amrex
{
  void null_deleter_sten(BaseStencil * a_sten)
 {}
  void null_deleter_ind(BaseIndex* a_sten)
 {}
  /************************************/
  void
  AggEBPWLFillPatch::
  setDefaultValues()
  {
    m_isDefined = false;
    m_refRat = -1;
    m_nComp = -1;
    m_radius = -1;
  }
  /************************************/
  AggEBPWLFillPatch::
  AggEBPWLFillPatch()
  {
    setDefaultValues();
  }
  /************************************/
  AggEBPWLFillPatch::
  ~AggEBPWLFillPatch()
  {
  }
  /************************************/
  void
  AggEBPWLFillPatch::
  define(const EBLevelGrid & a_eblgFine,
         const EBLevelGrid & a_eblgCoar,
         const int& a_nref,
         const int& a_nvar,
         const int& a_radius,
         const int& a_ghost,
         const bool& a_forceNoEBCF)
  {
    BL_PROFILE("AggEBPWLFillPatch::define");
    BL_ASSERT(a_nref > 0);                               //
    BL_ASSERT(a_nvar > 0);
    //first cell is interpolated, rest are extrapolated.
    //if this is 0, there is not even interpolation.
    BL_ASSERT(a_radius > 0);
    m_forceNoEBCF = a_forceNoEBCF;
    m_radius = a_radius;
    m_isDefined = true;
    m_refRat = a_nref;
    m_nComp = a_nvar;
    m_ghost = a_ghost;
    m_coarDomain = a_eblgCoar.getDomain();
     
    m_eblgFine = a_eblgFine;
    m_eblgCoar = a_eblgFine;
    //setup the non-EB PWLFillPatch is going to have to wait for fillpatch to gain generality
               
    //if(!m_forceNoEBCF)
    {
      if (m_eblgFine.getMaxCoarseningRatio() < m_refRat)
      {
        m_eblgFine.setMaxCoarseningRatio(m_refRat);
      }

      //needs to be extra because the stencil of
      //the derivs has an extra coarse cell
      const int coarse_slope_radius =
        (m_radius + m_refRat - 1) / m_refRat;
      const int coarse_ghost_radius = coarse_slope_radius + 1;
               
      m_coarGhostRad = coarse_ghost_radius;
      int ghost = m_coarGhostRad;
               
      coarsen(m_eblgCoFi, m_eblgFine, m_refRat);

      //  m_coarsenedFineEBISL.setMaxCoarseningRatio(m_refRat);
      EBCellFactory ebcellfact(m_eblgCoFi.getEBISL());
      m_coarOnFDataOld.define( m_eblgCoFi.getDBL(),
                               m_eblgCoFi.getDM(), m_nComp,
                               ghost, MFInfo(), ebcellfact);
      m_coarOnFDataNew.define( m_eblgCoFi.getDBL(),
                               m_eblgCoFi.getDM(), m_nComp,
                               ghost, MFInfo(), ebcellfact);
      makeStencils();
    }
  }
  /************************************/
  /************************************/
  void
  AggEBPWLFillPatch::
  makeStencils()
  {
    BL_PROFILE("AggEBPWLFillPatch::makeStencil");
//    BL_ASSERT(!m_forceNoEBCF);
    LayoutData< IntVectSet > irregRegionsFine;
    LayoutData< IntVectSet > irregRegionsCoFi;
    LayoutData< Vector<VolIndex> > srcVoFs;
    LayoutData< IntVectSet >  coarCeInterp[SpaceDim];
    LayoutData< IntVectSet >  coarLoInterp[SpaceDim];
    LayoutData< IntVectSet >  coarHiInterp[SpaceDim];
    LayoutData< Vector <VoFStencil> >  hiStencils[SpaceDim];
    LayoutData< Vector <VoFStencil> >  loStencils[SpaceDim];
    irregRegionsFine.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    irregRegionsCoFi.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    srcVoFs.define(         m_eblgFine.getDBL(), m_eblgFine.getDM());
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      coarCeInterp[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      coarLoInterp[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      coarHiInterp[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
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
    BL_PROFILE("AggEBPWLFillPatch::defineAggStencils");
//    BL_ASSERT(!m_forceNoEBCF);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_stenLo[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      m_stenHi[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
               
      for (MFIter mfi(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM()); mfi.isValid(); ++mfi)
      {
        Vector<std::shared_ptr<BaseIndex   > > baseindice(a_srcVoFs[mfi].size());
        Vector<std::shared_ptr<BaseStencil > > basestenlo(a_srcVoFs[mfi].size());
        Vector<std::shared_ptr<BaseStencil > > basestenhi(a_srcVoFs[mfi].size());
               
        for (int ivof= 0; ivof < a_srcVoFs[mfi].size(); ivof++)
        {
          baseindice[ivof] =   std::shared_ptr<BaseIndex>((BaseIndex*)(&a_srcVoFs[mfi][ivof]), &null_deleter_ind);
          basestenlo[ivof] = std::shared_ptr<BaseStencil>(    &a_loStencils[idir][mfi][ivof] , &null_deleter_sten);
          basestenhi[ivof] = std::shared_ptr<BaseStencil>(    &a_hiStencils[idir][mfi][ivof] , &null_deleter_sten);
        }
        //all the slopes are the same size so we can use any of them really
        m_stenLo[idir][mfi] = std::shared_ptr< AggStencil <EBCellFAB, BaseIVFAB<Real> > >
          (new AggStencil <EBCellFAB, BaseIVFAB<Real> >(baseindice, basestenlo, m_coarOnFDataOld[mfi], m_slopeLoOld[idir][mfi]));
        m_stenHi[idir][mfi] = std::shared_ptr< AggStencil <EBCellFAB, BaseIVFAB<Real> > >
          (new AggStencil <EBCellFAB, BaseIVFAB<Real> >(baseindice, basestenhi, m_coarOnFDataOld[mfi], m_slopeHiOld[idir][mfi]));
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
      m_slopeLoOld[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      m_slopeHiOld[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      m_slopeCeOld[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      m_slopeLoNew[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      m_slopeHiNew[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      m_slopeCeNew[idir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
               
      for (MFIter mfi(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM()); mfi.isValid(); ++mfi)
      {
        const EBGraph& ebgraph = m_eblgCoFi.getEBISL()[mfi].getEBGraph();
        const IntVectSet& ivs = a_irregRegionsCoFi[mfi];
        m_slopeLoOld[idir][mfi].define(ivs, ebgraph, m_nComp);
        m_slopeHiOld[idir][mfi].define(ivs, ebgraph, m_nComp);
        m_slopeCeOld[idir][mfi].define(ivs, ebgraph, m_nComp);
        m_slopeLoNew[idir][mfi].define(ivs, ebgraph, m_nComp);
        m_slopeHiNew[idir][mfi].define(ivs, ebgraph, m_nComp);
        m_slopeCeNew[idir][mfi].define(ivs, ebgraph, m_nComp);
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
    //first get the coarse offsets for the slopes and so on
    m_coarOffsets.define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());

    for (MFIter mfi(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM()); mfi.isValid(); ++mfi)
    {
      m_coarOffsets[mfi].resize(a_srcVoFsCoar[mfi].size());
      for (int ivof = 0; ivof < a_srcVoFsCoar[mfi].size(); ivof++)
      {
        //all the slope containers are the same size so I can
        //use any of them
        const VolIndex& vof = a_srcVoFsCoar[mfi][ivof];
        m_coarOffsets[mfi][ivof].slop_access.offset =     m_slopeLoOld[0][mfi].offset(  vof, 0);
        m_coarOffsets[mfi][ivof].slop_access.dataID =     m_slopeLoOld[0][mfi].dataType(vof);
        m_coarOffsets[mfi][ivof].data_access.offset =    m_coarOnFDataOld[mfi].offset(  vof, 0);
        m_coarOffsets[mfi][ivof].data_access.dataID =    m_coarOnFDataOld[mfi].dataType(vof);
               
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_coarOffsets[mfi][ivof].has_lo[idir] = (a_loStencils[idir][mfi][ivof].size() != 0);
          m_coarOffsets[mfi][ivof].has_hi[idir] = (a_hiStencils[idir][mfi][ivof].size() != 0);
               
          m_coarOffsets[mfi][ivof].ivs_lo[idir] = (a_coarLoInte[idir][mfi].contains(vof.gridIndex()));
          m_coarOffsets[mfi][ivof].ivs_hi[idir] = (a_coarHiInte[idir][mfi].contains(vof.gridIndex()));
          m_coarOffsets[mfi][ivof].ivs_ce[idir] = (a_coarCeInte[idir][mfi].contains(vof.gridIndex()));
        }
      }
    }
               
    //now get the offsets for the fine data so we can put the answers in
    //the proper place.  the dummy data holder is so we can get the
    //offsets
    m_fineOffsets.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    EBCellFactory ebcellfact(m_eblgFine.getEBISL());

    FabArray<EBCellFAB> dummyfine(m_eblgFine.getDBL(), 
                                  m_eblgFine.getDM(), 1, 
                                  m_ghost, MFInfo(), ebcellfact);

    for (MFIter mfi(dummyfine); mfi.isValid(); ++mfi)
    {
      const IntVectSet& ivs = a_irregRegionsFine[mfi];
      VoFIterator vofitfine(ivs, m_eblgFine.getEBISL()[mfi].getEBGraph());
      const Vector<VolIndex>& vofVec = vofitfine.getVector();
      m_fineOffsets[mfi].resize(vofVec.size());
      for (int ivof = 0; ivof < vofVec.size(); ivof++)
      {
        const VolIndex& vof = vofVec[ivof];
        m_fineOffsets[mfi][ivof].dest_access.offset = dummyfine[mfi].offset(  vof, 0);
        m_fineOffsets[mfi][ivof].dest_access.dataID = dummyfine[mfi].dataType(vof);
        //now find which coarse VoF this find vof coarsens to
        //(so we can use the slopes there)
        int coarind = -1;
        bool found = false;
        VolIndex coarVoF =  m_eblgFine.getEBISL().coarsen(vof, m_refRat, mfi);
        for (int icoar = 0; icoar < a_srcVoFsCoar[mfi].size(); icoar++)
        {
          const VolIndex& listVoF = a_srcVoFsCoar[mfi][icoar];
          if (listVoF == coarVoF)
          {
            found = true;
            coarind = icoar;
            break;
          }
        }
        if (!found)
        {
          amrex::Error("coarse-fine inconsistency in AggEBPWLFP");
        }
               
        m_fineOffsets[mfi][ivof].slop_index = coarind;
        m_fineOffsets[mfi][ivof].coariv = coarVoF.gridIndex();
        m_fineOffsets[mfi][ivof].fineiv =     vof.gridIndex();
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
    BL_PROFILE("AggEBPWLFillPatch::getIVS");
    Box domFine = refine(m_coarDomain, m_refRat);
//    BL_ASSERT(!m_forceNoEBCF);
    for (MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      IntVectSet&    localIrregF = a_irregRegionsFine[mfi];
      const EBISBox&   ebisBoxCF = m_eblgCoFi.getEBISL()[mfi];
      //const Box&        regionCF = ebisBoxCF.getRegion();

      //NOTE::this needs to be done for ALL boxes that are "near" the EB interface.
      //If one uses an isRegular check (as we did before) for this, you end up with inconsistent (yet accurate)
      //interpolations. Inconsistent (non-EB vs EB) interpolations cause problems for some codes that
      //rely on the fact that ghost cells are filled in exactly the same way.
      //This was a problem for the mac-projections because the advection algorithms' face velocities depend on
      //fillpatched ghost cells that must contain the same data where ghost cells overlap from one box to the next.
               
      //make localIrreg for all c/f ghost ivs within m_radius
      Box bigBox= m_eblgFine.getDBL()[mfi.index()];
      bigBox.grow(m_radius);
      bigBox &= domFine;
      localIrregF = IntVectSet(bigBox);
      for(int ibox = 0; ibox < m_eblgFine.getDBL().size(); ibox++)
      {
        localIrregF -= m_eblgFine.getDBL()[ibox];
      }

      //keep only fine ivs that are within 1 coarse radius of irregular coarse ivs
      //the remainder of fine ivs are done with the non-EB FillPatch
      ////turning this off because the non-eb stuff not ready
      //IntVectSet  irregCF = ebisBoxCF.getIrregIVS(regionCF);
      //irregCF.grow(1);
      //         
      //IntVectSet irregRefinedCF = irregCF;
      //irregRefinedCF.refine(m_refRat);
      //localIrregF &= irregRefinedCF;
               
      IntVectSet irregCoar = localIrregF;
      irregCoar.coarsen(m_refRat);
               
      a_irregRegionsCoar[mfi] = irregCoar;
      VoFIterator vofit(irregCoar, ebisBoxCF.getEBGraph());
      a_srcVoFs[mfi] = vofit.getVector();
               
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
    BL_PROFILE("AggEBPWLFillPatch::getLoHiCenIVS");
    //now create the coarse lo hi and center intvectsets--
    //these sets determine where you have to do one sided differences
    //because of the configuration of the coarse layout.
//    BL_ASSERT(!m_forceNoEBCF);

    for (MFIter mfi(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM()); mfi.isValid(); ++mfi) // 
    {
      const Box& fineBox = m_eblgFine.getDBL()[mfi];
      Box coarsenedFineBox = coarsen(grow(fineBox, m_radius), m_refRat);
      coarsenedFineBox &= m_coarDomain;
               
      IntVectSet coarsenedFineInterp(coarsenedFineBox);
               
      // Iterate over boxes in coarsened fine domain, and subtract off
      // from the set of coarse cells from which the fine ghost cells
      // will be interpolated.

      for(int ibox = 0; ibox < m_eblgCoFi.getDBL().size(); ibox++)
      {
        const Box& otherCoarsenedBox = m_eblgCoFi.getDBL()[ibox];
        coarsenedFineInterp -= otherCoarsenedBox;
      }
               
      // Now that we have the coarsened cells required for interpolation,
      // construct IntvectSets specifying the one-sided and centered
      // stencil locations.
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        IntVectSet& coarseCeInterp = a_coarCeInterp[idir][mfi];
        IntVectSet& coarseLoInterp = a_coarLoInterp[idir][mfi];
        IntVectSet& coarseHiInterp = a_coarHiInterp[idir][mfi];
               
        coarseCeInterp = coarsenedFineInterp;
               
        coarseLoInterp = coarseCeInterp;
        coarseLoInterp.shift(BASISV(idir));
               
        coarseHiInterp = coarseCeInterp;
        coarseHiInterp.shift(-BASISV(idir));
               
        // We iterate over the coarse grids and subtract them off of the
        // one-sided stencils.
        for(int ibox = 0; ibox < m_eblgCoar.getDBL().size(); ibox++)
        {
          const Box& bx = m_eblgCoar.getDBL()[ibox];
          coarseLoInterp -= bx;
          coarseHiInterp -= bx;
        }
               
        coarseLoInterp.shift(-BASISV(idir));
        coarseHiInterp.shift(BASISV(idir));
               
        coarseCeInterp -= coarseLoInterp;
        coarseCeInterp -= coarseHiInterp;
               
      }
    }
  }
               
  /************************************/
  /************************************/
  void
  AggEBPWLFillPatch::
  getSten(LayoutData<Vector<VoFStencil> >  a_loStencils  [SpaceDim],
          LayoutData<Vector<VoFStencil> >  a_hiStencils  [SpaceDim],
          LayoutData<IntVectSet>           a_coarLoInterp[SpaceDim],
          LayoutData<IntVectSet>           a_coarHiInterp[SpaceDim],
          LayoutData<IntVectSet>           a_coarCeInterp[SpaceDim],
          const LayoutData< Vector<VolIndex> >&   a_srcVoFs)
  {
    BL_PROFILE("AggEBPWLFillPatch::getSten");               //
//    BL_ASSERT(!m_forceNoEBCF);
               
    /////////////////////////
    for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      a_hiStencils[derivDir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
      a_loStencils[derivDir].define(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM());
               
      for (MFIter mfi(m_eblgCoFi.getDBL(), m_eblgCoFi.getDM()); mfi.isValid(); ++mfi) // 
      {
        const EBISBox& ebisBox  =  m_eblgCoFi.getEBISL()[mfi];
        a_loStencils[derivDir][mfi].resize(a_srcVoFs[mfi].size());
        a_hiStencils[derivDir][mfi].resize(a_srcVoFs[mfi].size());
               
        for (int ivof = 0; ivof < a_srcVoFs[mfi].size(); ivof++)
        {
          const VolIndex& vof = a_srcVoFs[mfi][ivof];
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
            a_loStencils[derivDir][mfi][ivof].clear();
            BL_ASSERT(stenVoFs.size() == stenWeig.size());
            for (int isten = 0; isten < stenVoFs.size(); isten++)
            {
              a_loStencils[derivDir][mfi][ivof].add(stenVoFs[isten], stenWeig[isten]);
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
            a_hiStencils[derivDir][mfi][ivof].clear();
            BL_ASSERT(stenVoFs.size() == stenWeig.size());
            for (int isten = 0; isten < stenVoFs.size(); isten++)
            {
              a_hiStencils[derivDir][mfi][ivof].add(stenVoFs[isten], stenWeig[isten]);
            }
          }
        } //end loop over vofs
      } //end loop over grids
    } //end loop over derivative directions
  }
  /************************************/
  void
  AggEBPWLFillPatch::
  interpolate(FabArray<EBCellFAB>& a_fineData,
              const FabArray<EBCellFAB>& a_coarDataOld,
              const FabArray<EBCellFAB>& a_coarDataNew,
              const Real& a_coarTimeOld,
              const Real& a_coarTimeNew,
              const Real& a_fineTime,
              int idst, int inco) const
  {
    BL_PROFILE("AggEBPWLFillPatch::interpolate");
    BL_ASSERT(isDefined());
    BL_ASSERT(a_coarTimeNew >= a_coarTimeOld);
    BL_ASSERT(a_coarTimeNew >= a_fineTime);
    BL_ASSERT(a_fineTime >= a_coarTimeOld);
               
    //do non-EB PWLFillPatch
    {
      //this bit will have to wait until fillpatch gets generalized
    }
               
     FabArray<EBCellFAB>& castOld= const_cast<FabArray<EBCellFAB>& >(a_coarDataOld); 
     FabArray<EBCellFAB>& castNew= const_cast<FabArray<EBCellFAB>& >(a_coarDataNew);
     castOld.FillBoundary();
     castNew.FillBoundary();

//    if(!m_forceNoEBCF)
    //now we have to do this everywhere
    {
      int srcGhost = 0;
      int dstGhost = m_coarGhostRad;

      m_coarOnFDataOld.copy(a_coarDataOld, idst, idst, inco, srcGhost, dstGhost);
      m_coarOnFDataNew.copy(a_coarDataNew, idst, idst, inco, srcGhost, dstGhost);
               

      for (MFIter mfi(m_coarOnFDataOld); mfi.isValid(); ++mfi) 
      {
        //Box fineBox = m_fineGrids[mfi];
        //Box coarsenedFineBox = m_coarsenedFineGrids[mfi];
        EBCellFAB&  fineData = a_fineData[mfi];
        const EBCellFAB&  coarDataOld = m_coarOnFDataOld[mfi];
        const EBCellFAB&  coarDataNew = m_coarOnFDataNew[mfi];
        // interpolateFAB interpolates from an entire coarse grid onto an
        // entire fine grids coarse-fine regions.
        interpolateFAB(fineData,
                       coarDataOld,
                       coarDataNew,
                       a_coarTimeOld,
                       a_coarTimeNew,
                       a_fineTime, mfi,
                       idst, inco);
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
      deltaminmod = rsign*std::min(std::abs(deltalo), std::abs(deltahi));
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
            const MFIter& a_dit,
            int idst, int inco) const
  {
    BL_PROFILE("AggEBPWLFillPatch::getSlopes");
//    BL_ASSERT(!m_forceNoEBCF);
    int ibeg = idst;  int isiz = inco;
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
        BL_ASSERT(m_slopeLoOld[idir][a_dit].numDataTypes() == 1);
        BL_ASSERT(m_slopeHiOld[idir][a_dit].numDataTypes() == 1);
        BL_ASSERT(m_slopeCeOld[idir][a_dit].numDataTypes() == 1);
        BL_ASSERT(m_slopeLoNew[idir][a_dit].numDataTypes() == 1);
        BL_ASSERT(m_slopeHiNew[idir][a_dit].numDataTypes() == 1);
        BL_ASSERT(m_slopeCeNew[idir][a_dit].numDataTypes() == 1);
               
               
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
                 const MFIter& a_dit,
                 int idst, int inco) const
  {
    BL_PROFILE("AggEBPWLFillPatch::interpolateFAB");
    getSlopes(a_coarDataOld, a_coarDataNew, a_dit, idst, inco);
    //EBCellFAB has two data types, baseivfab has one
    BL_ASSERT(a_coarDataOld.numDataTypes() == 2);
    BL_ASSERT(a_coarDataNew.numDataTypes() == 2);
    BL_ASSERT(a_fine.numDataTypes() == 2);
    //interpolation factor
    Real factor = 0.0;
    if ((a_coarTimeNew - a_coarTimeOld) > 1.0e-8)
      factor = (a_fineTime - a_coarTimeOld)/(a_coarTimeNew - a_coarTimeOld);
    int endcomp = idst+inco-1;
    for (int icomp = idst; icomp <= endcomp; icomp++)
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
}
               
