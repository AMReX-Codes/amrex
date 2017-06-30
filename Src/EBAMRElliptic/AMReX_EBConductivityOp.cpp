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

#include "NWOEBConductivityOp.H"
#include "AMReX_EBEllipticFort_F.H"
namespace amrex
{

  //-----------------------------------------------------------------------
  NWOEBConductivityOp::
  NWOEBConductivityOp(const EBLevelGrid &                          a_eblgFine,
                      const EBLevelGrid &                          a_eblg,
                      const EBLevelGrid &                          a_eblgCoar,
                      const EBLevelGrid &                          a_eblgCoarMG,
                      const shared_ptr<ConductivityBaseDomainBC>&  a_domainBC,
                      const shared_ptr<ConductivityBaseEBBC>&      a_ebBC,
                      const Real    &                              a_dx,
                      const Real    &                              a_dxCoar,
                      const int&                                   a_refToFine,
                      const int&                                   a_refToCoar,
                      const bool&                                  a_hasFine,
                      const bool&                                  a_hasCoar,
                      const bool&                                  a_hasCoarMG,
                      const Real&                                  a_alpha,
                      const Real&                                  a_beta,
                      const shared_ptr<FabArray<EBCellFAB> >&      a_acoef,
                      const shared_ptr<FabArray<EBFluxFAB> >&      a_bcoef,
                      const int    &                               a_ghostPhi,
                      const int    &                               a_ghostRHS)
  : LevelTGAHelmOp<FabArray<EBCellFAB>, EBFluxFAB>(false), // is time-independent
    m_ghostPhi(a_ghostPhi),
    m_ghostRHS(a_ghostRHS),
    m_eblg(a_eblg),
    m_eblgFine(),
    m_eblgCoar(),
    m_eblgCoarMG(),
    m_domainBC(a_domainBC),
    m_ebBC(a_ebBC),
    m_dxFine(),
    m_dx(a_dx),
    m_dxCoar(a_dxCoar),
    m_acoef(a_acoef),
    m_bcoef(a_bcoef),
    m_alpha(a_alpha),
    m_beta(a_beta),
    m_alphaDiagWeight(),
    m_betaDiagWeight(),
    m_refToFine(a_refToFine),
    m_refToCoar(a_refToCoar),
    m_hasFine(a_hasFine),
    m_hasCoar(a_hasCoar),
    m_hasCoarMG(a_hasCoarMG),
    m_ebAverage(),
    m_ebInterp(),
    m_opEBStencil(),
    m_relCoef(),
    m_vofIterIrreg(),
    m_vofIterMulti(),
    m_vofIterDomLo(),
    m_vofIterDomHi(),
    m_fastFR(),
    m_ebAverageMG(),
    m_ebInterpMG(),
    m_dblCoarMG(),
    m_ebislCoarMG(),
    m_domainCoarMG(),
    m_colors()
  {
    BL_PROFILE("NWOEBConductivityOp::ConductivityOp");
    int ncomp = 1;

    m_turnOffBCs = false; //REALLY needs to default to false

    if (m_hasFine)
    {
      m_eblgFine       = a_eblgFine;
      m_dxFine         = m_dx/a_refToFine;
    }

    EBCellFactory fact(m_eblg.getEBISL());
    m_relCoef.define(m_eblg.getDBL(), 1, IntVect::Zero, fact);
    if (m_hasCoar)
    {
      m_eblgCoar       = a_eblgCoar;
      m_ebInterp.define(  m_eblg, m_eblgCoar,  mgRef, m_ghostPhi);
      m_ebAverage.define( m_eblg, m_eblgCoar,  mgRef, m_ghostPhi);
    }

    if (m_hasCoarMG)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg, m_eblgCoarMG,  mgRef, m_ghostPhi);
      m_ebAverageMG.define(m_eblg, m_eblgCoarMG,  mgRef, m_ghostRHS);
    }

    //define stencils for the operator
    defineStencils();

  }
  //-----------------------------------------------------------------------
  Real
  NWOEBConductivityOp::
  getSafety()
  {
    Real safety = 1.0;
    return safety;
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  calculateAlphaWeight()
  {
    DataIterator dit = m_eblg.getDBL().dataIterator(); 
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      VoFIterator& vofit = m_vofIterIrreg[dit[mybox]];
      for (vofit.reset(); vofit.ok(); ++vofit)
      {
        const VolIndex& VoF = vofit();
        Real volFrac = m_eblg.getEBISL()[mfi].volFrac(VoF);
        Real alphaWeight = (*m_acoef)[mfi](VoF, 0);
        alphaWeight *= volFrac;

        m_alphaDiagWeight[mfi](VoF, 0) = alphaWeight;
      }
    }
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getDivFStencil(VoFStencil&      a_vofStencil,
                 const VolIndex&  a_vof,
                 const DataIndex& a_dit)
  {
    BL_PROFILE("NWOEBConductivityOp::getDivFStencil");
    const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit];
    a_vofStencil.clear();
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
        for (int iface = 0; iface < faces.size(); iface++)
        {
          VoFStencil fluxStencil;
          getFluxStencil(fluxStencil, faces[iface], a_dit);
          Real areaFrac = ebisBox.areaFrac(faces[iface]);
          fluxStencil *= Real(isign)*areaFrac/m_dx;
          a_vofStencil += fluxStencil;
        }
      }
    }
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getFluxStencil(VoFStencil&      a_fluxStencil,
                 const FaceIndex& a_face,
                 const MFIter   & a_mfi)
  {
    /// stencil for flux computation.   the truly ugly part of this computation
    /// beta and eta are multiplied in here

    BL_PROFILE("NWOEBConductivityOp::getFluxStencil");
    //need to do this by interpolating to centroids
    //so get the stencil at each face center and add with
    //interpolation weights
    FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                       (*m_eblg.getCFIVS())[a_dit],
                                                       m_eblg.getEBISL()[a_dit],
                                                       m_eblg.getDomain());

    a_fluxStencil.clear();
    for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, face, a_dit);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                             const FaceIndex& a_face,
                             const DataIndex& a_dit)
  {
    BL_PROFILE("NWOEBConductivityOp::getFaceCenteredFluxStencil");
    //face centered gradient is just a centered diff
    int faceDir= a_face.direction();
    a_fluxStencil.clear();

    if (!a_face.isBoundary())
    {
      a_fluxStencil.add(a_face.getVoF(Side::Hi),  1.0/m_dx, 0);
      a_fluxStencil.add(a_face.getVoF(Side::Lo), -1.0/m_dx, 0);
      a_fluxStencil *= (*m_bcoef)[a_dit][faceDir](a_face,0);
    }
    else
    {
      //the boundary condition handles this one.
    }
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  setAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
  {
    BL_PROFILE("NWOEBConductivityOp::setAlphaAndBeta");
    m_alpha = a_alpha;
    m_beta  = a_beta;
    calculateAlphaWeight(); //need to do this because the a coef has probably been changed under us
    calculateRelaxationCoefficient();
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  kappaScale(FabArray<EBCellFAB> & a_rhs)
  {
    BL_PROFILE("NWOEBConductivityOp::kappaScale");
    EBLevelDataOps::kappaWeight(a_rhs);
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  diagonalScale(FabArray<EBCellFAB> & a_rhs,
                bool a_kappaWeighted)
  {

    BL_PROFILE("NWOEBConductivityOp::diagonalScale");
    if (a_kappaWeighted)
      EBFabArrayOps::kappaWeight(a_rhs);

    //also have to weight by the identity coefficient
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mf)
    {
      a_rhs[mfi] *= (*m_acoef)[mfi];
    }
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  divideByIdentityCoef(FabArray<EBCellFAB> & a_rhs)
  {

    BL_PROFILE("NWOEBConductivityOp::divideByIdentityCoef");

    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mf)
    {
      a_rhs[mfi] /= (*m_acoef)[dit[mfi]];
    }
  }
  //-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  calculateRelaxationCoefficient()
  {
    BL_PROFILE("ebco::calculateRelCoef");
    // define regular relaxation coefficent
    Real safety = getSafety();
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mf)
    {
        const Box& grid = m_eblg.getDBL().get(dit[mybox]);
      
        // For time-independent acoef, initialize lambda = alpha * acoef.
        const EBCellFAB& acofab = (*m_acoef)[dit[mybox]];
        m_relCoef[dit[mybox]].setVal(0.);
        m_relCoef[dit[mybox]].plus(acofab, 0, 0, 1);
        m_relCoef[dit[mybox]]*= m_alpha;
      
        // Compute the relaxation coefficient in regular cells.
        BaseFab<Real>& regRel = m_relCoef[dit[mybox]].getSingleValuedFAB();
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          BaseFab<Real>& regBCo = (*m_bcoef)[dit[mybox]][idir].getSingleValuedFAB();
          ebefnd_decrinvrelcoefebco(BL_TO_FORTRAN_FAB(regRel),
                                    BL_TO_FORTRAN_FAB(regBCo,0),
                                    BL_TO_FORTRAN_BOX(grid),
                                    &m_beta, &m_dx, &idir);
        }
      
        //now invert so lambda = stable lambda for variable coef lapl
        //(according to phil, this is the correct lambda)
        FORT_INVERTLAMBDAEBCO(BL_TO_FORTRAN_FAB(regRel),
                              BL_TO_FORTRAN_BOX(grid), &safety);
      
        // Now go over the irregular cells.
        VoFIterator& vofit = m_vofIterIrreg[dit[mybox]];
        for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& VoF = vofit();
          
          Real alphaWeight = m_alphaDiagWeight[dit[mybox]](VoF, 0);
          Real  betaWeight =  m_betaDiagWeight[dit[mybox]](VoF, 0);
          alphaWeight *= m_alpha;
          betaWeight  *= m_beta;
          
          Real diagWeight  = alphaWeight + betaWeight;
          
          // Set the irregular relaxation coefficients.
          if (Abs(diagWeight) > 1.0e-9)
          {
            m_relCoef[dit[mybox]](VoF, 0) = safety/diagWeight;
          }
          else
          {
            m_relCoef[dit[mybox]](VoF, 0) = 0.;
          }
        }
    }
    }
  }
 
  void null_deleter_ebco_vof(VolIndex *)
 {}
  void null_deleter_ebco_sten(VolIndex *)
 {}

//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  defineStencils()
  {
    BL_PROFILE("NWOEBConductivityOp::defineStencils");
    // create ebstencil for irregular applyOp
    m_opEBStencil.define(m_eblg.getDBL(), m_eblg.getDM());
    // create vofstencils for applyOp and

    Real fakeBeta = 1;
    m_domainBC->setCoef(m_eblg,   fakeBeta ,      m_bcoef   );
    m_ebBC->setCoef(    m_eblg,   fakeBeta ,      m_bcoIrreg);

    Real dxScale = 1.0/m_dx;
    m_ebBC->define((*m_eblg.getCFIVS()), dxScale); //has to happen AFTER coefs are set
    LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);

    m_vofIterIrreg.define(     m_eblg.getDBL(), m_eblg.getDM()); // vofiterator cache
    m_vofIterMulti.define(     m_eblg.getDBL(), m_eblg.getDM()); // vofiterator cache
    m_alphaDiagWeight.define(  m_eblg.getDBL(), m_eblg.getDM());
    m_betaDiagWeight.define(   m_eblg.getDBL(), m_eblg.getDM());
    Box sideBoxLo[SpaceDim];
    Box sideBoxHi[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box domainBox = m_eblg.getDomain().domainBox();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofIterDomLo[idir].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofIterDomHi[idir].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }
    EBArith::getMultiColors(m_colors);

    EBCellFactory fact(m_eblg.getEBISL());              
    FabArray<EBCellFAB> phiProxy(m_eblg.getDBL(), , m_eblg.getDM(), SpaceDim, m_ghostPhi, MFInfo(), fact);
    FabArray<EBCellFAB> rhsProxy(m_eblg.getDBL(), , m_eblg.getDM(), SpaceDim, m_ghostRHS, MFInfo(), fact);

    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const Box& curBox = m_eblg.getDBL()[mfi];
      const EBISBox& ebisBox = m_eblg.getEBISL()[mfi];
      const EBGraph& ebgraph = ebisBox.getEBGraph();

      IntVectSet irregIVS = ebisBox.getIrregIVS(curBox);
      IntVectSet multiIVS = ebisBox.getMultiCells(curBox);

      //cache the vofIterators
      m_alphaDiagWeight[mfi].define(irregIVS,ebgraph, 1);
      m_betaDiagWeight [mfi].define(irregIVS,ebgraph, 1);
      m_vofIterIrreg   [mfi].define(irregIVS,ebgraph);
      m_vofIterMulti   [mfi].define(multiIVS,ebgraph);

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        IntVectSet loIrreg = irregIVS;
        IntVectSet hiIrreg = irregIVS;
        loIrreg &= sideBoxLo[idir];
        hiIrreg &= sideBoxHi[idir];
        m_vofIterDomLo[idir][mfi].define(loIrreg,ebgraph);
        m_vofIterDomHi[idir][mfi].define(hiIrreg,ebgraph);
      }

      VoFIterator &  vofit = m_vofIterIrreg[mfi];
      const vector<VolIndex>& vofvec = vofit.getVector();
      vector<VoFStencil> stenvec(vofvec.size());

      // cast from VolIndex to BaseIndex
      vector<shared_ptr<BaseIndex> >    dstVoF(vofvec.size());
      // cast from VoFStencil to BaseStencil
      vector<shared_ptr<BaseStencil> > stencil(vofvec.size());

      for(int ivec = 0; ivec < vofvec.size(); ivec++)
      {
        const VolIndex& vof = vofvec[ivec];
        VoFStencil& vofsten = stenvec[ivec];

        //bcoef is included here in the flux consistent
        //with the regular

        getDivFStencil(vofsten,vof, mfi);
        if (fluxStencil != NULL)
        {
          BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[mfi];
          //this fills the stencil with the gradient*beta*bcoef
          VoFStencil  fluxStencilPt = fluxStencilBaseIVFAB(vof,0);
          vofsten += fluxStencilPt;
        }
        //beta weight is the weight for the divergence term. 
        Real betaWeight = EBArith::getDiagWeight(vofsten, vof);
        const IntVect& iv = vof.gridIndex();
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box loSide = bdryLo(m_eblg.getDomain(),idir);
          loSide.shiftHalf(idir,1);
          Real adjust = 0;
          if (loSide.contains(iv))
          {
            Real faceAreaFrac = 0.0;
            Real weightedAreaFrac = 0;
            Vector<FaceIndex> faces = ebisBox.getFaces(vof,idir,Side::Lo);
            for (int i = 0; i < faces.size(); i++)
            {
              weightedAreaFrac = ebisBox.areaFrac(faces[i]);
              weightedAreaFrac *= (*m_bcoef)[mfi][idir](faces[i],0);
              faceAreaFrac +=  weightedAreaFrac;
            }
            adjust += -weightedAreaFrac /(m_dx*m_dx);
          }
          Box hiSide = bdryHi(m_eblg.getDomain(),idir);
          hiSide.shiftHalf(idir,-1);
          if (hiSide.contains(iv))
          {
            Real faceAreaFrac = 0.0;
            Real weightedAreaFrac = 0;
            Vector<FaceIndex> faces = ebisBox.getFaces(vof,idir,Side::Hi);
            for (int i = 0; i < faces.size(); i++)
            {
              weightedAreaFrac = ebisBox.areaFrac(faces[i]);
              weightedAreaFrac *= (*m_bcoef)[mfi][idir](faces[i],0);
              faceAreaFrac +=  weightedAreaFrac;
            }
            adjust += -weightedAreaFrac /(m_dx*m_dx);
          }
          betaWeight += adjust;
        }
            
        dstVoF[ivec]  = shared_ptr<BaseIndex  >(&vof    , &null_deleter_ebco_vof);
        stencil[ivec] = shared_ptr<BaseStencil>(&vofsten, &null_deleter_ebco_sten);

        //weight for in identity term
        Real volFrac = ebisBox.volFrac(vof);
        Real alphaWeight = (*m_acoef)[mfi](vof, 0);
        alphaWeight *= volFrac;

        m_alphaDiagWeight[mfi](vof, 0) = alphaWeight;
        m_betaDiagWeight[ mfi](vof, 0)  = betaWeight;
      }

      m_opEBStencil[mfi] = shared_ptr<VCAggStencil>
        (new VCAggStencil(dstVoF, stencil, phiProxy[mfi], rhsProxy[mfi], m_relCoef[mfi], m_alphaDiagWeight[mfi],1));

    }//mfi

    calculateAlphaWeight();
    calculateRelaxationCoefficient();

    if (m_hasFine)
    {
      int ncomp = 1;
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp, s_forceNoEBCF);
    }

  }

//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  residual(FabArray<EBCellFAB>&       a_residual,
           const FabArray<EBCellFAB>& a_phi,
           const FabArray<EBCellFAB>& a_rhs,
           bool                        a_homogeneousPhysBC)
  {
    BL_PROFILE("NWOEBConductivityOp::residual");
    //this is a multigrid operator so only homogeneous CF BC
    //and null coar level
    BL_ASSERT(a_residual.nGrow() == m_ghostRHS);
    BL_ASSERT(a_phi.nGrow() == m_ghostPhi);
    applyOp(a_residual,a_phi,NULL, a_homogeneousPhysBC, true);
    axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  preCond(FabArray<EBCellFAB>&       a_lhs,
          const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("NWOEBConductivityOp::preCond");
    EBFabArrayOps::assign(a_lhs, a_rhs);
    EBFabArrayOps::scale(a_lhs, m_relCoef);

    relax(a_lhs, a_rhs, 40);
  }

//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  applyOp(FabArray<EBCellFAB>&                    a_lhs,
          const FabArray<EBCellFAB>&              a_phi,
          const FabArray<EBCellFAB>* const        a_phiCoar,
          const bool&                              a_homogeneousPhysBC,
          const bool&                              a_homogeneousCFBC)
  {
    FabArray<EBCellFAB>& phi = (FabArray<EBCellFAB>&)a_phi;
    if (m_hasCoar)
    {
      if(a_homogeneousCFBC)
      {
        m_interpWithCoarser->coarseFineInterpH(phi,  0, 0, 1);
      }
      else
      {
        m_interpWithCoarser->coarseFineInterp(phi, *a_phiCoar,  0, 0, 1);
        //          dumpEBLevelGhost(&phi);
      }
    }
    applyOp(a_lhs, a_phi, a_homogeneousPhysBC);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  applyOp(FabArray<EBCellFAB>&                    a_lhs,
          const FabArray<EBCellFAB>&              a_phi,
          bool                                     a_homogeneous)
  {
    BL_PROFILE("ebco::applyOp");
    FabArray<EBCellFAB>&  phi = (FabArray<EBCellFAB>&) a_phi;
    {
      BL_PROFILE("ghostcell fills");
      fillPhiGhost(a_phi, a_homogeneous);
      phi.exchange(m_exchangeCopier);
    }

    /**
       This is not done anymore because alpha and a are now part of the stencil
       and part of regular apply
       EBFabArrayOps::setToZero(a_lhs);
       incr( a_lhs, a_phi, m_alpha);
    **/
    {
      BL_PROFILE("applying op without bcs");
      DataIterator dit = m_eblg.getDBL().dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {

        applyOpRegular(  a_lhs[dit[mybox]], a_phi[dit[mybox]], a_homogeneous, dit[mybox]);
        applyOpIrregular(a_lhs[dit[mybox]], a_phi[dit[mybox]], a_homogeneous, dit[mybox]);
      }
    }
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  fillPhiGhost(const FabArray<EBCellFAB>& a_phi, bool a_homog) const
  {
    BL_PROFILE("nwoebco::fillghostLD");
    DataIterator dit = m_eblg.getDBL().dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
      fillPhiGhost(a_phi[dit[mybox]], dit[mybox], a_homog);
    }
    FabArray<EBCellFAB>& phi = (FabArray<EBCellFAB>&)(a_phi);
    phi.exchange(m_exchangeCopier);
  }
  void
  NWOEBConductivityOp::
  fillPhiGhost(const EBCellFAB& a_phi, const DataIndex& a_datInd, bool a_homog) const
  {
    BL_PROFILE("nwoebco::fillghostfab");

    EBCellFAB& phi = (EBCellFAB&) a_phi;
    ConductivityBaseDomainBC* condbc = dynamic_cast<ConductivityBaseDomainBC*>(&(*m_domainBC));
    Box grid = m_eblg.getDBL()[a_datInd];
    Box domBox = m_eblg.getDomain().domainBox();
    if(condbc == NULL)
    {
      MayDay::Error("dynamic cast failed");
    }
    if (!s_turnOffBCs)
    {
      FArrayBox& fab = phi.getFArrayBox();
      condbc->fillPhiGhost(fab, grid, domBox, m_dx, a_homog);
    }
    else
    {
      Box valid = m_eblg.getDBL()[a_datInd];
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          FArrayBox& fab = phi.getFArrayBox();
          ExtrapolateBC(fab, valid,  m_dx, idir, sit());
        }
        //grow so that we hit corners
        valid.grow(idir, 1);
      }
    }

  }
/*****/
  void
  NWOEBConductivityOp::
  applyOpRegular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const DataIndex&       a_datInd)
  {
    BL_PROFILE("nwoebco::applyopRegular");
    //assumes ghost cells already filled

    const Box& grid = m_eblg.getDBL()[a_datInd];
    BaseFab<Real>      & lphfab  =  a_lhs.getSingleValuedFAB();
    const BaseFab<Real>& phifab  =  a_phi.getSingleValuedFAB();
    const BaseFab<Real>& acofab  =(*m_acoef)[a_datInd].getSingleValuedFAB();
    const EBFluxFAB& bco  =   (*m_bcoef)[a_datInd];

    Vector<const BaseFab<Real>* > bcoside(3, &(bco[0].getSingleValuedFAB()));
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      bcoside[idir] = &(bco[idir].getSingleValuedFAB());
    }

    FORT_APPLYOPEBCONDNOBCS(CHF_FRA1(lphfab,0),
                            CHF_CONST_FRA1(phifab,0),
                            CHF_CONST_FRA1(acofab,0),
                            CHF_CONST_FRA1((*bcoside[0]), 0),
                            CHF_CONST_FRA1((*bcoside[1]), 0),
                            CHF_CONST_FRA1((*bcoside[2]), 0),
                            CHF_CONST_REAL(m_dx),
                            CHF_CONST_REAL(m_alpha),
                            CHF_CONST_REAL(m_beta),
                            CHF_BOX(grid));

  }
//-----------------------------------------------------------------------

  void
  NWOEBConductivityOp::
  applyOpIrregular(EBCellFAB&             a_lhs,
                   const EBCellFAB&       a_phi,
                   const bool&            a_homogeneous,
                   const DataIndex&       a_datInd)
  {
    BL_PROFILE("nwoebco::applyOpIrregular");
    RealVect vectDx = m_dx*RealVect::Unit;

    m_opEBStencil[a_datInd]->apply(a_lhs, a_phi, m_alphaDiagWeight[a_datInd], m_alpha, m_beta, 0, false);
  
    if (!a_homogeneous)
    {
      const Real factor = m_beta/m_dx;
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIterIrreg[a_datInd], (*m_eblg.getCFIVS()),
                          a_datInd, RealVect::Zero, vectDx, factor,
                          a_homogeneous, 0.0);
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      int comp = 0;
      for (m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
      {
        Real flux;
        const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
        m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                RealVect::Zero,vectDx,idir,Side::Lo, a_datInd, 0.0,
                                a_homogeneous);

        //area gets multiplied in by bc operator
        a_lhs(vof,comp) -= flux*m_beta/m_dx;
      }
      for (m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
      {
        Real flux;
        const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
        m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                RealVect::Zero,vectDx,idir,Side::Hi,a_datInd,0.0,
                                a_homogeneous);

        //area gets multiplied in by bc operator
        a_lhs(vof,comp) += flux*m_beta/m_dx;
      }
    }

  }

//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  create(FabArray<EBCellFAB>&       a_lhs,
         const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("ebco::create");
    int ncomp = a_rhs.nComp();
    EBCellFactory ebcellfact(m_eblg.getEBISL());
    a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.nGrow(), ebcellfact);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  createCoarsened(FabArray<EBCellFAB>&       a_lhs,
                  const FabArray<EBCellFAB>& a_rhs,
                  const int &                 a_refRat)
  {
    BL_PROFILE("ebco::createCoar");
    int ncomp = a_rhs.nComp();
    IntVect ghostVect = a_rhs.nGrow();

    BL_ASSERT(m_eblg.getDBL().coarsenable(a_refRat));

    //fill ebislayout
    DisjointBoxLayout dblCoarsenedFine;
    coarsen(dblCoarsenedFine, m_eblg.getDBL(), a_refRat);

    EBISLayout ebislCoarsenedFine;
    IntVect ghostVec = a_rhs.nGrow();

    Box coarDom = coarsen(m_eblg.getDomain(), a_refRat);
    m_eblg.getEBIS()->fillEBISLayout(ebislCoarsenedFine, dblCoarsenedFine, coarDom , ghostVec[0]);
    if (m_refToCoar > 1)
    {
      ebislCoarsenedFine.setMaxRefinementRatio(m_refToCoar, m_eblg.getEBIS());
    }

    //create coarsened data
    EBCellFactory ebcellfactCoarsenedFine(ebislCoarsenedFine);
    a_lhs.define(dblCoarsenedFine, ncomp,ghostVec, ebcellfactCoarsenedFine);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  assign(FabArray<EBCellFAB>&       a_lhs,
         const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("ebco::assign");
    EBFabArrayOps::assign(a_lhs,a_rhs);
  }
//-----------------------------------------------------------------------
  Real
  NWOEBConductivityOp::
  dotProduct(const FabArray<EBCellFAB>& a_1,
             const FabArray<EBCellFAB>& a_2)
  {
    BL_PROFILE("ebco::dotProd");
    Box domain;
    Real volume;

    return EBFabArrayOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  incr(FabArray<EBCellFAB>&       a_lhs,
       const FabArray<EBCellFAB>& a_x,
       Real                        a_scale)
  {
    BL_PROFILE("ebco::incr");
    EBFabArrayOps::incr(a_lhs,a_x,a_scale);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  axby(FabArray<EBCellFAB>&       a_lhs,
       const FabArray<EBCellFAB>& a_x,
       const FabArray<EBCellFAB>& a_y,
       Real                        a_a,
       Real                        a_b)
  {
    BL_PROFILE("ebco::axby");
    EBFabArrayOps::axby(a_lhs,a_x,a_y,a_a,a_b);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  scale(FabArray<EBCellFAB>& a_lhs,
        const Real&           a_scale)
  {
    BL_PROFILE("ebco::scale");
    EBFabArrayOps::scale(a_lhs,a_scale);
  }
//-------------------------------
  Real 
  NWOEBConductivityOp::
  norm(const FabArray<EBCellFAB>& a_rhs,
       int                         a_ord)
  {
    BL_PROFILERS("NWOEBConductivityOp::norm");
    BL_PROFILER("mpi_allreduce",t1);

    Real maxNorm = 0.0;

    maxNorm = localMaxNorm(a_rhs);

    CH_START(t1);
#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL,
                               MPI_MAX, Chombo_MPI::comm);
    if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communcation error on norm");
    }
    maxNorm = tmp;
#endif
    //  Real volume=1.;
    //  EBFabArrayOps::gatherBroadCast(maxNorm, volume, 0);
    CH_STOP(t1);

    return maxNorm;
  }

  Real 
  NWOEBConductivityOp::
  localMaxNorm(const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("EBAMRPoissonOp::localMaxNorm");
    return  EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);
    //  Box domain;
    //  Real volume;
    //  return EBFabArrayOps::kappaNorm(volume,a_rhs,EBLEVELDATAOPS_ALLVOFS,domain,0);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  setToZero(FabArray<EBCellFAB>& a_lhs)
  {
    BL_PROFILE("ebco::setToZero");
    EBFabArrayOps::setToZero(a_lhs);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  setVal(FabArray<EBCellFAB>& a_lhs, const Real& a_value)
  {
    BL_PROFILE("ebco::setVal");
    EBFabArrayOps::setVal(a_lhs, a_value);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  createCoarser(FabArray<EBCellFAB>&       a_coar,
                const FabArray<EBCellFAB>& a_fine,
                bool                        a_ghosted)
  {
    BL_PROFILE("ebco::createCoarser");
    BL_ASSERT(a_fine.nComp() == 1);
    const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
    Box coarDom = coarsen(m_eblg.getDomain(), 2);

    int nghost = a_fine.nGrow()[0];
    EBISLayout coarEBISL;

    const EBIndexSpace* const ebisPtr = m_eblg.getEBIS();
    ebisPtr->fillEBISLayout(coarEBISL,
                            dbl, coarDom, nghost);

    EBCellFactory ebcellfact(coarEBISL);
    a_coar.define(dbl, 1,a_fine.nGrow(),ebcellfact);
  }
//------------------------------------------------------------------------
/*****/
  void
  NWOEBConductivityOp::
  relax(FabArray<EBCellFAB>&       a_phi,
        const FabArray<EBCellFAB>& a_rhs,
        int                         a_iterations)
  {
    BL_PROFILE("nwoebco::relax");
    BL_ASSERT(a_phi.isDefined());
    BL_ASSERT(a_rhs.isDefined());
    BL_ASSERT(a_phi.nGrow() >= IntVect::Unit);
    BL_ASSERT(a_phi.nComp() == a_rhs.nComp());

    //change to false to get a lot of printouts
    static bool printedStuff = true;

    if(!printedStuff)
    {
      pout() << "rhs: " << endl;;
      FabArray<EBCellFAB>* rhs = (FabArray<EBCellFAB>*)(&a_rhs);
      printMaxMinLDCell(rhs);
    }
    // do first red, then black passes
    for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
      {
        //after this lphi = L(phi)
        //this call contains bcs and exchange
        if ((icolor == 0) || (!s_doLazyRelax))
        {
          BL_PROFILE("ghostfill and homogcfinterp");
          fillPhiGhost(a_phi, true);
          homogeneousCFInterp(a_phi);
        }
        if(!printedStuff)
        {
          pout() << "iter = "<< whichIter << ", icolor = " << icolor << endl;
          pout() << "phi: " ;
          printMaxMinLDCell(&a_phi);
        }
        gsrbColor(a_phi, a_rhs, m_colors[icolor]);
      }
    }

    //MayDay::Abort("leaving after first relax");
    if(!printedStuff)
    {
      pout() << "rhs: " << endl;;
      FabArray<EBCellFAB>* rhs = (FabArray<EBCellFAB>*)(&a_rhs);
      printMaxMinLDCell(rhs);
      pout() << "phi coming out of relax: "  << endl;
      printMaxMinLDCell(&a_phi);
    }
    printedStuff = true;
  }

  void
  NWOEBConductivityOp::
  homogeneousCFInterp(FabArray<EBCellFAB>&   a_phif)
  {
    BL_PROFILE("nwoebco::homog_cfinterp");
    if (m_hasCoar)
    {
      m_interpWithCoarser->coarseFineInterpH(a_phif,  0, 0, 1);
    }
  }

//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  gsrbColor(FabArray<EBCellFAB>&       a_phi,
            const FabArray<EBCellFAB>& a_rhs,
            const IntVect&              a_color)
  {

    BL_PROFILE("nwoebco::gsrbColor");
    const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
    DataIterator dit = m_eblg.getDBL().dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
      //first do the the regular stuff
      Box dblBox  = dbl.get(dit[mybox]);
      const DataIndex datInd = dit[mybox];
      BaseFab<Real>&       regPhi =     a_phi[datInd].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[datInd].getSingleValuedFAB();
      const BaseFab<Real>& regRel = m_relCoef[datInd].getSingleValuedFAB();


      //assumes ghost cells already filled

      const BaseFab<Real>& acofab  =(*m_acoef)[datInd].getSingleValuedFAB();
      const EBFluxFAB& bco  =       (*m_bcoef)[datInd];

      Vector<const BaseFab<Real>* > bcoside(3, &(bco[0].getSingleValuedFAB()));
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        bcoside[idir] = &(bco[idir].getSingleValuedFAB());
      }

      IntVect loIV = dblBox.smallEnd();
      IntVect hiIV = dblBox.bigEnd();
        
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (loIV[idir] % 2 != a_color[idir])
        {
          loIV[idir]++;
        }
      }
        
      m_opEBStencil[datInd]->cachePhi(a_phi[datInd]);

      if (loIV <= hiIV)
      {
        Box coloredBox(loIV, hiIV);
        FORT_GSCOLOREBCONDOP(CHF_FRA1(regPhi,0),
                             CHF_CONST_FRA1(regRhs,0),
                             CHF_CONST_FRA1(regRel,0),
                             CHF_CONST_FRA1(acofab,0),
                             CHF_CONST_FRA1((*bcoside[0]), 0),
                             CHF_CONST_FRA1((*bcoside[1]), 0),
                             CHF_CONST_FRA1((*bcoside[2]), 0),
                             CHF_CONST_REAL(m_dx),
                             CHF_CONST_REAL(m_alpha),
                             CHF_CONST_REAL(m_beta),
                             CHF_BOX(coloredBox));
      }

      m_opEBStencil[datInd]->uncachePhi(a_phi[datInd]);


      m_opEBStencil[datInd]->relax(a_phi[datInd], a_rhs[datInd], 
                                   m_relCoef[datInd], m_alphaDiagWeight[datInd],
                                   m_alpha, m_beta, 0, a_color);
    }
  }

//-----------------------------------------------------------------------
  void NWOEBConductivityOp::
  restrictResidual(FabArray<EBCellFAB>&       a_resCoar,
                   FabArray<EBCellFAB>&       a_phiThisLevel,
                   const FabArray<EBCellFAB>& a_rhsThisLevel)
  {
    BL_PROFILE("NWOEBConductivityOp::restrictResidual");

    BL_ASSERT(a_resCoar.nComp() == 1);
    BL_ASSERT(a_phiThisLevel.nComp() == 1);
    BL_ASSERT(a_rhsThisLevel.nComp() == 1);

    FabArray<EBCellFAB> resThisLevel;
    bool homogeneous = true;

    EBCellFactory ebcellfactTL(m_eblg.getEBISL());
    IntVect ghostVec = a_rhsThisLevel.nGrow();

    resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

    // Get the residual on the fine grid
    residual(resThisLevel,a_phiThisLevel,a_rhsThisLevel,homogeneous);

    // now use our nifty averaging operator
    Interval variables(0, 0);
    if (m_layoutChanged)
    {
      m_ebAverageMG.average(a_resCoar, resThisLevel, variables);
    }
    else
    {
      m_ebAverageMG.averageMG(a_resCoar, resThisLevel, variables);
    }
  }
//-----------------------------------------------------------------------
  void NWOEBConductivityOp::
  prolongIncrement(FabArray<EBCellFAB>&       a_phiThisLevel,
                   const FabArray<EBCellFAB>& a_correctCoar)
  {
    BL_PROFILE("NWOEBConductivityOp::prolongIncrement");
    Interval vars(0, 0);
    if (m_layoutChanged)
    {
      m_ebInterpMG.pwcInterp(a_phiThisLevel, a_correctCoar, vars);
    }
    else
    {
      m_ebInterpMG.pwcInterpMG(a_phiThisLevel, a_correctCoar, vars);
    }
  }
//-----------------------------------------------------------------------
  int NWOEBConductivityOp::
  refToCoarser()
  {
    return m_refToCoar;
  }
//-----------------------------------------------------------------------
  int NWOEBConductivityOp::
  refToFiner()
  {
    return m_refToFine;
  }
//-----------------------------------------------------------------------
  void NWOEBConductivityOp::
  AMRResidual(FabArray<EBCellFAB>&       a_residual,
              const FabArray<EBCellFAB>& a_phiFine,
              const FabArray<EBCellFAB>& a_phi,
              const FabArray<EBCellFAB>& a_phiCoar,
              const FabArray<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC,
              AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILERS("NWOEBConductivityOp::AMRResidual");
    BL_PROFILER("AMROperator", t1);
    BL_PROFILER("axby", t2);
    BL_ASSERT(a_residual.nGrow() == m_ghostRHS);
    BL_ASSERT(a_rhs.nGrow() == m_ghostRHS);
    BL_ASSERT(a_residual.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);
    BL_ASSERT(a_rhs.nComp() == 1);

    CH_START(t1);
    AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar,
                a_homogeneousPhysBC, a_finerOp);
    CH_STOP(t1);

    //multiply by -1 so a_residual now holds -L(phi)
    //add in rhs so a_residual = rhs - L(phi)
    CH_START(t2);
    axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
    CH_STOP(t2);
  }
//-----------------------------------------------------------------------
  void NWOEBConductivityOp::
  AMROperator(FabArray<EBCellFAB>&       a_LofPhi,
              const FabArray<EBCellFAB>& a_phiFine,
              const FabArray<EBCellFAB>& a_phi,
              const FabArray<EBCellFAB>& a_phiCoar,
              bool a_homogeneousPhysBC,
              AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILERS("NWOEBConductivityOp::AMROperator");
    BL_PROFILER("applyOp", t1);
    BL_PROFILER("reflux", t2);
    BL_ASSERT(a_LofPhi.nGrow() == m_ghostRHS);
    BL_ASSERT(a_LofPhi.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);

    //apply the operator between this and the next coarser level.
    CH_START(t1);
    applyOp(a_LofPhi, a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);
    CH_STOP(t1);

    //now reflux to enforce flux-matching from finer levels
    if (m_hasFine)
    {
      BL_ASSERT(a_finerOp != NULL);
      CH_START(t2);

      reflux(a_LofPhi, a_phiFine, a_phi, a_finerOp);

      CH_STOP(t2);
    }
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  reflux(FabArray<EBCellFAB>& a_residual,
         const FabArray<EBCellFAB>& a_phiFine,
         const FabArray<EBCellFAB>& a_phi,
         AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILERS("NWOEBConductivityOp::fastReflux");
    BL_PROFILER("setToZero",t2);
    BL_PROFILER("incrementCoar",t3);
    BL_PROFILER("incrementFine",t4);
    BL_PROFILER("reflux_from_reg",t5);
    Interval interv(0,0);

    CH_START(t2);
    m_fastFR.setToZero();
    CH_STOP(t2);
    CH_START(t3);
    incrementFRCoar(m_fastFR, a_phiFine, a_phi);
    CH_STOP(t3);

    CH_START(t4);
    incrementFRFine(m_fastFR, a_phiFine, a_phi, a_finerOp);
    CH_STOP(t4);
    CH_START(t5);

    Real scale = 1.0/m_dx;
    m_fastFR.reflux(a_residual, interv, scale);

    CH_STOP(t5);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  incrementFRCoar(EBFastFR&             a_fluxReg,
                  const FabArray<EBCellFAB>& a_phiFine,
                  const FabArray<EBCellFAB>& a_phi)
  {
    BL_PROFILE("NWOEBConductivityOp::incrementFRCoar");
    BL_ASSERT(a_phiFine.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);

    int ncomp = 1;
    Interval interv(0,0);

    DataIterator dit = m_eblg.getDBL().dataIterator(); 
    int nbox = dit.size();
#pragma omp parallel for
    for(int mybox=0; mybox<nbox; mybox++)
    {

      const EBCellFAB& coarfab = a_phi[dit[mybox]];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit[mybox]];
      const Box&  box = m_eblg.getDBL().get(dit[mybox]);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        //no boundary faces here.

        Box ghostedBox = box;
        ghostedBox.grow(1);
        ghostedBox.grow(idir,-1);
        ghostedBox &= m_eblg.getDomain();

        EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);

        //old way
        //getFlux(coarflux, coarfab, ghostedBox, box, m_eblg.getDomain(), ebisBox, m_dx, dit[mybox], idir);

        // new way
        fillPhiGhost(a_phi[dit[mybox]], dit[mybox], false);
        getFluxRegOnly(coarflux, coarfab, ghostedBox, m_dx, dit[mybox], idir);

        Real scale = 1.0; //beta and bcoef already in flux
        for (SideIterator sit; sit.ok(); ++sit)
        {
          a_fluxReg.incrementCoarseBoth(coarflux, scale, dit[mybox], interv, idir, sit());
        }
      }
    }
  }

//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getFlux(EBFluxFAB&                    a_flux,
          const FabArray<EBCellFAB>&   a_data,
          const Box&                    a_grid,
          const DataIndex&              a_dit,
          Real                          a_scale)
  {
    BL_PROFILE("ebco::getflux1");
    a_flux.define(m_eblg.getEBISL()[a_dit], a_grid, 1);
    a_flux.setVal(0.);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box ghostedBox = a_grid;
      ghostedBox.grow(1);
      ghostedBox.grow(idir,-1);
      ghostedBox &= m_eblg.getDomain();
      fillPhiGhost(a_data[a_dit], a_dit, false);

      getFlux(a_flux[idir], a_data[a_dit], ghostedBox, a_grid,
              m_eblg.getDomain(), m_eblg.getEBISL()[a_dit], m_dx, a_dit, idir);

    }
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getFlux(EBFaceFAB&                    a_fluxCentroid,
          const EBCellFAB&              a_phi,
          const Box&                    a_ghostedBox,
          const Box&                    a_fabBox,
          const Box&          a_domain,
          const EBISBox&                a_ebisBox,
          const Real&                   a_dx,
          const DataIndex&              a_datInd,
          const int&                    a_idir)
  {
    BL_PROFILE("ebco::getFlux2");
    //has some extra cells so...
    a_fluxCentroid.setVal(0.);
    int ncomp = a_phi.nComp();
    BL_ASSERT(ncomp == a_fluxCentroid.nComp());
    Box cellBox = a_ghostedBox;
    //want only interior faces
    cellBox.grow(a_idir, 1);
    cellBox &= a_domain;
    cellBox.grow(a_idir,-1);

    Box faceBox = surroundingNodes(cellBox, a_idir);
    EBFaceFAB fluxCenter(a_ebisBox, a_ghostedBox, a_idir,1);

    //make a EBFaceFAB (including ghost cells) that will hold centered gradients
    BaseFab<Real>& regFlux  = fluxCenter.getSingleValuedFAB();
    const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
    const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
    const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());

    FORT_GETFLUXEBCO(CHF_FRA1(regFlux,0),
                     CHF_CONST_FRA1(regBCo,0),
                     CHF_CONST_FRA1(regPhi, 0),
                     CHF_BOX(faceBox),
                     CHF_CONST_REAL(a_dx),
                     CHF_CONST_INT(a_idir));


    a_fluxCentroid.copy(fluxCenter);

    // VolIndex vofdeblo1(IntVect(D_DECL(223, 247,0)), 0);
    // VolIndex vofdebhi1(IntVect(D_DECL(223, 248,0)), 0);
    // VolIndex vofdeblo2(IntVect(D_DECL(224, 247,0)), 0);
    // VolIndex vofdebhi2(IntVect(D_DECL(224, 248,0)), 0);
    // FaceIndex facedeb1(vofdeblo1, vofdebhi1, 1);
    // FaceIndex facedeb2(vofdeblo2, vofdebhi2, 1);

    IntVectSet ivsCell = a_ebisBox.getIrregIVS(cellBox);
    if (!ivsCell.isEmpty())
    {
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;

      for (FaceIterator faceit(ivsCell, a_ebisBox.getEBGraph(), a_idir,stopCrit);
           faceit.ok(); ++faceit)
      {
        const FaceIndex& face = faceit();
        Real phiHi = a_phi(face.getVoF(Side::Hi), 0);
        Real phiLo = a_phi(face.getVoF(Side::Lo), 0);
        Real fluxFace = bcoebff(face, 0)*(phiHi - phiLo)/a_dx;
        //          if (EBCellFAB::s_verbose && ((face==facedeb1) || (face==facedeb2)))
        //            {
        //              pout() << "NWOEBConductivityOp::getFlux at "<< face ;
        //              pout() << ", phiHi, phiLo, flux = " << phiHi << ", " << phiLo << ", "<< fluxFace << endl;
        //            }
        fluxCenter(face, 0) = fluxFace;
      }
      //interpolate from face centers to face centroids
      Box cellBox = a_fluxCentroid.getCellRegion();
      EBArith::interpolateFluxToCentroids(a_fluxCentroid,
                                          fluxCenter,
                                          a_fabBox,
                                          a_ebisBox,
                                          a_domain,
                                          a_idir);
    }

    a_fluxCentroid *= m_beta;
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
                 const EBCellFAB&              a_phi,
                 const Box&                    a_ghostedBox,
                 const Real&                   a_dx,
                 const DataIndex&              a_datInd,
                 const int&                    a_idir)
  {
    BL_PROFILE("ebco::getFluxRegOnly");
    const Box& domain = m_eblg.getDomain();

    //has some extra cells so...
    a_fluxCentroid.setVal(0.);
    int ncomp = a_phi.nComp();
    BL_ASSERT(ncomp == a_fluxCentroid.nComp());
    Box cellBox = a_ghostedBox;
    //want only interior faces
    cellBox.grow(a_idir, 1);
    cellBox &= domain;
    cellBox.grow(a_idir,-1);

    Box faceBox = surroundingNodes(cellBox, a_idir);

    //make a EBFaceFAB (including ghost cells) that will hold centered gradients
    BaseFab<Real>& regFlux = a_fluxCentroid.getSingleValuedFAB();
    const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
    const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
    const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());

    FORT_GETFLUXEBCO(CHF_FRA1(regFlux,0),
                     CHF_CONST_FRA1(regBCo,0),
                     CHF_CONST_FRA1(regPhi, 0),
                     CHF_BOX(faceBox),
                     CHF_CONST_REAL(a_dx),
                     CHF_CONST_INT(a_idir));

    a_fluxCentroid *= m_beta;
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  incrementFRFine(EBFastFR&             a_fluxReg,
                  const FabArray<EBCellFAB>& a_phiFine,
                  const FabArray<EBCellFAB>& a_phi,
                  AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("NWOEBConductivityOp::incrementFRFine");
    BL_ASSERT(a_phiFine.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);
    BL_ASSERT(m_hasFine);
    int ncomp = 1;
    Interval interv(0,0);
    NWOEBConductivityOp& finerEBAMROp = (NWOEBConductivityOp& )(*a_finerOp);

    //ghost cells of phiFine need to be filled
    FabArray<EBCellFAB>& phiFine = (FabArray<EBCellFAB>&) a_phiFine;
    finerEBAMROp.m_interpWithCoarser->coarseFineInterp(phiFine, a_phi, 0, 0, 1);
    phiFine.exchange(finerEBAMROp.m_exchangeCopier);

    DataIterator ditf = a_phiFine.dataIterator();
    int nbox = ditf.size();
#pragma omp parallel for
    for(int mybox=0; mybox<nbox; mybox++)
    {
      const Box&     boxFine = m_eblgFine.getDBL().get(ditf[mybox]);
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[ditf[mybox]];
      const EBCellFAB& phiFine = a_phiFine[ditf[mybox]];

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        for (SideIterator sit; sit.ok(); sit.next())
        {
          Box fabBox = adjCellBox(boxFine, idir, sit(), 1);
          fabBox.shift(idir, -sign(sit()));

          Box ghostedBox = fabBox;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblgFine.getDomain();

          EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);
          finerEBAMROp.getFlux(fluxFine, phiFine, ghostedBox, fabBox, m_eblgFine.getDomain(),
                               ebisBoxFine, m_dxFine, ditf[mybox], idir);

          Real scale = 1.0; //beta and bcoef already in flux

          a_fluxReg.incrementFineBoth(fluxFine, scale, ditf[mybox], interv, idir, sit());
        }
      }
    }
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  getFlux(FArrayBox&                    a_flux,
          const FArrayBox&              a_phi,
          const Box&                    a_faceBox,
          const int&                    a_idir,
          const Real&                   a_dx,
          const DataIndex&              a_datInd)
  {
    BL_PROFILE("ebco::getflux3");
    const EBFaceFAB& bcoebff = (*m_bcoef)[a_datInd][a_idir];
    const FArrayBox& regBCo = (const FArrayBox&)(bcoebff.getSingleValuedFAB());
    FORT_GETFLUXEBCO(CHF_FRA1(a_flux,0),
                     CHF_CONST_FRA1(regBCo,0),
                     CHF_CONST_FRA1(a_phi, 0),
                     CHF_BOX(a_faceBox),
                     CHF_CONST_REAL(a_dx),
                     CHF_CONST_INT(a_idir));

    a_flux *= m_beta;
  }
//-----------------------------------------------------------------------

  void
  NWOEBConductivityOp::
  AMRResidualNC(FabArray<EBCellFAB>&       a_residual,
                const FabArray<EBCellFAB>& a_phiFine,
                const FabArray<EBCellFAB>& a_phi,
                const FabArray<EBCellFAB>& a_rhs,
                bool a_homogeneousPhysBC,
                AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    //dummy. there is no coarse when this is called
    BL_ASSERT(a_residual.nGrow() == m_ghostRHS);
    BL_ASSERT(a_rhs.nGrow() == m_ghostRHS);
    FabArray<EBCellFAB> phiC;
    AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousPhysBC, a_finerOp);
  }
//-----------------------------------------------------------------------

  void
  NWOEBConductivityOp::
  AMRResidualNF(FabArray<EBCellFAB>&       a_residual,
                const FabArray<EBCellFAB>& a_phi,
                const FabArray<EBCellFAB>& a_phiCoar,
                const FabArray<EBCellFAB>& a_rhs,
                bool a_homogeneousPhysBC)
  {
    BL_PROFILE("ebco::amrresNF");
    BL_ASSERT(a_residual.nGrow() == m_ghostRHS);
    BL_ASSERT(a_rhs.nGrow() == m_ghostRHS);

    AMROperatorNF(a_residual, a_phi, a_phiCoar,
                  a_homogeneousPhysBC);
    axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  }
//-----------------------------------------------------------------------

  void
  NWOEBConductivityOp::
  AMROperatorNC(FabArray<EBCellFAB>&       a_LofPhi,
                const FabArray<EBCellFAB>& a_phiFine,
                const FabArray<EBCellFAB>& a_phi,
                bool a_homogeneousPhysBC,
                AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("ebco::amrOpNC");
    //dummy. there is no coarse when this is called
    BL_ASSERT(a_LofPhi.nGrow() == m_ghostRHS);
    FabArray<EBCellFAB> phiC;
    AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
                a_homogeneousPhysBC, a_finerOp);
  }
//-----------------------------------------------------------------------

  void
  NWOEBConductivityOp::
  AMROperatorNF(FabArray<EBCellFAB>&       a_LofPhi,
                const FabArray<EBCellFAB>& a_phi,
                const FabArray<EBCellFAB>& a_phiCoar,
                bool a_homogeneousPhysBC)
  {
    BL_ASSERT(a_LofPhi.nGrow() == m_ghostRHS);

    applyOp(a_LofPhi,a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);

  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  AMRRestrict(FabArray<EBCellFAB>&       a_resCoar,
              const FabArray<EBCellFAB>& a_residual,
              const FabArray<EBCellFAB>& a_correction,
              const FabArray<EBCellFAB>& a_coarCorrection, 
              bool a_skip_res )
  {
    BL_PROFILE("NWOEBConductivityOp::AMRRestrict");
    BL_ASSERT(a_residual.nGrow() == m_ghostRHS);
    BL_ASSERT(a_correction.nGrow() == m_ghostPhi);
    BL_ASSERT(a_coarCorrection.nGrow() == m_ghostPhi);
    BL_ASSERT(!a_skip_res);

    BL_ASSERT(a_residual.nComp() == 1);
    BL_ASSERT(a_resCoar.nComp() == 1);
    BL_ASSERT(a_correction.nComp() == 1);

    FabArray<EBCellFAB> resThisLevel;
    bool homogeneousPhys = true;
    bool homogeneousCF =   false;

    EBCellFactory ebcellfactTL(m_eblg.getEBISL());
    IntVect ghostVec = a_residual.nGrow();

    resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);
    EBFabArrayOps::setVal(resThisLevel, 0.0);

    //API says that we must average(a_residual - L(correction, coarCorrection))
    applyOp(resThisLevel, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);
    incr(resThisLevel, a_residual, -1.0);
    scale(resThisLevel,-1.0);

    //use our nifty averaging operator
    Interval variables(0, 0);
    m_ebAverage.average(a_resCoar, resThisLevel, variables);

  }
//-----------------------------------------------------------------------
  Real
  NWOEBConductivityOp::
  AMRNorm(const FabArray<EBCellFAB>& a_coarResid,
          const FabArray<EBCellFAB>& a_fineResid,
          const int& a_refRat,
          const int& a_ord)

  {
    // compute norm over all cells on coarse not covered by finer
    BL_PROFILE("NWOEBConductivityOp::AMRNorm");
    MayDay::Error("never called");
    //return norm of temp
    return norm(a_coarResid, a_ord);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  AMRProlong(FabArray<EBCellFAB>&       a_correction,
             const FabArray<EBCellFAB>& a_coarCorrection)
  {
    BL_PROFILE("NWOEBConductivityOp::AMRProlong");
    //use cached interpolation object
    Interval variables(0, 0);
    m_ebInterp.pwcInterp(a_correction, a_coarCorrection, variables);
  }
//-----------------------------------------------------------------------
  void
  NWOEBConductivityOp::
  AMRUpdateResidual(FabArray<EBCellFAB>&       a_residual,
                    const FabArray<EBCellFAB>& a_correction,
                    const FabArray<EBCellFAB>& a_coarCorrection)
  {
    BL_PROFILE("NWOEBConductivityOp::AMRUpdateResidual");
    BL_ASSERT(a_residual.nGrow()   == m_ghostRHS);
    BL_ASSERT(a_correction.nGrow() == m_ghostPhi);
    BL_ASSERT(a_coarCorrection.nGrow() == m_ghostPhi);

    FabArray<EBCellFAB> lcorr;
    bool homogeneousPhys = true;
    bool homogeneousCF   = false;

    EBCellFactory ebcellfactTL(m_eblg.getEBISL());
    IntVect ghostVec = a_residual.nGrow();

    lcorr.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

    applyOp(lcorr, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);

    incr(a_residual, lcorr, -1);
  }
//-----------------------------------------------------------------------

}
