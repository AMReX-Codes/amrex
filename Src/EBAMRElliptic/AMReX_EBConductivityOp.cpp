#include "AMReX_EBConductivityOp.H"
#include "AMReX_EBEllipticFort_F.H"
#include "AMReX_IrregFABFactory.H"
#include <AMReX_ParallelDescriptor.H>
namespace amrex
{

  //-------------------------------------------------------------------------------
  EBConductivityOp::
  EBConductivityOp(const EBLevelGrid &                          a_eblgFine,
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
                   const int    &                               a_ghost)
  : AMRLevelOp<FabArray<EBCellFAB> >(),
    m_ghost(a_ghost),
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
    m_colors()
  {
    BL_PROFILE("EBConductivityOp::ConductivityOp");

    m_turnOffBCs = false; //REALLY needs to default to false

    if (m_hasFine)
    {
      m_eblgFine       = a_eblgFine;
      m_dxFine         = m_dx/a_refToFine;
    }

    EBCellFactory fact(m_eblg.getEBISL());
    m_relCoef.define(m_eblg.getDBL(), m_eblg.getDM(), 1, 0, MFInfo(), fact);
    if (m_hasCoar)
    {
      m_eblgCoar       = a_eblgCoar;
      bool useKappaWeighting = false; //we are averaging stuff that is already kappa weighted
      m_ebInterp.define(  m_eblg, m_eblgCoar,  m_refToCoar, m_ghost);
      m_ebAverage.define( m_eblg, m_eblgCoar,  m_refToCoar, useKappaWeighting, m_ghost);
      int numGhostToFill= 1; //no EB/CF crossing so only need one cell
      m_cfInterp.define(  m_eblg, m_eblgCoar,  m_refToCoar, m_ghost, numGhostToFill);
    }

    if (m_hasCoarMG)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg, m_eblgCoarMG,  mgRef, m_ghost);
      bool useKappaWeighting = false; //we are averaging stuff that is already kappa weighted
      m_ebAverageMG.define(m_eblg, m_eblgCoarMG,  mgRef, useKappaWeighting, m_ghost);
    }

    //define stencils for the operator
    defineStencils();

  }
  //-------------------------------------------------------------------------------
  Real
  EBConductivityOp::
  getSafety()
  {
    Real safety = 1.0;
    return safety;
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  calculateAlphaWeight()
  {
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      VoFIterator& vofit = m_vofIterIrreg[mfi];
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
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  getDivFStencil(VoFStencil&      a_vofStencil,
                 const VolIndex&  a_vof,
                 const MFIter  &  a_mfi)
  {
    BL_PROFILE("EBConductivityOp::getDivFStencil");
    const EBISBox& ebisBox = m_eblg.getEBISL()[a_mfi];
    a_vofStencil.clear();
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
        for (int iface = 0; iface < faces.size(); iface++)
        {
          VoFStencil fluxStencil;
          getFluxStencil(fluxStencil, faces[iface], a_mfi);
          Real areaFrac = ebisBox.areaFrac(faces[iface]);
          fluxStencil *= Real(isign)*areaFrac/m_dx;
          a_vofStencil += fluxStencil;
        }
      }
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  getFluxStencil(VoFStencil&      a_fluxStencil,
                 const FaceIndex& a_face,
                 const MFIter   & a_mfi)
  {
    /// stencil for flux computation.   the truly ugly part of this computation
    /// beta and eta are multiplied in here

    BL_PROFILE("EBConductivityOp::getFluxStencil");
    //need to do this by interpolating to centroids
    //so get the stencil at each face center and add with
    //interpolation weights
    FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                       (*m_eblg.getCFIVS())[a_mfi],
                                                       m_eblg.getEBISL()[a_mfi],
                                                       m_eblg.getDomain());

    a_fluxStencil.clear();
    for (int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, face, a_mfi);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                             const FaceIndex& a_face,
                             const MFIter   & a_mfi)
  {
    BL_PROFILE("EBConductivityOp::getFaceCenteredFluxStencil");
    //face centered gradient is just a centered diff
    int faceDir= a_face.direction();
    a_fluxStencil.clear();

    if (!a_face.isBoundary())
    {
      a_fluxStencil.add(a_face.getVoF(Side::Hi),  1.0/m_dx, 0);
      a_fluxStencil.add(a_face.getVoF(Side::Lo), -1.0/m_dx, 0);
      a_fluxStencil *= (*m_bcoef)[a_mfi][faceDir](a_face,0);
    }
    else
    {
      //the boundary condition handles this one.
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  setAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
  {
    BL_PROFILE("EBConductivityOp::setAlphaAndBeta");
    m_alpha = a_alpha;
    m_beta  = a_beta;
    calculateAlphaWeight(); //need to do this because the a coef has probably been changed under us
    calculateRelaxationCoefficient();
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  kappaScale(FabArray<EBCellFAB> & a_rhs)
  {
    BL_PROFILE("EBConductivityOp::kappaScale");
    EBLevelDataOps::kappaWeight(a_rhs);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  diagonalScale(FabArray<EBCellFAB> & a_rhs,
                bool a_kappaWeighted)
  {

    BL_PROFILE("EBConductivityOp::diagonalScale");
    if (a_kappaWeighted)
      EBLevelDataOps::kappaWeight(a_rhs);

    //also have to weight by the identity coefficient
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      a_rhs[mfi] *= (*m_acoef)[mfi];
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  divideByIdentityCoef(FabArray<EBCellFAB> & a_rhs)
  {

    BL_PROFILE("EBConductivityOp::divideByIdentityCoef");

    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      a_rhs[mfi] /= (*m_acoef)[mfi];
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  calculateRelaxationCoefficient()
  {
    BL_PROFILE("ebco::calculateRelCoef");
    // define regular relaxation coefficent
    Real safety = getSafety();
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const Box& grid = m_eblg.getDBL()[mfi];
      
      // For time-independent acoef, initialize lambda = alpha * acoef.
      const EBCellFAB& acofab = (*m_acoef)[mfi];
      m_relCoef[mfi].setVal(0.);
      m_relCoef[mfi].plus(acofab, 0, 0, 1);
      m_relCoef[mfi]*= m_alpha;
      
      // Compute the relaxation coefficient in regular cells.
      BaseFab<Real>& regRel = m_relCoef[mfi].getSingleValuedFAB();
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        BaseFab<Real>& regBCo = (*m_bcoef)[mfi][idir].getSingleValuedFAB();
        ebefnd_decrinvrelcoefebco(BL_TO_FORTRAN_FAB(regRel),
                                  BL_TO_FORTRAN_FAB(regBCo),
                                  BL_TO_FORTRAN_BOX(grid),
                                  &m_beta, &m_dx, &idir);
      }
      
      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      ebefnd_invertlambdaebco(BL_TO_FORTRAN_FAB(regRel),
                              BL_TO_FORTRAN_BOX(grid), &safety);
      
      // Now go over the irregular cells.
      VoFIterator& vofit = m_vofIterIrreg[mfi];
      for (vofit.reset(); vofit.ok(); ++vofit)
      {
        const VolIndex& VoF = vofit();
          
        Real alphaWeight = m_alphaDiagWeight[mfi](VoF, 0);
        Real  betaWeight =  m_betaDiagWeight[mfi](VoF, 0);
        alphaWeight *= m_alpha;
        betaWeight  *= m_beta;
          
        Real diagWeight  = alphaWeight + betaWeight;
          
        // Set the irregular relaxation coefficients.
        if (std::abs(diagWeight) > 1.0e-9)
        {
          m_relCoef[mfi](VoF, 0) = safety/diagWeight;
        }
        else
        {
          m_relCoef[mfi](VoF, 0) = 0.;
        }
      }
    }
  }

  //-------------------------------------------------------------------------------
 
  void null_deleter_ebco_vof(VolIndex *)
  {}
  void null_deleter_ebco_sten(VoFStencil *)
  {}

  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  defineStencils()
  {
    BL_PROFILE("EBConductivityOp::defineStencils");
    // create ebstencil for irregular applyOp
    m_opEBStencil.define(m_eblg.getDBL(), m_eblg.getDM());
    // create vofstencils for applyOp and

    Real fakeBeta = 1;
    m_domainBC->setCoef(m_eblg,   fakeBeta ,  m_dx, RealVect::Zero,    m_bcoef);
    m_ebBC->setCoef(    m_eblg,   fakeBeta ,  m_dx, RealVect::Zero,    m_bcoef);

    LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);
    IrregFABFactory irregFact(m_eblg.getEBISL());
    m_vofIterIrreg.define(    m_eblg.getDBL(), m_eblg.getDM()); // vofiterator cache
    m_vofIterMulti.define(    m_eblg.getDBL(), m_eblg.getDM()); // vofiterator cache
    m_alphaDiagWeight.define( m_eblg.getDBL(), m_eblg.getDM(), 1, 0, MFInfo(), irregFact);
    m_betaDiagWeight.define(  m_eblg.getDBL(), m_eblg.getDM(), 1, 0, MFInfo(), irregFact);
    Box sideBoxLo[SpaceDim];
    Box sideBoxHi[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box domainBox = m_eblg.getDomain();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_vofIterDomLo[idir].define( m_eblg.getDBL(), m_eblg.getDM()); // vofiterator cache for domain lo
      m_vofIterDomHi[idir].define( m_eblg.getDBL(), m_eblg.getDM()); // vofiterator cache for domain hi
    }
    EBArith::getMultiColors(m_colors);

    EBCellFactory fact(m_eblg.getEBISL());              
    FabArray<EBCellFAB> phiProxy(m_eblg.getDBL(), m_eblg.getDM(), SpaceDim, m_ghost, MFInfo(), fact);
    FabArray<EBCellFAB> rhsProxy(m_eblg.getDBL(), m_eblg.getDM(), SpaceDim, m_ghost, MFInfo(), fact);

    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const Box& curBox = m_eblg.getDBL()[mfi];
      const EBISBox& ebisBox = m_eblg.getEBISL()[mfi];
      const EBGraph& ebgraph = ebisBox.getEBGraph();

      IntVectSet irregIVS = ebisBox.getIrregIVS(curBox);
      IntVectSet multiIVS = ebisBox.getMultiCells(curBox);

      //cache the vofIterators
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
        VolIndex  & vof     = (VolIndex &)(vofvec[ivec]);
        VoFStencil& vofsten = stenvec[ivec];

        //bcoef is included here in the flux consistent
        //with the regular
//begin debug        
        IntVect ivdebug(D_DECL(15,12,0));
        int idebug = 0;
        if(vof.gridIndex() == ivdebug)
        {
          idebug = 1;
        }
//end debug
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
            vector<FaceIndex> faces = ebisBox.getFaces(vof,idir,Side::Lo);
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
            vector<FaceIndex> faces = ebisBox.getFaces(vof,idir,Side::Hi);
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
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp);
    }

  }

  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  residual(FabArray<EBCellFAB>&       a_residual,
           const FabArray<EBCellFAB>& a_phi,
           const FabArray<EBCellFAB>& a_rhs,
           bool                        a_homogeneousPhysBC)
  {
    BL_PROFILE("EBConductivityOp::residual");
    //this is a EBMultiGrid operator so only homogeneous CF BC
    //and null coar level
    BL_ASSERT(a_residual.nGrow() == m_ghost);
    BL_ASSERT(a_phi.nGrow() == m_ghost);
    applyOp(a_residual,a_phi,NULL, a_homogeneousPhysBC, true);
    axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  preCond(FabArray<EBCellFAB>&       a_lhs,
          const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("EBConductivityOp::preCond");
    EBLevelDataOps::assign(a_lhs, a_rhs, 1.0);
    EBLevelDataOps::scale(a_lhs, m_relCoef);

    relax(a_lhs, a_rhs, 40);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
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
        m_cfInterp.coarseFineInterpH(phi,  0, 0, 1);
      }
      else
      {
        m_cfInterp.coarseFineInterp(phi, *a_phiCoar,  0, 0, 1);
      }
    }
    applyOp(a_lhs, a_phi, a_homogeneousPhysBC);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  applyOp(FabArray<EBCellFAB>&                    a_lhs,
          const FabArray<EBCellFAB>&              a_phi,
          bool                                    a_homogeneous)
  {
    BL_PROFILE("ebco::applyOp");
    FabArray<EBCellFAB>&  phi = (FabArray<EBCellFAB>&) a_phi;
    {
      BL_PROFILE("ghostcell fills");
      fillPhiGhost(phi, a_homogeneous);
    }

    {
      BL_PROFILE("applying op without bcs");
      for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
      {
        applyOpRegular(  a_lhs[mfi], a_phi[mfi], a_homogeneous, mfi);
        applyOpIrregular(a_lhs[mfi], a_phi[mfi], a_homogeneous, mfi);
      }
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  fillPhiGhost(const FabArray<EBCellFAB>& a_phi, bool a_homog) const
  {
    BL_PROFILE("nwoebco::fillghostLD");
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      fillPhiGhost(a_phi[mfi], mfi, a_homog);
    }
    FabArray<EBCellFAB>& phi = (FabArray<EBCellFAB>&)(a_phi);
    phi.FillBoundary();
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  fillPhiGhost(const EBCellFAB& a_phi, const MFIter   & a_datInd, bool a_homog) const
  {
    BL_PROFILE("nwoebco::fillghostfab");

    EBCellFAB& phi = (EBCellFAB&) a_phi;
    ConductivityBaseDomainBC* condbc = dynamic_cast<ConductivityBaseDomainBC*>(&(*m_domainBC));
    Box grid = m_eblg.getDBL()[a_datInd];
    Box domBox = m_eblg.getDomain();
    if(condbc == NULL)
    {
      amrex::Error("dynamic cast failed");
    }
    if (!m_turnOffBCs)
    {
      FArrayBox& fab = phi.getFArrayBox();
      condbc->fillPhiGhost(fab, grid, a_homog);
    }
    else
    {
      Box valid = m_eblg.getDBL()[a_datInd];
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        for(SideIterator sit; sit.ok(); ++sit)
        {
          BaseFab<Real>& fab = phi.getSingleValuedFAB();
          EBArith::ExtrapolateBC(fab, valid,  m_dx, idir, sit());
        }
        //grow so that we hit corners
        valid.grow(idir, 1);
      }
    }

  }
/*****/
  void
  EBConductivityOp::
  applyOpRegular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const MFIter   &       a_datInd)
  {
    BL_PROFILE("nwoebco::applyopRegular");
    //assumes ghost cells already filled

    const Box& grid = m_eblg.getDBL()[a_datInd];
    BaseFab<Real>      & lphfab  =  a_lhs.getSingleValuedFAB();
    const BaseFab<Real>& phifab  =  a_phi.getSingleValuedFAB();
    const BaseFab<Real>& acofab  =(*m_acoef)[a_datInd].getSingleValuedFAB();
    const EBFluxFAB& bco  =   (*m_bcoef)[a_datInd];

    vector<const BaseFab<Real>* > bcoside(3, &(bco[0].getSingleValuedFAB()));
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      bcoside[idir] = &(bco[idir].getSingleValuedFAB());
    }

    
    ebefnd_applyop_ebcond_nobcs(BL_TO_FORTRAN_FAB(lphfab),
                                BL_TO_FORTRAN_FAB(phifab),
                                BL_TO_FORTRAN_FAB(acofab),
                                BL_TO_FORTRAN_FAB((*bcoside[0])),
                                BL_TO_FORTRAN_FAB((*bcoside[1])),
                                BL_TO_FORTRAN_FAB((*bcoside[2])),
                                BL_TO_FORTRAN_BOX(grid),
                                &m_dx, &m_alpha, &m_beta);
  }
  //-------------------------------------------------------------------------------

  void
  EBConductivityOp::
  applyOpIrregular(EBCellFAB&             a_lhs,
                   const EBCellFAB&       a_phi,
                   const bool&            a_homogeneous,
                   const MFIter   &       a_datInd)
  {
    BL_PROFILE("nwoebco::applyOpIrregular");
    RealVect vectDx = m_dx*RealVect::Unit;

    m_opEBStencil[a_datInd]->apply(a_lhs, a_phi, m_alphaDiagWeight[a_datInd], m_alpha, m_beta, 0, false);
  
    if (!a_homogeneous)
    {
      const Real factor = m_beta/m_dx;
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIterIrreg[a_datInd].getVector(), 
                          a_datInd, factor, a_homogeneous);
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      int comp = 0;
      for (m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
      {
        Real flux;
        const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
        m_domainBC->getFaceFlux(flux,vof, a_datInd, a_phi,
                                idir,Side::Lo, a_homogeneous);

        //area gets multiplied in by bc operator
        a_lhs(vof,comp) -= flux*m_beta/m_dx;
      }
      for (m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
      {
        Real flux;
        const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
        m_domainBC->getFaceFlux(flux,vof, a_datInd, a_phi,
                                idir,Side::Hi, a_homogeneous);

        //area gets multiplied in by bc operator
        a_lhs(vof,comp) += flux*m_beta/m_dx;
      }
    }

  }

  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  create(FabArray<EBCellFAB>&       a_lhs,
         const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("ebco::create");
    int ncomp = a_rhs.nComp();
    EBCellFactory ebcellfact(m_eblg.getEBISL());
    a_lhs.clear();
    a_lhs.define(m_eblg.getDBL(), m_eblg.getDM(), ncomp, a_rhs.nGrow(), MFInfo(), ebcellfact);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  createCoarsened(FabArray<EBCellFAB>&       a_lhs,
                  const FabArray<EBCellFAB>& a_rhs,
                  const int &                a_refRat)
  {
    BL_PROFILE("ebco::createCoar");
    int ncomp = a_rhs.nComp();
    int ghost = a_rhs.nGrow();

    BL_ASSERT(m_eblg.getDBL().coarsenable(a_refRat));

    EBLevelGrid eblgCoFi;
    coarsen(eblgCoFi, m_eblg, a_refRat);

    EBCellFactory ebcellfactCoFi(eblgCoFi.getEBISL());
    a_lhs.define(eblgCoFi.getDBL(), eblgCoFi.getDM(), ncomp, ghost, MFInfo(), ebcellfactCoFi);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  assign(FabArray<EBCellFAB>&       a_lhs,
         const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("ebco::assign");
    EBLevelDataOps::assign(a_lhs,a_rhs, 1.0);
  }
  //-------------------------------------------------------------------------------
  Real
  EBConductivityOp::
  dotProduct(const FabArray<EBCellFAB>& a_1,
             const FabArray<EBCellFAB>& a_2)
  {
    BL_PROFILE("ebco::dotProd");
    Box domain;
    Real volume;

    return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,m_eblg);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  incr(FabArray<EBCellFAB>&       a_lhs,
       const FabArray<EBCellFAB>& a_x,
       Real                        a_scale)
  {
    BL_PROFILE("ebco::incr");
    EBLevelDataOps::incr(a_lhs,a_x,a_scale);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  axby(FabArray<EBCellFAB>&       a_lhs,
       const FabArray<EBCellFAB>& a_x,
       const FabArray<EBCellFAB>& a_y,
       Real                        a_a,
       Real                        a_b)
  {
    BL_PROFILE("ebco::axby");
    EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  scale(FabArray<EBCellFAB>& a_lhs,
        const Real&           a_scale)
  {
    BL_PROFILE("ebco::scale");
    EBLevelDataOps::scale(a_lhs,a_scale);
  }
  //---------------------------------------
  Real 
  EBConductivityOp::
  norm(const FabArray<EBCellFAB>& a_rhs,
       int                         a_ord)
  {
    BL_PROFILE("EBConductivityOp::norm");

    Real maxNorm = 0.0;

    maxNorm = localMaxNorm(a_rhs);

#ifdef BL_USE_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&maxNorm, &tmp, 1,
                               ParallelDescriptor::Mpi_typemap<Real>::type(),
                               MPI_MAX, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS)
    { //bark!!!
      amrex::Error("sorry, but I had a communcation error on norm");
    }
    maxNorm = tmp;
#endif
    //  Real volume=1.;
    //  EBLevelDataOps::gatherBroadCast(maxNorm, volume, 0);

    return maxNorm;
  }

  Real 
  EBConductivityOp::
  localMaxNorm(const FabArray<EBCellFAB>& a_rhs)
  {
    BL_PROFILE("EBAMRPoissonOp::localMaxNorm");
    Real maxval, minval;
    EBLevelDataOps::getMaxMin(maxval, minval, a_rhs, 0, true);
    return maxval;
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  setToZero(FabArray<EBCellFAB>& a_lhs)
  {
    BL_PROFILE("ebco::setToZero");
    EBLevelDataOps::setVal(a_lhs, 0.0);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  setVal(FabArray<EBCellFAB>& a_lhs, const Real& a_value)
  {
    BL_PROFILE("ebco::setVal");
    EBLevelDataOps::setVal(a_lhs, a_value);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  createCoarser(FabArray<EBCellFAB>&       a_coar,
                const FabArray<EBCellFAB>& a_fine,
                bool                        a_ghosted)
  {
    BL_PROFILE("ebco::createCoarser");

    EBCellFactory ebcellfact(m_eblgCoarMG.getEBISL());
    a_coar.define(m_eblgCoarMG.getDBL(), m_eblgCoarMG.getDM(), a_fine.nComp(), a_fine.nGrow(), MFInfo(), ebcellfact);
  }
  //--------------------------------------------------------------------------------
/*****/
  void
  EBConductivityOp::
  relax(FabArray<EBCellFAB>&       a_phi,
        const FabArray<EBCellFAB>& a_rhs,
        int                         a_iterations)
  {
    BL_PROFILE("nwoebco::relax");
    BL_ASSERT(a_phi.nGrow() == m_ghost);
    BL_ASSERT(a_rhs.nGrow() == m_ghost);
    BL_ASSERT(a_phi.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == a_rhs.nComp());

    // do first red, then black passes
    for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for (int icolor = 0; icolor < m_colors.size(); icolor++)
      {
        //after this lphi = L(phi)
        //this call contains bcs and exchange
        BL_PROFILE("ghostfill and homogcfinterp");
        fillPhiGhost(a_phi, true);
        homogeneousCFInterp(a_phi);

        gsrbColor(a_phi, a_rhs, m_colors[icolor]);
      }
    }
  }

  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  homogeneousCFInterp(FabArray<EBCellFAB>&   a_phif)
  {
    BL_PROFILE("nwoebco::homog_cfinterp");
    if (m_hasCoar)
    {
      m_cfInterp.coarseFineInterpH(a_phif,  0, 0, 1);
    }
  }

  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  gsrbColor(FabArray<EBCellFAB>&       a_phi,
            const FabArray<EBCellFAB>& a_rhs,
            const IntVect&             a_color)
  {

    BL_PROFILE("nwoebco::gsrbColor");

    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      //first do the the regular stuff
      Box dblBox  =   m_eblg.getDBL()[mfi];
      BaseFab<Real>&       regPhi =     a_phi[mfi].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[mfi].getSingleValuedFAB();
      const BaseFab<Real>& regRel = m_relCoef[mfi].getSingleValuedFAB();


      //assumes ghost cells already filled

      const BaseFab<Real>& acofab  =(*m_acoef)[mfi].getSingleValuedFAB();
      const EBFluxFAB& bco  =       (*m_bcoef)[mfi];

      vector<const BaseFab<Real>* > bcoside(3, &(bco[0].getSingleValuedFAB()));
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
        
      m_opEBStencil[mfi]->cachePhi(a_phi[mfi]);

      if (loIV <= hiIV)
      {
        Box coloredBox(loIV, hiIV);
        ebefnd_gscolor_ebcond(BL_TO_FORTRAN_FAB(regPhi),
                              BL_TO_FORTRAN_FAB(regRhs),
                              BL_TO_FORTRAN_FAB(regRel),
                              BL_TO_FORTRAN_FAB(acofab),
                              BL_TO_FORTRAN_FAB((*bcoside[0])),
                              BL_TO_FORTRAN_FAB((*bcoside[1])),
                              BL_TO_FORTRAN_FAB((*bcoside[2])),
                              BL_TO_FORTRAN_BOX(coloredBox),
                              &m_dx, &m_alpha, &m_beta);
      }

      m_opEBStencil[mfi]->uncachePhi(a_phi[mfi]);


      m_opEBStencil[mfi]->relax(a_phi[mfi], a_rhs[mfi], 
                                   m_relCoef[mfi], m_alphaDiagWeight[mfi],
                                   m_alpha, m_beta, 0, a_color);
    }
  }

  //-------------------------------------------------------------------------------
  void EBConductivityOp::
  restrictResidual(FabArray<EBCellFAB>&       a_resCoar,
                   FabArray<EBCellFAB>&       a_phiThisLevel,
                   const FabArray<EBCellFAB>& a_rhsThisLevel)
  {
    BL_PROFILE("EBConductivityOp::restrictResidual");

    BL_ASSERT(a_resCoar.nComp() == 1);
    BL_ASSERT(a_phiThisLevel.nComp() == 1);
    BL_ASSERT(a_rhsThisLevel.nComp() == 1);

    FabArray<EBCellFAB> resThisLevel;
    bool homogeneous = true;

    EBCellFactory ebcellfactTL(m_eblg.getEBISL());
    int ghost = a_rhsThisLevel.nGrow();

    resThisLevel.define(m_eblg.getDBL(), m_eblg.getDM(), 1, ghost, MFInfo(), ebcellfactTL);

    // Get the residual on the fine grid
    residual(resThisLevel,a_phiThisLevel,a_rhsThisLevel,homogeneous);

    // now use our nifty averaging operator
    m_ebAverageMG.average(a_resCoar, resThisLevel, 0, 0, 1);

  }
  //-------------------------------------------------------------------------------
  void EBConductivityOp::
  prolongIncrement(FabArray<EBCellFAB>&       a_phiThisLevel,
                   const FabArray<EBCellFAB>& a_correctCoar)
  {
    BL_PROFILE("EBConductivityOp::prolongIncrement");
    m_ebInterpMG.interpolate(a_phiThisLevel, a_correctCoar, 0, 0, 1);
  }
  //-------------------------------------------------------------------------------
  int EBConductivityOp::
  refToCoarser()
  {
    return m_refToCoar;
  }
  //-------------------------------------------------------------------------------
  int EBConductivityOp::
  refToFiner()
  {
    return m_refToFine;
  }
  //-------------------------------------------------------------------------------
  void EBConductivityOp::
  AMRResidual(FabArray<EBCellFAB>&       a_residual,
              const FabArray<EBCellFAB>& a_phiFine,
              const FabArray<EBCellFAB>& a_phi,
              const FabArray<EBCellFAB>& a_phiCoar,
              const FabArray<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC,
              AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("EBConductivityOp::AMRResidual");
    BL_ASSERT(a_residual.nGrow() == m_ghost);
    BL_ASSERT(a_rhs.nGrow() == m_ghost);
    BL_ASSERT(a_residual.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);
    BL_ASSERT(a_rhs.nComp() == 1);

    AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar,
                a_homogeneousPhysBC, a_finerOp);

    //multiply by -1 so a_residual now holds -L(phi)
    //add in rhs so a_residual = rhs - L(phi)
    axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  }
  //-------------------------------------------------------------------------------
  void EBConductivityOp::
  AMROperator(FabArray<EBCellFAB>&       a_LofPhi,
              const FabArray<EBCellFAB>& a_phiFine,
              const FabArray<EBCellFAB>& a_phi,
              const FabArray<EBCellFAB>& a_phiCoar,
              bool a_homogeneousPhysBC,
              AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("EBConductivityOp::AMROperator");
    BL_ASSERT(a_LofPhi.nGrow() == m_ghost);
    BL_ASSERT(a_LofPhi.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);

    //apply the operator between this and the next coarser level.
    applyOp(a_LofPhi, a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);


    //now reflux to enforce flux-matching from finer levels
    if (m_hasFine)
    {
      BL_ASSERT(a_finerOp != NULL);
      reflux(a_LofPhi, a_phiFine, a_phi, a_finerOp);
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  reflux(FabArray<EBCellFAB>& a_residual,
         const FabArray<EBCellFAB>& a_phiFine,
         const FabArray<EBCellFAB>& a_phi,
         AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("EBConductivityOp::reflux");

    m_fastFR.setToZero();

    incrementFRCoar(m_fastFR, a_phiFine, a_phi);

    incrementFRFine(m_fastFR, a_phiFine, a_phi, a_finerOp);


    Real scale = 1.0/m_dx;
    m_fastFR.reflux(a_residual, scale, 0, 0, 1);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  incrementFRCoar(EBFastFR&             a_fluxReg,
                  const FabArray<EBCellFAB>& a_phiFine,
                  const FabArray<EBCellFAB>& a_phi)
  {
    BL_PROFILE("EBConductivityOp::incrementFRCoar");
    BL_ASSERT(a_phiFine.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);

    int ncomp = 1;
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {

      const EBCellFAB & coarfab = a_phi[mfi];
      const EBISBox   & ebisBox = m_eblg.getEBISL()[mfi];
      const Box          &  box = m_eblg.getDBL()  [mfi];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        //no boundary faces here.

        Box ghostedBox = box;
        ghostedBox.grow(1);
        ghostedBox.grow(idir,-1);
        ghostedBox &= m_eblg.getDomain();

        EBFaceFAB coarflux(ebisBox, ghostedBox, idir, ncomp);

        // new way (because we know no EB/CF)
        fillPhiGhost(a_phi[mfi], mfi, false);
        getFluxRegOnly(coarflux, coarfab, ghostedBox, mfi, idir);

        Real scale = 1.0; //beta and bcoef already in flux
        for (SideIterator sit; sit.ok(); ++sit)
        {
          a_fluxReg.incrementCoarse(coarflux, scale, mfi, 0, 0, 1);
        }
      }
    }
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
                 const EBCellFAB&              a_phi,
                 const Box&                    a_ghostedBox,
                 const MFIter   &              a_mfi,
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
    BaseFab<Real>& regFlux      = a_fluxCentroid.getSingleValuedFAB();
    const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
    const EBFaceFAB& bcoebff    = (*m_bcoef)[a_mfi][a_idir];
    const BaseFab<Real>& regBCo = bcoebff.getSingleValuedFAB();

    ebefnd_getflux_ebco(BL_TO_FORTRAN_FAB(regFlux),
                        BL_TO_FORTRAN_FAB(regPhi),
                        BL_TO_FORTRAN_FAB(regBCo),
                        BL_TO_FORTRAN_BOX(faceBox),
                        &m_beta, &m_dx, &a_idir);

  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  incrementFRFine(EBFastFR&             a_fluxReg,
                  const FabArray<EBCellFAB>& a_phiFine,
                  const FabArray<EBCellFAB>& a_phi,
                  AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("EBConductivityOp::incrementFRFine");
    BL_ASSERT(a_phiFine.nComp() == 1);
    BL_ASSERT(a_phi.nComp() == 1);
    BL_ASSERT(m_hasFine);
    int ncomp = 1;

    EBConductivityOp& finerEBAMROp = (EBConductivityOp& )(*a_finerOp);

    //ghost cells of phiFine need to be filled
    FabArray<EBCellFAB>& phiFine = (FabArray<EBCellFAB>&) a_phiFine;
    finerEBAMROp.m_cfInterp.coarseFineInterp(phiFine, a_phi, 0, 0, 1);
    phiFine.FillBoundary();

    for(MFIter mfi_f(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi_f.isValid(); ++mfi_f)
    {
      const Box&     boxFine     = m_eblgFine.getDBL()  [mfi_f];
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[mfi_f];
      const EBCellFAB& phiFine = a_phiFine[mfi_f];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
//        for (SideIterator sit; sit.ok(); sit.next())
//        {
        Box ghostedBox = boxFine;
        ghostedBox.grow(1);
//        ghostedBox.grow(idir,-1);
          ghostedBox &= m_eblgFine.getDomain();

          EBFaceFAB fluxFine(ebisBoxFine, ghostedBox, idir, ncomp);


          //again, no EB/CF
          finerEBAMROp.getFluxRegOnly(fluxFine, phiFine, ghostedBox, mfi_f, idir);

          Real scale = 1.0; //beta and bcoef already in flux

          a_fluxReg.incrementFine(fluxFine, scale, mfi_f, 0, 0, 1);
//        }
      }
    }
  }
  //-------------------------------------------------------------------------------

  void
  EBConductivityOp::
  AMRResidualNC(FabArray<EBCellFAB>&       a_residual,
                const FabArray<EBCellFAB>& a_phiFine,
                const FabArray<EBCellFAB>& a_phi,
                const FabArray<EBCellFAB>& a_rhs,
                bool a_homogeneousPhysBC,
                AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    //dummy. there is no coarse when this is called
    BL_ASSERT(a_residual.nGrow() == m_ghost);
    BL_ASSERT(a_rhs.nGrow() == m_ghost);
    FabArray<EBCellFAB> phiC;
    AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousPhysBC, a_finerOp);
  }
  //-------------------------------------------------------------------------------

  void
  EBConductivityOp::
  AMRResidualNF(FabArray<EBCellFAB>&       a_residual,
                const FabArray<EBCellFAB>& a_phi,
                const FabArray<EBCellFAB>& a_phiCoar,
                const FabArray<EBCellFAB>& a_rhs,
                bool a_homogeneousPhysBC)
  {
    BL_PROFILE("ebco::amrresNF");
    BL_ASSERT(a_residual.nGrow() == m_ghost);
    BL_ASSERT(a_rhs.nGrow() == m_ghost);

    AMROperatorNF(a_residual, a_phi, a_phiCoar,
                  a_homogeneousPhysBC);
    axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
  }
  //-------------------------------------------------------------------------------

  void
  EBConductivityOp::
  AMROperatorNC(FabArray<EBCellFAB>&       a_LofPhi,
                const FabArray<EBCellFAB>& a_phiFine,
                const FabArray<EBCellFAB>& a_phi,
                bool a_homogeneousPhysBC,
                AMRLevelOp<FabArray<EBCellFAB> >* a_finerOp)
  {
    BL_PROFILE("ebco::amrOpNC");
    //dummy. there is no coarse when this is called
    BL_ASSERT(a_LofPhi.nGrow() == m_ghost);
    FabArray<EBCellFAB> phiC;
    AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
                a_homogeneousPhysBC, a_finerOp);
  }
  //-------------------------------------------------------------------------------

  void
  EBConductivityOp::
  AMROperatorNF(FabArray<EBCellFAB>&       a_LofPhi,
                const FabArray<EBCellFAB>& a_phi,
                const FabArray<EBCellFAB>& a_phiCoar,
                bool a_homogeneousPhysBC)
  {
    BL_ASSERT(a_LofPhi.nGrow() == m_ghost);

    applyOp(a_LofPhi,a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);

  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  AMRRestrict(FabArray<EBCellFAB>&       a_resCoar,
              const FabArray<EBCellFAB>& a_residual,
              const FabArray<EBCellFAB>& a_correction,
              const FabArray<EBCellFAB>& a_coarCorrection)
  {
    BL_PROFILE("EBConductivityOp::AMRRestrict");
    BL_ASSERT(a_residual.nGrow() == m_ghost);
    BL_ASSERT(a_correction.nGrow() == m_ghost);
    BL_ASSERT(a_coarCorrection.nGrow() == m_ghost);

    BL_ASSERT(a_residual.nComp() == 1);
    BL_ASSERT(a_resCoar.nComp() == 1);
    BL_ASSERT(a_correction.nComp() == 1);

    FabArray<EBCellFAB> resThisLevel;
    bool homogeneousPhys = true;
    bool homogeneousCF =   false;

    EBCellFactory ebcellfactTL(m_eblg.getEBISL());
    int ghostVec = a_residual.nGrow();

    resThisLevel.define(m_eblg.getDBL(), m_eblg.getDM(), 1, ghostVec, MFInfo(), ebcellfactTL);
    EBLevelDataOps::setVal(resThisLevel, 0.0);

    //API says that we must average(a_residual - L(correction, coarCorrection))
    applyOp(resThisLevel, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);
    incr(resThisLevel, a_residual, -1.0);
    scale(resThisLevel,-1.0);

    //use our nifty averaging operator
    m_ebAverage.average(a_resCoar, resThisLevel, 0, 0, 1);

  }
  //-------------------------------------------------------------------------------
  Real
  EBConductivityOp::
  AMRNorm(const FabArray<EBCellFAB>& a_coarResid,
          const FabArray<EBCellFAB>& a_fineResid,
          const int& a_refRat,
          const int& a_ord)

  {
    // compute norm over all cells on coarse not covered by finer
    BL_PROFILE("EBConductivityOp::AMRNorm");
    amrex::Error("never called");
    //return norm of temp
    return norm(a_coarResid, a_ord);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  AMRProlong(FabArray<EBCellFAB>&       a_correction,
             const FabArray<EBCellFAB>& a_coarCorrection)
  {
    BL_PROFILE("EBConductivityOp::AMRProlong");
    //use cached interpolation object
    m_ebInterp.interpolate(a_correction, a_coarCorrection, 0, 0, 1);
  }
  //-------------------------------------------------------------------------------
  void
  EBConductivityOp::
  AMRUpdateResidual(FabArray<EBCellFAB>&       a_residual,
                    const FabArray<EBCellFAB>& a_correction,
                    const FabArray<EBCellFAB>& a_coarCorrection)
  {
    BL_PROFILE("EBConductivityOp::AMRUpdateResidual");
    BL_ASSERT(a_residual.nGrow()   == m_ghost);
    BL_ASSERT(a_correction.nGrow() == m_ghost);
    BL_ASSERT(a_coarCorrection.nGrow() == m_ghost);

    FabArray<EBCellFAB> lcorr;
    bool homogeneousPhys = true;
    bool homogeneousCF   = false;

    EBCellFactory ebcellfactTL(m_eblg.getEBISL());
    int ghostVec = a_residual.nGrow();

    lcorr.define(m_eblg.getDBL(), m_eblg.getDM(), 1, ghostVec, MFInfo(), ebcellfactTL);

    applyOp(lcorr, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);

    incr(a_residual, lcorr, -1);
  }
  //-------------------------------------------------------------------------------

}
