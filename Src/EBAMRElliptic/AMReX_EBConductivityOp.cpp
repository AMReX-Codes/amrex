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
#include "EBAMRPoissonOp.H"
#include "NWOEBConductivityOp.H"
#include "EBQuadCFInterp.H"
#include "EBConductivityOpF_F.H"
#include "InterpF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "BCFunc.H"
#include "CH_Timer.H"
#include "BCFunc.H"
#include "EBLevelGrid.H"
#include "EBAMRPoissonOp.H"
#include "EBAMRPoissonOpF_F.H"
#include "EBAlias.H" 
#include "BCFunc.H"
#include "EBCoarseAverage.H"
#include "EBDebugOut.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"
bool NWOEBConductivityOp::s_doLazyRelax = false;

bool NWOEBConductivityOp::s_turnOffBCs = false; //REALLY needs to default to false
bool NWOEBConductivityOp::s_forceNoEBCF = false; //REALLY needs to default to false

//-----------------------------------------------------------------------
NWOEBConductivityOp::
NWOEBConductivityOp(const EBLevelGrid &                                  a_eblgFine,
                    const EBLevelGrid &                                  a_eblg,
                    const EBLevelGrid &                                  a_eblgCoar,
                    const EBLevelGrid &                                  a_eblgCoarMG,
                    const RefCountedPtr<NWOEBQuadCFInterp>&              a_quadCFI,
                    const RefCountedPtr<ConductivityBaseDomainBC>&       a_domainBC,
                    const RefCountedPtr<ConductivityBaseEBBC>&           a_ebBC,
                    const Real    &                                      a_dx,
                    const Real    &                                      a_dxCoar,
                    const int&                                           a_refToFine,
                    const int&                                           a_refToCoar,
                    const bool&                                          a_hasFine,
                    const bool&                                          a_hasCoar,
                    const bool&                                          a_hasMGObjects,
                    const bool&                                          a_layoutChanged,
                    const Real&                                          a_alpha,
                    const Real&                                          a_beta,
                    const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
                    const IntVect&                                       a_ghostPhi,
                    const IntVect&                                       a_ghostRHS,
                    const int&                                           a_relaxType)
: LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // is time-independent
  m_relaxType(a_relaxType),
  m_ghostPhi(a_ghostPhi),
  m_ghostRHS(a_ghostRHS),
  m_interpWithCoarser(a_quadCFI),
  m_eblg(a_eblg),
  m_eblgFine(),
  m_eblgCoar(),
  m_eblgCoarMG(),
  m_eblgCoarsenedFine(),
  m_domainBC(a_domainBC),
  m_ebBC(a_ebBC),
  m_dxFine(),
  m_dx(a_dx),
  m_dxCoar(a_dxCoar),
  m_acoef(a_acoef),
  m_bcoef(a_bcoef),
  m_bcoIrreg(a_bcoIrreg),
  m_alpha(a_alpha),
  m_beta(a_beta),
  m_alphaDiagWeight(),
  m_betaDiagWeight(),
  m_refToFine(a_hasFine ? a_refToFine : 1),
  m_refToCoar(a_hasCoar ? a_refToCoar : 1),
  m_hasFine(a_hasFine),
  m_hasInterpAve(false),
  m_hasCoar(a_hasCoar),
  m_ebAverage(),
  m_ebInterp(),
  m_opEBStencil(),
  m_relCoef(),
  m_vofIterIrreg(),
  m_vofIterMulti(),
  m_vofIterDomLo(),
  m_vofIterDomHi(),
  m_loCFIVS(),
  m_hiCFIVS(),
  m_fastFR(),
  m_hasMGObjects(a_hasMGObjects),
  m_layoutChanged(a_layoutChanged),
  m_ebAverageMG(),
  m_ebInterpMG(),
  m_dblCoarMG(),
  m_ebislCoarMG(),
  m_domainCoarMG(),
  m_colors()
{
  CH_TIME("NWOEBConductivityOp::ConductivityOp");
  int ncomp = 1;

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

      DataIterator dit = m_eblg.getDBL().dataIterator(); 
      int nbox = dit.size();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_loCFIVS[idir].define(m_eblg.getDBL());
          m_hiCFIVS[idir].define(m_eblg.getDBL());
#pragma omp parallel
          {
#pragma omp for
            for (int mybox=0;mybox<nbox;mybox++)
              {
                //                pout() << "doing lo cfivs for box " << mybox << endl;
                m_loCFIVS[idir][dit[mybox]].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit[mybox]),
                                                   m_eblg.getDBL(), idir,Side::Lo);
                //                pout() << "doing hi cfivs for box " << mybox << endl;
                m_hiCFIVS[idir][dit[mybox]].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit[mybox]),
                                                   m_eblg.getDBL(), idir,Side::Hi);
              }
          }
        }

      //if this fails, then the AMR grids violate proper nesting.
      ProblemDomain domainCoarsenedFine;
      DisjointBoxLayout dblCoarsenedFine;

      int maxBoxSize = 32;
      bool dumbool;
      bool hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarsenedFine,
                                                          domainCoarsenedFine,
                                                          m_eblg.getDBL(),
                                                          m_eblg.getEBISL(),
                                                          m_eblg.getDomain(),
                                                          m_refToCoar,
                                                          m_eblg.getEBIS(),
                                                          maxBoxSize, dumbool);

      //should follow from coarsenable
      if (hasCoarser)
        {
          m_eblgCoarsenedFine = EBLevelGrid(dblCoarsenedFine, domainCoarsenedFine, 4, m_eblg.getEBIS());
          m_hasInterpAve = true;
          m_ebInterp.define( m_eblg.getDBL(),     m_eblgCoar.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoar.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             a_ghostPhi);
          m_ebAverage.define(m_eblg.getDBL(),     m_eblgCoarsenedFine.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoarsenedFine.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, ncomp, m_eblg.getEBIS(),
                             a_ghostRHS);
        }
    }

  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  if (m_hasMGObjects)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain(), mgRef, ncomp, m_eblg.getEBIS(),
                           a_ghostPhi);
      m_ebAverageMG.define(m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, ncomp, m_eblg.getEBIS(),
                           a_ghostRHS);

    }

  m_exchangeCopier.define(m_eblg.getDBL(), m_eblg.getDBL(), m_ghostPhi,  true);

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
  int nbox = dit.size();
#pragma omp parallel
  {
#pragma omp for
    for(int mybox=0; mybox<nbox; mybox++)
      {

        VoFIterator& vofit = m_vofIterIrreg[dit[mybox]];
        for (vofit.reset(); vofit.ok(); ++vofit)
          {
            const VolIndex& VoF = vofit();
            Real volFrac = m_eblg.getEBISL()[dit[mybox]].volFrac(VoF);
            Real alphaWeight = (*m_acoef)[dit[mybox]](VoF, 0);
            alphaWeight *= volFrac;

            m_alphaDiagWeight[dit[mybox]](VoF, 0) = alphaWeight;
          }
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
  CH_TIME("NWOEBConductivityOp::getDivFStencil");
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
               const DataIndex& a_dit)
{
  /// stencil for flux computation.   the truly ugly part of this computation
  /// beta and eta are multiplied in here

  CH_TIME("NWOEBConductivityOp::getFluxStencil");
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
  CH_TIME("NWOEBConductivityOp::getFaceCenteredFluxStencil");
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
  CH_TIME("NWOEBConductivityOp::setAlphaAndBeta");
  m_alpha = a_alpha;
  m_beta  = a_beta;
  calculateAlphaWeight(); //need to do this because the a coef has probably been changed under us
  calculateRelaxationCoefficient();
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
kappaScale(LevelData<EBCellFAB> & a_rhs)
{
  CH_TIME("NWOEBConductivityOp::kappaScale");
  EBLevelDataOps::kappaWeight(a_rhs);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
diagonalScale(LevelData<EBCellFAB> & a_rhs,
              bool a_kappaWeighted)
{

  CH_TIME("NWOEBConductivityOp::diagonalScale");
  if (a_kappaWeighted)
    EBLevelDataOps::kappaWeight(a_rhs);

  //also have to weight by the coefficient
  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox = dit.size();
#pragma omp parallel
  {
#pragma omp for
    for(int mybox=0; mybox<nbox; mybox++)
      {
        a_rhs[dit[mybox]] *= (*m_acoef)[dit[mybox]];
      }
  }

}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs)
{

  CH_TIME("NWOEBConductivityOp::divideByIdentityCoef");

  DataIterator dit = m_eblg.getDBL().dataIterator(); 
  int nbox = dit.size();
#pragma omp parallel
  {
#pragma omp for
    for(int mybox=0; mybox<nbox; mybox++)
      {
        a_rhs[dit[mybox]] /= (*m_acoef)[dit[mybox]];
      }
  }
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
calculateRelaxationCoefficient()
{
  CH_TIME("ebco::calculateRelCoef");
  // define regular relaxation coefficent
  Real safety = getSafety();
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
#pragma omp parallel
  {
#pragma omp for
    for (int mybox=0;mybox<nbox;mybox++)
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
            FORT_DECRINVRELCOEFEBCO(CHF_FRA1(regRel,0),
                                    CHF_FRA1(regBCo,0),
                                    CHF_CONST_REAL(m_beta),
                                    CHF_BOX(grid),
                                    CHF_REAL(m_dx),
                                    CHF_INT(idir));
          }
      
        //now invert so lambda = stable lambda for variable coef lapl
        //(according to phil, this is the correct lambda)
        FORT_INVERTLAMBDAEBCO(CHF_FRA1(regRel,0),
                              CHF_REAL(safety),
                              CHF_BOX(grid));
      
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
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
defineStencils()
{
  CH_TIME("NWOEBConductivityOp::defineStencils");
  // create ebstencil for irregular applyOp
  m_opEBStencil.define(m_eblg.getDBL());
  // create vofstencils for applyOp and

  Real fakeBeta = 1;
  m_domainBC->setCoef(m_eblg,   fakeBeta ,      m_bcoef   );
  m_ebBC->setCoef(    m_eblg,   fakeBeta ,      m_bcoIrreg);

  Real dxScale = 1.0/m_dx;
  m_ebBC->define((*m_eblg.getCFIVS()), dxScale); //has to happen AFTER coefs are set
  LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(0);

  m_vofIterIrreg.define(     m_eblg.getDBL()); // vofiterator cache
  m_vofIterMulti.define(     m_eblg.getDBL()); // vofiterator cache
  m_alphaDiagWeight.define(  m_eblg.getDBL());
  m_betaDiagWeight.define(   m_eblg.getDBL());
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
  LevelData<EBCellFAB> phiProxy(m_eblg.getDBL(), SpaceDim, m_ghostPhi, fact);
  LevelData<EBCellFAB> rhsProxy(m_eblg.getDBL(), SpaceDim, m_ghostRHS, fact);
  DataIterator dit = m_eblg.getDBL().dataIterator(); 

  int nbox = dit.size();

  // this breaks  and it is a bad idea anyway (lots of dynamic memory allocation)
  //#pragma omp parallel for
  for(int mybox=0; mybox<nbox; mybox++)
    {
      const DataIndex& datind = dit[mybox];
      const Box& curBox = m_eblg.getDBL().get(datind);
      const EBISBox& ebisBox = m_eblg.getEBISL()[datind];
      const EBGraph& ebgraph = ebisBox.getEBGraph();

      IntVectSet irregIVS = ebisBox.getIrregIVS(curBox);
      IntVectSet multiIVS = ebisBox.getMultiCells(curBox);

      //cache the vofIterators
      m_alphaDiagWeight[datind].define(irregIVS,ebgraph, 1);
      m_betaDiagWeight [datind].define(irregIVS,ebgraph, 1);
      m_vofIterIrreg   [datind].define(irregIVS,ebgraph);
      m_vofIterMulti   [datind].define(multiIVS,ebgraph);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = irregIVS;
          IntVectSet hiIrreg = irregIVS;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofIterDomLo[idir][datind].define(loIrreg,ebgraph);
          m_vofIterDomHi[idir][datind].define(hiIrreg,ebgraph);
        }

      VoFIterator &  vofit = m_vofIterIrreg[datind];
      const Vector<VolIndex>& vofvec = vofit.getVector();
      // cast from VolIndex to BaseIndex
      Vector<RefCountedPtr<BaseIndex> >    dstVoF(vofvec.size());
      Vector<RefCountedPtr<BaseStencil> > stencil(vofvec.size());
      for(int ivec = 0; ivec < vofvec.size(); ivec++)
        {
          const VolIndex& vof = vofvec[ivec];
          VoFStencil vofsten;

          //bcoef is included here in the flux consistent
          //with the regular

          getDivFStencil(vofsten,vof, datind);
          if (fluxStencil != NULL)
            {
              BaseIVFAB<VoFStencil>& fluxStencilBaseIVFAB = (*fluxStencil)[datind];
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
                      weightedAreaFrac *= (*m_bcoef)[datind][idir](faces[i],0);
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
                      weightedAreaFrac *= (*m_bcoef)[datind][idir](faces[i],0);
                      faceAreaFrac +=  weightedAreaFrac;
                    }
                  adjust += -weightedAreaFrac /(m_dx*m_dx);
                }
              betaWeight += adjust;
            }
            
          dstVoF[ivec]  = RefCountedPtr<BaseIndex  >(new  VolIndex(vof));
          stencil[ivec] = RefCountedPtr<BaseStencil>(new VoFStencil(vofsten));

          //weight for in identity term
          Real volFrac = ebisBox.volFrac(vof);
          Real alphaWeight = (*m_acoef)[datind](vof, 0);
          alphaWeight *= volFrac;

          m_alphaDiagWeight[datind](vof, 0) = alphaWeight;
          m_betaDiagWeight[datind]( vof, 0)  = betaWeight;
        }

      m_opEBStencil[datind] = RefCountedPtr<VCAggStencil>
        (new VCAggStencil(dstVoF, stencil, phiProxy[datind], rhsProxy[datind], m_relCoef[datind], m_alphaDiagWeight[datind],1));

    }//dit

  calculateAlphaWeight();
  calculateRelaxationCoefficient();

  if (m_hasFine)
    {

      int ncomp = 1;
      m_fastFR.define(m_eblgFine, m_eblg, m_refToFine, ncomp, s_forceNoEBCF);
      m_hasEBCF = m_fastFR.hasEBCF();
    }
  defineEBCFStencils();

}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
defineEBCFStencils()
{
  ///this routine is ugly and complicated.
  //I will attempt to comment it but I fear it is a lost cause
  //because the algorithm is so arcane.
  //We are attempting to only do stuff at the very specific
  //points where there is an embedded boundary crossing a coarse-fine
  //interface.   We happen to know that EBFastFR already has done this
  //choreography and we want to leverage it.

  //EBFastFR has data structures in it that serve as buffers and so on
  //that we will (thankfully) be able to leave alone.   We are only
  //going to access the data structures wherein it identifies which
  //coarse cells are associated with the coarse-fine interface
  // where the EB crosses and use this list to build up  which faces
  // need to be cal

  //important factoid: beta gets multiplied in at the last minute
  //(on evaluation) because it can change as diffusion solvers progress.
  if (m_hasFine && m_hasEBCF)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              //coarse fine stuff is between me and next finer level
              //fine stuff lives over m_eblgfine
              //coar stuff lives over m_eblg
              int index = m_fastFR.index(idir, sit());
              m_stencilCoar[index].define(m_eblg.getDBL());
              m_faceitCoar [index].define(m_eblg.getDBL());

              DataIterator dit = m_eblg.getDBL().dataIterator(); 
              int nbox = dit.size();

#pragma omp parallel for
              for(int mybox=0; mybox<nbox; mybox++)
                {
                  Vector<FaceIndex>& facesEBCFCoar =  m_faceitCoar[index][dit[mybox]];
                  Vector<VoFStencil>& stencEBCFCoar= m_stencilCoar[index][dit[mybox]];
                  Vector<VoFIterator>& vofitlist = m_fastFR.getVoFItCoar(dit[mybox], idir, sit());
                  //first build up the list of the faces
                  for (int ivofit = 0; ivofit < vofitlist.size(); ivofit++)
                    {
                      VoFIterator& vofit = vofitlist[ivofit];
                      for (vofit.reset(); vofit.ok(); ++vofit)
                        {
                          //on the coarse side of the CF interface we are
                          //looking in the flip direction because we look
                          //back at the interface
                          Vector<FaceIndex> facespt = m_eblg.getEBISL()[dit[mybox]].getFaces(vofit(), idir, flip(sit()));
                          facesEBCFCoar.append(facespt);
                        }
                    }

                  stencEBCFCoar.resize(facesEBCFCoar.size());
                  for (int iface = 0; iface < stencEBCFCoar.size(); iface++)
                    {
                      IntVectSet cfivs; //does not apply here
                      getFluxStencil(stencEBCFCoar[iface], facesEBCFCoar[iface], dit[mybox]);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousPhysBC)
{
  CH_TIME("NWOEBConductivityOp::residual");
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_phi.ghostVect() == m_ghostPhi);
  applyOp(a_residual,a_phi,NULL, a_homogeneousPhysBC, true);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
preCond(LevelData<EBCellFAB>&       a_lhs,
        const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("NWOEBConductivityOp::preCond");
  EBLevelDataOps::assign(a_lhs, a_rhs);
  EBLevelDataOps::scale(a_lhs, m_relCoef);

  relax(a_lhs, a_rhs, 40);
}

//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
applyOp(LevelData<EBCellFAB>&                    a_lhs,
        const LevelData<EBCellFAB>&              a_phi,
        const LevelData<EBCellFAB>* const        a_phiCoar,
        const bool&                              a_homogeneousPhysBC,
        const bool&                              a_homogeneousCFBC)
{
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)a_phi;
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
applyOp(LevelData<EBCellFAB>&                    a_lhs,
        const LevelData<EBCellFAB>&              a_phi,
        bool                                     a_homogeneous)
{
  CH_TIME("ebco::applyOp");
  LevelData<EBCellFAB>&  phi = (LevelData<EBCellFAB>&) a_phi;
  {
    CH_TIME("ghostcell fills");
    fillPhiGhost(a_phi, a_homogeneous);
    phi.exchange(m_exchangeCopier);
  }

  /**
     This is not done anymore because alpha and a are now part of the stencil
     and part of regular apply
     EBLevelDataOps::setToZero(a_lhs);
     incr( a_lhs, a_phi, m_alpha);
  **/
  {
    CH_TIME("applying op without bcs");
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
fillPhiGhost(const LevelData<EBCellFAB>& a_phi, bool a_homog) const
{
  CH_TIME("nwoebco::fillghostLD");
  DataIterator dit = m_eblg.getDBL().dataIterator();
  int nbox=dit.size();
#pragma omp parallel for
  for (int mybox=0;mybox<nbox; mybox++)
    {
      fillPhiGhost(a_phi[dit[mybox]], dit[mybox], a_homog);
    }
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)(a_phi);
  phi.exchange(m_exchangeCopier);
}
void
NWOEBConductivityOp::
fillPhiGhost(const EBCellFAB& a_phi, const DataIndex& a_datInd, bool a_homog) const
{
  CH_TIME("nwoebco::fillghostfab");

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
  CH_TIME("nwoebco::applyopRegular");
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
  CH_TIME("nwoebco::applyOpIrregular");
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
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("ebco::create");
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                const LevelData<EBCellFAB>& a_rhs,
                const int &                 a_refRat)
{
  CH_TIME("ebco::createCoar");
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  CH_assert(m_eblg.getDBL().coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, m_eblg.getDBL(), a_refRat);

  EBISLayout ebislCoarsenedFine;
  IntVect ghostVec = a_rhs.ghostVect();

  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), a_refRat);
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
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("ebco::assign");
  EBLevelDataOps::assign(a_lhs,a_rhs);
}
//-----------------------------------------------------------------------
Real
NWOEBConductivityOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  CH_TIME("ebco::dotProd");
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  CH_TIME("ebco::incr");
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  CH_TIME("ebco::axby");
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  CH_TIME("ebco::scale");
  EBLevelDataOps::scale(a_lhs,a_scale);
}
//-------------------------------
Real 
NWOEBConductivityOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
  CH_TIMERS("NWOEBConductivityOp::norm");
  CH_TIMER("mpi_allreduce",t1);

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
  //  EBLevelDataOps::gatherBroadCast(maxNorm, volume, 0);
  CH_STOP(t1);

  return maxNorm;
}

Real 
NWOEBConductivityOp::
localMaxNorm(const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBAMRPoissonOp::localMaxNorm");
  return  EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);
  //  ProblemDomain domain;
  //  Real volume;
  //  return EBLevelDataOps::kappaNorm(volume,a_rhs,EBLEVELDATAOPS_ALLVOFS,domain,0);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  CH_TIME("ebco::setToZero");
  EBLevelDataOps::setToZero(a_lhs);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("ebco::setVal");
  EBLevelDataOps::setVal(a_lhs, a_value);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
createCoarser(LevelData<EBCellFAB>&       a_coar,
              const LevelData<EBCellFAB>& a_fine,
              bool                        a_ghosted)
{
  CH_TIME("ebco::createCoarser");
  CH_assert(a_fine.nComp() == 1);
  const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), 2);

  int nghost = a_fine.ghostVect()[0];
  EBISLayout coarEBISL;

  const EBIndexSpace* const ebisPtr = m_eblg.getEBIS();
  ebisPtr->fillEBISLayout(coarEBISL,
                          dbl, coarDom, nghost);

  EBCellFactory ebcellfact(coarEBISL);
  a_coar.define(dbl, 1,a_fine.ghostVect(),ebcellfact);
}
//------------------------------------------------------------------------
/*****/
void
NWOEBConductivityOp::
relax(LevelData<EBCellFAB>&       a_phi,
      const LevelData<EBCellFAB>& a_rhs,
      int                         a_iterations)
{
  CH_TIME("nwoebco::relax");
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  //change to false to get a lot of printouts
  static bool printedStuff = true;

  if(!printedStuff)
    {
      pout() << "rhs: " << endl;;
      LevelData<EBCellFAB>* rhs = (LevelData<EBCellFAB>*)(&a_rhs);
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
              CH_TIME("ghostfill and homogcfinterp");
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
      LevelData<EBCellFAB>* rhs = (LevelData<EBCellFAB>*)(&a_rhs);
      printMaxMinLDCell(rhs);
      pout() << "phi coming out of relax: "  << endl;
      printMaxMinLDCell(&a_phi);
    }
  printedStuff = true;
}

void
NWOEBConductivityOp::
homogeneousCFInterp(LevelData<EBCellFAB>&   a_phif)
{
  CH_TIME("nwoebco::homog_cfinterp");
  if (m_hasCoar)
    {
      m_interpWithCoarser->coarseFineInterpH(a_phif,  0, 0, 1);
    }
}

//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
gsrbColor(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_rhs,
          const IntVect&              a_color)
{

  CH_TIME("nwoebco::gsrbColor");
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
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_rhsThisLevel)
{
  CH_TIME("NWOEBConductivityOp::restrictResidual");

  CH_assert(a_resCoar.nComp() == 1);
  CH_assert(a_phiThisLevel.nComp() == 1);
  CH_assert(a_rhsThisLevel.nComp() == 1);

  LevelData<EBCellFAB> resThisLevel;
  bool homogeneous = true;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_rhsThisLevel.ghostVect();

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
prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                 const LevelData<EBCellFAB>& a_correctCoar)
{
  CH_TIME("NWOEBConductivityOp::prolongIncrement");
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
AMRResidual(LevelData<EBCellFAB>&       a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            const LevelData<EBCellFAB>& a_rhs,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("NWOEBConductivityOp::AMRResidual");
  CH_TIMER("AMROperator", t1);
  CH_TIMER("axby", t2);
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostRHS);
  CH_assert(a_residual.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

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
AMROperator(LevelData<EBCellFAB>&       a_LofPhi,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoar,
            bool a_homogeneousPhysBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("NWOEBConductivityOp::AMROperator");
  CH_TIMER("applyOp", t1);
  CH_TIMER("reflux", t2);
  CH_assert(a_LofPhi.ghostVect() == m_ghostRHS);
  CH_assert(a_LofPhi.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

  //apply the operator between this and the next coarser level.
  CH_START(t1);
  applyOp(a_LofPhi, a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);
  CH_STOP(t1);

  //now reflux to enforce flux-matching from finer levels
  if (m_hasFine)
    {
      CH_assert(a_finerOp != NULL);
      CH_START(t2);

      reflux(a_LofPhi, a_phiFine, a_phi, a_finerOp);

      CH_STOP(t2);
    }
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
reflux(LevelData<EBCellFAB>& a_residual,
       const LevelData<EBCellFAB>& a_phiFine,
       const LevelData<EBCellFAB>& a_phi,
       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIMERS("NWOEBConductivityOp::fastReflux");
  CH_TIMER("setToZero",t2);
  CH_TIMER("incrementCoar",t3);
  CH_TIMER("incrementFine",t4);
  CH_TIMER("reflux_from_reg",t5);
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
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("NWOEBConductivityOp::incrementFRCoar");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);

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
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex>*  faceit;
              Vector<VoFStencil>* stencil;
              int index = EBFastFR::index(idir, sit());
              if (m_hasEBCF)
                {
                  faceit  = &( m_faceitCoar[index][dit[mybox]]);
                  stencil = &(m_stencilCoar[index][dit[mybox]]);
                }
              getFluxEBCF(coarflux, coarfab, ghostedBox, *faceit, *stencil);
            }

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
getFluxEBCF(EBFaceFAB&                    a_flux,
            const EBCellFAB&              a_phi,
            const Box&                    a_ghostedBox,
            Vector<FaceIndex>&            a_faceitEBCF,
            Vector<VoFStencil>&           a_stenEBCF)
{
  CH_TIME("NWOEBConductivityOp::getFluxEBCF");

  //only do the evil stuff if you have a coarse-fine /  EB crossing situation

  if (m_hasEBCF)
    {
      CH_TIME("EBCF stuff");
      for (int iface = 0; iface < a_faceitEBCF.size(); iface++)
        {
          const FaceIndex& face =     a_faceitEBCF[iface];
          const VoFStencil& stencil   = a_stenEBCF[iface];
          Real fluxval = 0;
          for (int isten = 0; isten < stencil.size(); isten++)
            {
              fluxval += stencil.weight(isten)*(a_phi(stencil.vof(isten), 0));
            }
          //note the last minute beta
          a_flux(face, 0) = m_beta*fluxval;
        }
    }
}

//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
getFlux(EBFluxFAB&                    a_flux,
        const LevelData<EBCellFAB>&   a_data,
        const Box&                    a_grid,
        const DataIndex&              a_dit,
        Real                          a_scale)
{
  CH_TIME("ebco::getflux1");
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
        const ProblemDomain&          a_domain,
        const EBISBox&                a_ebisBox,
        const Real&                   a_dx,
        const DataIndex&              a_datInd,
        const int&                    a_idir)
{
  CH_TIME("ebco::getFlux2");
  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());
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
  CH_TIME("ebco::getFluxRegOnly");
  const ProblemDomain& domain = m_eblg.getDomain();

  //has some extra cells so...
  a_fluxCentroid.setVal(0.);
  int ncomp = a_phi.nComp();
  CH_assert(ncomp == a_fluxCentroid.nComp());
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
                const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi,
                AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("NWOEBConductivityOp::incrementFRFine");
  CH_assert(a_phiFine.nComp() == 1);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(m_hasFine);
  int ncomp = 1;
  Interval interv(0,0);
  NWOEBConductivityOp& finerEBAMROp = (NWOEBConductivityOp& )(*a_finerOp);

  //ghost cells of phiFine need to be filled
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;
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
  CH_TIME("ebco::getflux3");
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
AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //dummy. there is no coarse when this is called
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostRHS);
  LevelData<EBCellFAB> phiC;
  AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousPhysBC, a_finerOp);
}
//-----------------------------------------------------------------------

void
NWOEBConductivityOp::
AMRResidualNF(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousPhysBC)
{
  CH_TIME("ebco::amrresNF");
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostRHS);

  AMROperatorNF(a_residual, a_phi, a_phiCoar,
                a_homogeneousPhysBC);
  axby(a_residual,a_residual,a_rhs,-1.0, 1.0);
}
//-----------------------------------------------------------------------

void
NWOEBConductivityOp::
AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              bool a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  CH_TIME("ebco::amrOpNC");
  //dummy. there is no coarse when this is called
  CH_assert(a_LofPhi.ghostVect() == m_ghostRHS);
  LevelData<EBCellFAB> phiC;
  AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
              a_homogeneousPhysBC, a_finerOp);
}
//-----------------------------------------------------------------------

void
NWOEBConductivityOp::
AMROperatorNF(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoar,
              bool a_homogeneousPhysBC)
{
  CH_assert(a_LofPhi.ghostVect() == m_ghostRHS);

  applyOp(a_LofPhi,a_phi, &a_phiCoar,  a_homogeneousPhysBC, false);

}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
AMRRestrict(LevelData<EBCellFAB>&       a_resCoar,
            const LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_correction,
            const LevelData<EBCellFAB>& a_coarCorrection, 
            bool a_skip_res )
{
  CH_TIME("NWOEBConductivityOp::AMRRestrict");
  CH_assert(a_residual.ghostVect() == m_ghostRHS);
  CH_assert(a_correction.ghostVect() == m_ghostPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostPhi);
  CH_assert(!a_skip_res);

  CH_assert(a_residual.nComp() == 1);
  CH_assert(a_resCoar.nComp() == 1);
  CH_assert(a_correction.nComp() == 1);

  LevelData<EBCellFAB> resThisLevel;
  bool homogeneousPhys = true;
  bool homogeneousCF =   false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  resThisLevel.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);
  EBLevelDataOps::setVal(resThisLevel, 0.0);

  //API says that we must average(a_residual - L(correction, coarCorrection))
  applyOp(resThisLevel, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);
  incr(resThisLevel, a_residual, -1.0);
  scale(resThisLevel,-1.0);

  //use our nifty averaging operator
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebAverage.average(a_resCoar, resThisLevel, variables);

}
//-----------------------------------------------------------------------
Real
NWOEBConductivityOp::
AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
        const LevelData<EBCellFAB>& a_fineResid,
        const int& a_refRat,
        const int& a_ord)

{
  // compute norm over all cells on coarse not covered by finer
  CH_TIME("NWOEBConductivityOp::AMRNorm");
  MayDay::Error("never called");
  //return norm of temp
  return norm(a_coarResid, a_ord);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
AMRProlong(LevelData<EBCellFAB>&       a_correction,
           const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("NWOEBConductivityOp::AMRProlong");
  //use cached interpolation object
  Interval variables(0, 0);
  CH_assert(m_hasInterpAve);
  m_ebInterp.pwcInterp(a_correction, a_coarCorrection, variables);
}
//-----------------------------------------------------------------------
void
NWOEBConductivityOp::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarCorrection)
{
  CH_TIME("NWOEBConductivityOp::AMRUpdateResidual");
  CH_assert(a_residual.ghostVect()   == m_ghostRHS);
  CH_assert(a_correction.ghostVect() == m_ghostPhi);
  CH_assert(a_coarCorrection.ghostVect() == m_ghostPhi);

  LevelData<EBCellFAB> lcorr;
  bool homogeneousPhys = true;
  bool homogeneousCF   = false;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  lcorr.define(m_eblg.getDBL(), 1, ghostVec, ebcellfactTL);

  applyOp(lcorr, a_correction, &a_coarCorrection, homogeneousPhys, homogeneousCF);

  incr(a_residual, lcorr, -1);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
