#include "AMReX_EBCFInterp.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBLevelDataOps.H"

namespace amrex
{
  void null_deleter_nwo_sten(BaseStencil * a_sten)
 {}
  void null_deleter_nwo_ind(BaseIndex* a_sten)
 {}
  /***********************/
  void
  EBCFInterp::
  define(const EBLevelGrid &       a_eblgFine,
         const EBLevelGrid &       a_eblgCoar,
         const int         &       a_refRat,
         const int         &       a_ghostCellsInData,
         const int         &       a_ghostCellsToFill,
         int  a_orderOfPolynomial,
         bool a_slowMode)
  {
    BL_ASSERT(a_refRat > 0);
    BL_ASSERT(a_ghostCellsInData >= 0);
    BL_ASSERT(a_ghostCellsToFill >= 0);
    BL_ASSERT(a_ghostCellsInData >= a_ghostCellsToFill);
    BL_ASSERT(a_eblgFine.coarsenable(a_refRat));
    m_orderOfPolynomial = a_orderOfPolynomial;
    m_slowMode = a_slowMode;
    m_isDefined = true;
    m_eblgFine  = a_eblgFine;
    m_eblgCoar  = a_eblgCoar;
    m_refRat    = a_refRat;
    m_dataGhost = a_ghostCellsInData;
    m_fillGhost = a_ghostCellsToFill;
    
    //need this many for the stencil
    m_buffGhost = 2;
    defineInternals();
  }
  /***********************/
  void 
  EBCFInterp::
  defineInternals()
  {
    BL_PROFILE("NWOEBCFI::defineInternals");
    coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
    m_eblgCoFi.setMaxRefinementRatio(m_refRat);

    EBCellFactory factFine(m_eblgFine.getEBISL());
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    //variable number does not matter here.
    int nvar = 1;
    FabArray<EBCellFAB> proxyLevel(m_eblgFine.getDBL(),m_eblgFine.getDM(), nvar, m_dataGhost, MFInfo(), factFine);
    FabArray<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_buffGhost, MFInfo(), factCoFi);

    if(m_slowMode)
    {
      m_slowStencils.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
      m_slowVoFs.define(    m_eblgFine.getDBL(), m_eblgFine.getDM());
      
    }
    else
    {
      m_stencil.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    }
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      Vector< std::shared_ptr<BaseIndex  > > baseDstVoFs;
      Vector< std::shared_ptr<BaseStencil> > baseSten;
      IntVectSet cfivs = getCFIVS(mfi);
      
      const EBISBox  & ebisFine =   m_eblgFine.getEBISL()[ mfi];
      const EBISBox  & ebisCoFi =   m_eblgCoFi.getEBISL()[ mfi];

      VoFIterator vofit(cfivs, ebisFine.getEBGraph());
      const Vector<VolIndex>& volvec = vofit.getVector();
      baseDstVoFs.resize(volvec.size());
      baseSten.resize(   volvec.size());
      Vector<VoFStencil> allsten(volvec.size());
      for(int ivec = 0; ivec < volvec.size(); ivec++)
      {
        getStencil(allsten[ivec],  volvec[ivec], ebisFine, ebisCoFi);
        baseSten    [ivec]  = std::shared_ptr<BaseStencil>(            &allsten[ivec] , &null_deleter_nwo_sten);
        baseDstVoFs [ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&volvec[ivec]), &null_deleter_nwo_ind);
      }

      EBCellFAB& coarProxy =   bufferCoFi[mfi];
      EBCellFAB& fineProxy =   proxyLevel[mfi];
      if(m_slowMode)
      {
        m_slowStencils[mfi] = allsten;
        m_slowVoFs    [mfi] = volvec;
      }
      else
      {
        m_stencil[mfi] = std::shared_ptr<AggStencil <EBCellFAB, EBCellFAB>  >
          (new AggStencil<EBCellFAB, EBCellFAB >(baseDstVoFs, baseSten, coarProxy, fineProxy));
      }
    }
  }
  /***********************/
  IntVectSet
  EBCFInterp::
  getCFIVS(const MFIter& a_mfi)
  {
    IntVectSet cfivs;
    if(m_fillGhost == 1)
    {
      cfivs   = (*m_eblgFine.getCFIVS())[a_mfi];
    }
    else
    {
      Box grownBox = m_eblgFine.getDBL()[a_mfi];
      grownBox.grow(m_fillGhost);
      grownBox &= m_eblgFine.getDomain();
      cfivs = IntVectSet(grownBox);
      for(int ibox = 0; ibox < m_eblgFine.getDBL().size(); ibox++)
      {
        cfivs -= m_eblgFine.getDBL()[ibox];
      }
    }
    return cfivs;
  }
  /***********************/
  void
  EBCFInterp::
  getStencil(VoFStencil           & a_stencil,
             const VolIndex       & a_vofFine,
             const EBISBox        & a_ebisFine,
             const EBISBox        & a_ebisCoFi)
  {
    VolIndex fineVoF = a_vofFine;
    //the values of these do not matter as this is interpolation
    Real dxFine = 1.0;  Real dxCoar = m_refRat;
    a_stencil.clear();
    VolIndex coarVoF = a_ebisFine.coarsen(a_vofFine);
    RealVect coarLoc = EBArith::getVoFLocation(coarVoF, dxCoar, RealVect::Zero);
    RealVect fineLoc = EBArith::getVoFLocation(fineVoF, dxFine, RealVect::Zero);
    RealVect dist = fineLoc - coarLoc;
    EBArith::getExtrapolationStencil(a_stencil, dist, dxCoar*RealVect::Unit, coarVoF, a_ebisCoFi, m_orderOfPolynomial);
  }  
  /***********************/
  void 
  EBCFInterp::
  coarseFineInterp(FabArray<EBCellFAB>&       a_phif,
                   const FabArray<EBCellFAB>& a_phic,
                   int isrc, int idst, int inco)
  {
    BL_PROFILE("NWOEBCFI::coarseFineInterp");

    int nvar = idst + inco;
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    FabArray<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_buffGhost, MFInfo(), factCoFi);
    //need to copy ghost cell data 
    bufferCoFi.copy(a_phic, isrc, idst, inco, 0, m_buffGhost);
    int ibox = 0;
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {

      //breaking things down for debugging.
      EBCellFAB& phifab =     a_phif[mfi];
      EBCellFAB& buffab = bufferCoFi[mfi];
      if(m_slowMode)
      {

        Vector<VolIndex  >& vofs     = m_slowVoFs[mfi];
        Vector<VoFStencil>& stencils = m_slowStencils[mfi];
        for(int ivof = 0; ivof < vofs.size(); ivof++)
        {
          for(int icomp = 0; icomp < inco; icomp++)
          {
            int srccomp = isrc + icomp;
            int dstcomp = idst + icomp;

            VoFStencil vofsten = stencils[ivof];
            vofsten.setAllVariables(srccomp);
            Real value = applyVoFStencil(stencils[ivof], buffab);
            phifab(vofs[ivof], dstcomp)= value;
          }
        }
      }
      else
      {
        //false is for increment only
        m_stencil[mfi]->apply(phifab,buffab, isrc, idst, inco, false);
      }
      ibox++;
    }
  }

  /***********************/
  void 
  EBCFInterp::
  coarseFineInterpH(FabArray<EBCellFAB> & a_phif,
                    int isrc, int idst, int inco)
  {
    BL_PROFILE("NWOEBCFI::coarseFineInterpH");
    BL_ASSERT(!m_slowMode);
    int nvar = idst + inco;
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    FabArray<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_buffGhost, MFInfo(), factCoFi);
    EBLevelDataOps::setVal(bufferCoFi, 0.0);
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {

      //false is for increment only
      m_stencil[mfi]->apply(a_phif[mfi],
                            bufferCoFi[mfi],
                            isrc, idst, inco, false);
    }

  }
}

