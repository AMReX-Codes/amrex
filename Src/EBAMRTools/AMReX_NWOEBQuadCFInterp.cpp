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

#include "AMReX_NWOEBQuadCFInterp.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBLevelDataOps.H"

namespace amrex
{
  void null_deleter_nwo_vsten(BaseStencil * a_sten)
 {}
  void null_deleter_nwo_vof(BaseIndex* a_sten)
 {}
  /***********************/
  void
  NWOEBQuadCFInterp::
  define(const EBLevelGrid &       a_eblgFine,
         const EBLevelGrid &       a_eblgCoar,
         const int         &       a_refRat,
         const int         &       a_ghost)
  {
    BL_ASSERT(a_refRat > 0);
    BL_ASSERT(a_ghost >= 0);
    BL_ASSERT(a_eblgFine.coarsenable(a_refRat));
    m_isDefined = true;
    m_eblgFine  = a_eblgFine;
    m_eblgCoar  = a_eblgCoar;
    m_refRat    = a_refRat;
    m_ghost     = a_ghost;
    defineInternals();
  }
  /***********************/
  void 
  NWOEBQuadCFInterp::
  defineInternals()
  {
    BL_PROFILE("NWOEBCFI::defineInternals");
    coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
    m_eblgCoFi.setMaxRefinementRatio(m_refRat);

    EBCellFactory factFine(m_eblgFine.getEBISL());
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    //variable number does not matter here.
    int nvar = 1;
    FabArray<EBCellFAB> proxyLevel(m_eblgFine.getDBL(),m_eblgFine.getDM(), nvar, m_ghost, MFInfo(), factFine);
    m_bufferCoFi.define(           m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_ghost, MFInfo(), factCoFi);

    m_stencil.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      std::vector< std::shared_ptr<BaseIndex  > > baseDstVoFs;
      std::vector< std::shared_ptr<BaseStencil> > baseSten;

      const IntVectSet& cfivs   = (*m_eblgFine.getCFIVS())[mfi];
      const EBISBox  & ebisFine =   m_eblgFine.getEBISL()[ mfi];
      const EBISBox  & ebisCoFi =   m_eblgCoFi.getEBISL()[ mfi];

      VoFItrator vofit(cfivs, ebisFine.getEBGraph());
      const std::vector<VolIndex>& volvec = vofit.getVector();
      a_baseDstVoFs.resize(volvec.size());
      a_stencils.resize(   volvec.size());
      std::vector<VoFStencil> allsten(volvec.size());
      for(int ivec = 0; ivec < volvec.size(); ivec++)
      {
        getStencil(allsten[ivec],  volvec[ivec], ebisFine, ebisCoFi);
        stencils    [ivec]  = std::shared_ptr<BaseStencil>(            &allsten[ivec] , &null_deleter_nwo_sten);
        baseDstVoFs [ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&volvec[ivec]), &null_deleter_nwo_ind);
      }

      EBCellFAB& coarProxy = m_bufferCoFi[mfi];
      EBCellFAB& fineProxy =   proxyLevel[mfi];
      m_stencil[mfi] = std::shared_ptr<AggStencil <EBCellFAB, EBCellFAB>  >
        (new AggStencil<EBCellFAB, EBCellFAB >(baseDstVoFs, baseSten, coarProxy, fineProxy));
    }
  }
  /***********************/
  void
  NWOEBQuadCFInterp::
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
    EBArith::getExtrapolationStencil(a_stencil, dist, dxCoar*RealVect::Unit, coarVoF, a_ebisCoFi);
  }  
  /***********************/
  void 
  NWOEBQuadCFInterp::
  coarseFineInterp(FabArray<EBCellFAB>&       a_phif,
                   const FabArray<EBCellFAB>& a_phic,
                   int isrc, int idst, int inco)
  {
    BL_PROFILE("NWOEBCFI::coarseFineInterp");

    //do not need to copy ghost cell data 
    m_bufferCoFi.copy(a_phic, isrc, idst, inco, 0, 0);
    //but we might NEED ghost cell data (just cannot assume it was correct in input.
    m_bufferCoFi.FillBoundary();
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {

      //false is for increment only
      m_stencil[mfi]->apply(a_phif[mfi],
                            m_bufferCoFi[mfi],
                            isrc, idst, inco, false);
    }
  }

  /***********************/
  void 
  NWOEBQuadCFInterp::
  coarseFineInterpH(FabArray<EBCellFAB> & a_phif,
                    int isrc, int idst, int inco)
  {
    BL_PROFILE("NWOEBCFI::coarseFineInterpH");
    EBLevelDataOps::setVal(m_bufferCoFi, 0.0);
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {

      //false is for increment only
      m_stencil[mfi]->apply(a_phif[mfi],
                            m_bufferCoFi[mfi],
                            isrc, idst, inco, false);
    }

  }
}

