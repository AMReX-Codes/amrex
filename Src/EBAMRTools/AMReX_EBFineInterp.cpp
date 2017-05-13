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


#include "AMReX_EBFineInterp.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBLoHiCenter.H"
#include "AMReX_EBFortND_F.H"


namespace amrex
{
  void null_deleter_ebfi_sten(BaseStencil * a_sten)
  {}
  void null_deleter_ebfi_ind(BaseIndex* a_sten)
  {}
  /************************************/
  void
  EBFineInterp::
  define(const EBLevelGrid   & a_eblgFine,
         const EBLevelGrid   & a_eblgCoar,
         const int           & a_nref,
         const int           & a_ghostCellsInData,
         int a_orderOfPolynomial)
  {
    BL_ASSERT(a_nref > 0);
    BL_ASSERT(a_eblgFine.coarsenable(a_nref));
    
    m_isDefined = true;
    m_eblgFine          = a_eblgFine;
    m_eblgCoar          = a_eblgCoar;
    m_refRat            = a_nref;
    m_orderOfPolynomial = a_orderOfPolynomial;
    m_buffGhost = 2;
    m_dataGhost = a_ghostCellsInData;

    defineInternals();
  }
  /************************************/
  void
  EBFineInterp::
  defineInternals()
  {
    BL_PROFILE("NWOEBCFI::defineInternals");
    coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
    m_eblgFine.setMaxCoarseningRatio(m_refRat);

    EBCellFactory factFine(m_eblgFine.getEBISL());
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    //variable number does not matter here.
    int nvar = 1;
    FabArray<EBCellFAB> proxyLevel(m_eblgFine.getDBL(),m_eblgFine.getDM(), nvar, m_dataGhost, MFInfo(), factFine);
    FabArray<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_buffGhost, MFInfo(), factCoFi);

    m_stencil.define(m_eblgFine.getDBL(), m_eblgFine.getDM());

    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      std::vector< std::shared_ptr<BaseIndex  > > baseDstVoFs;
      std::vector< std::shared_ptr<BaseStencil> > baseSten;
      const EBISBox  & ebisFine =   m_eblgFine.getEBISL()[ mfi];
      const EBISBox  & ebisCoFi =   m_eblgCoFi.getEBISL()[ mfi];

      Box gridFine =m_eblgFine.getDBL()[mfi];
      IntVectSet ivsIrreg = ebisFine.getIrregIVS(gridFine);
      
      VoFIterator vofit(ivsIrreg, ebisFine.getEBGraph());
      const std::vector<VolIndex>& volvec = vofit.getVector();
      baseDstVoFs.resize(volvec.size());
      baseSten.resize(   volvec.size());
      std::vector<VoFStencil> allsten(volvec.size());
      for(int ivec = 0; ivec < volvec.size(); ivec++)
      {
        getStencil(allsten[ivec],  volvec[ivec], ebisFine, ebisCoFi, mfi);
        baseSten    [ivec]  = std::shared_ptr<BaseStencil>(            &allsten[ivec] , &null_deleter_ebfi_sten);
        baseDstVoFs [ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&volvec[ivec]), &null_deleter_ebfi_ind);
      }

      EBCellFAB& coarProxy =   bufferCoFi[mfi];
      EBCellFAB& fineProxy =   proxyLevel[mfi];

      m_stencil[mfi] = std::shared_ptr<AggStencil <EBCellFAB, EBCellFAB>  >
        (new AggStencil<EBCellFAB, EBCellFAB >(baseDstVoFs, baseSten, coarProxy, fineProxy));

    }
  }
  /***********************/
  void
  EBFineInterp::
  getStencil(VoFStencil           & a_stencil,
             const VolIndex       & a_vofFine,
             const EBISBox        & a_ebisFine,
             const EBISBox        & a_ebisCoFi,
             const MFIter         & a_mfi)
  {
    VolIndex fineVoF = a_vofFine;
    //the values of these do not matter as this is interpolation
    Real dxFine = 1.0;  Real dxCoar = m_refRat;
    a_stencil.clear();
    VolIndex coarVoF = m_eblgFine.getEBISL().coarsen(fineVoF, m_refRat, a_mfi);
    RealVect coarLoc = EBArith::getVoFLocation(coarVoF, dxCoar, RealVect::Zero);
    RealVect fineLoc = EBArith::getVoFLocation(fineVoF, dxFine, RealVect::Zero);
    RealVect dist = fineLoc - coarLoc;
    EBArith::getExtrapolationStencil(a_stencil, dist, dxCoar*RealVect::Unit, coarVoF, a_ebisCoFi, m_orderOfPolynomial);
  }  
  /************************************/
  void
  EBFineInterp::interpolate(FabArray<EBCellFAB>      & a_dataFine,
                            const FabArray<EBCellFAB>& a_dataCoar,
                            int isrc, int idst, int inco)
  {
    BL_ASSERT(isDefined());

    int nvar = isrc + inco;
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    FabArray<EBCellFAB> dataCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_buffGhost, MFInfo(), factCoFi);

    //need to copy into ghost cell data 
    dataCoFi.copy(a_dataCoar, isrc, idst, inco, 0, m_buffGhost);

    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      BaseFab<Real>      & regFine = a_dataFine[mfi].getSingleValuedFAB();
      const BaseFab<Real>& regCoar =   dataCoFi[mfi].getSingleValuedFAB();
      const Box& gridFine = m_eblgFine.getDBL()[mfi];
      if(m_orderOfPolynomial < 1)
      {
        ebfnd_pwcinterp(BL_TO_FORTRAN_FAB(regFine),
                        BL_TO_FORTRAN_FAB(regCoar),
                        BL_TO_FORTRAN_BOX(gridFine),
                        &m_refRat, &isrc, &idst, &inco);
      }
      else if(m_orderOfPolynomial < 2)
      {
#if 0
        Box loBox [BL_SPACEDIM];
        Box hiBox [BL_SPACEDIM];
        Box ceBox [BL_SPACEDIM];
        int hasLo [BL_SPACEDIM];
        int hasHi [BL_SPACEDIM];
        int nearAnyBoundary;
        Box inBoxCoar = m_eblgCoFi.getDBL()[mfi];
        Box domainCoar= m_eblgCoFi.getDomain();
        EBLoHiCenAllDirs(loBox ,
                         hiBox ,
                         ceBox,
                         hasLo ,
                         hasHi ,
                         nearAnyBoundary,
                         inBoxCoar,
                         domainCoar);
        if(nearAnyBoundary == 0)
        {
          ebfnd_pwlinterp_nobound(BL_TO_FORTRAN_FAB(regFine),
                                  BL_TO_FORTRAN_FAB(regCoar),
                                  BL_TO_FORTRAN_BOX(gridFine),
                                  &m_refRat, &isrc, &idst, &inco);
        }
        else
        {

          // first add in constant bit
          ebfnd_pwcinterp(BL_TO_FORTRAN_FAB(regFine),
                          BL_TO_FORTRAN_FAB(regCoar),
                          BL_TO_FORTRAN_BOX(gridFine),
                          &m_refRat, &isrc, &idst, &inco);
          Box refBox(IntVect::Zero, IntVect::Zero);
          refBox.refine(m_refRat);
          //now increment each direction
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            ebfnd_pwl_incr_at_bound(BL_TO_FORTRAN_FAB(regFine),
                                    BL_TO_FORTRAN_FAB(regCoar),
                                    &hasLo[idir], &hasHi[idir], &idir,
                                    BL_TO_FORTRAN_BOX(loBox[idir]),
                                    BL_TO_FORTRAN_BOX(hiBox[idir]),
                                    BL_TO_FORTRAN_BOX(ceBox[idir]),
                                    BL_TO_FORTRAN_BOX(refBox),
                                    &m_refRat, &isrc, &idst, &inco);
          }
        }
#endif
      }
      else
      {
        //because mixed derivatives suck, only doing first derivs at domain boundaries.
        amrex::Abort("case not implemented");
      }

      //false is for increment only
      m_stencil[mfi]->apply(a_dataFine[mfi],
                              dataCoFi[mfi],
                            isrc, idst, inco, false);
    }

  }
  /************************************/
}

