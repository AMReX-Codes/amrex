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


#include "AMReX_EBCoarseAverage.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBLevelDataOps.H"
#include "AMReX_MultiFabUtil.H"


namespace amrex
{
  void null_deleter_ebc_sten(BaseStencil * a_sten)
 {}
  void null_deleter_ebc_ind(BaseIndex* a_sten)
 {}

  /************************************/
  EBCoarseAverage::
  EBCoarseAverage(const EBLevelGrid & a_eblgFine,
                  const EBLevelGrid & a_eblgCoar,
                  const int         & a_nref,
                  const int         & a_ghost,
                  const bool        & a_useKappaWeightingInStencil)
  {
    define(a_eblgFine, a_eblgCoar, a_nref, a_ghost, a_useKappaWeightingInStencil);
  }
  /************************************/
  void
  EBCoarseAverage::
  define(const EBLevelGrid & a_eblgFine,
         const EBLevelGrid & a_eblgCoar,
         const int         & a_nref,
         const int         & a_ghost,
         const bool        & a_useKappaWeightingInStencil)
  {
    BL_PROFILE("EBCoarseAverage::define(EBLG)");
    BL_ASSERT(a_nref > 0);
    BL_ASSERT(a_eblgFine.getEBISL().getGhost() >= 2);
    BL_ASSERT(a_eblgCoar.getEBISL().getGhost() >= 2);
    BL_ASSERT(a_eblgFine.coarsenable(a_nref));

    m_isDefined = true;
    m_refRat    = a_nref;
    m_eblgCoar  = a_eblgCoar;
    m_eblgFine  = a_eblgFine;
    m_ghost    = a_ghost;
    m_useKappaWeightingInStencil = a_useKappaWeightingInStencil;
    coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
    m_eblgCoFi.setMaxRefinementRatio(a_nref);
    defineStencils();
  }
  /************************************/
  void
  EBCoarseAverage::
  defineStencils()
  {
    m_stencil.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    EBCellFactory factFine(m_eblgFine.getEBISL());
    EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
    int nvar = 1; //not important
    FabArray<EBCellFAB> dummyFine(m_eblgFine.getDBL(),m_eblgFine.getDM(), nvar, m_ghost, MFInfo(), factCoFi);
    FabArray<EBCellFAB> dummyCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_ghost, MFInfo(), factCoFi);
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      BL_PROFILE("vof stencil definition");
      const     Box& gridCoFi = m_eblgCoFi.getDBL()  [mfi];
      const EBISBox& ebisCoFi = m_eblgCoFi.getEBISL()[mfi];
 
      IntVectSet ivs = ebisCoFi.getIrregIVS(gridCoFi);
      VoFIterator vofit(ivs, ebisCoFi.getEBGraph());
      const std::vector<VolIndex>& vofvec = vofit.getVector();
      // cast from VolIndex to BaseIndex
      std::vector<std::shared_ptr<BaseIndex> >    dstVoF(vofvec.size());
      std::vector<std::shared_ptr<BaseStencil> > stencil(vofvec.size());
      std::vector<VoFStencil>  allsten(vofvec.size());
      // fill stencils for the vofs
      for(int ivec = 0; ivec < vofvec.size(); ivec++)
      {
        definePointStencil(allsten[ivec], vofvec[ivec], mfi);

        // another cast from VolIndex to BaseIndex
        dstVoF[ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&vofvec[ivec]), &null_deleter_ebc_ind);
        stencil[ivec] = std::shared_ptr<BaseStencil>(            &allsten[ivec] , &null_deleter_ebc_sten);
      }
      //the fine  is the source, the coarsened fine is the destination
      m_stencil[mfi] = std::shared_ptr<AggStencil<EBCellFAB, EBCellFAB > >
        (new AggStencil<EBCellFAB, EBCellFAB >(dstVoF, stencil, dummyFine[mfi], dummyCoFi[mfi]));

    }
    
  }
  /************************************/
  void
  EBCoarseAverage::
  definePointStencil(VoFStencil& a_sten, const VolIndex& a_vofCoFi, const MFIter& a_mfi)
  {
    std::vector<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(a_vofCoFi, m_refRat, a_mfi);
    a_sten.clear();
    Real sumKappa = 0;
    for(int ivof = 0; ivof < fineVoFs.size(); ivof++)
    {
      Real weight = 1; //correct answer unless we are using kappa weighting
      const VolIndex & fineVoF = fineVoFs[ivof];
      if(m_useKappaWeightingInStencil)
      {
        Real volFrac = m_eblgFine.getEBISL()[a_mfi].volFrac(fineVoF);
        weight = volFrac;
        sumKappa += volFrac;
      }
      a_sten.add(fineVoF, weight);
    }
    if(m_useKappaWeightingInStencil)
    {
      if(sumKappa > 1.0e-10)
      {
        a_sten *= (1.0/sumKappa);
      }
      else
      {
        int numFineCellsPerCoarse = D_TERM(m_refRat, *m_refRat, *m_refRat);
        a_sten *= (1.0/Real(numFineCellsPerCoarse));
      }
    }
  }
  /************************************/
  void
  EBCoarseAverage::
  average(FabArray<EBCellFAB>       & a_dataCoar,
          const FabArray<EBCellFAB> & a_dataFine,
          int isrc, int idst, int inco)
  {
    BL_PROFILE("EBCoarseAverage::average(LD<EBCellFAB>)");
    shared_ptr<MultiFab> regCoarData, regFineData;

    EBLevelDataOps::aliasIntoMF(regCoarData, a_dataCoar, m_eblgCoar);
    EBLevelDataOps::aliasIntoMF(regFineData, a_dataFine, m_eblgFine);
    
    //the dx here better not matter so make a fake realbox
    RealVect rvlo = RealVect::Zero;
    RealVect rvhi = RealVect::Unit;
    RealBox rb(rvlo.dataPtr(), rvhi.dataPtr());

    Geometry geomCoar(m_eblgCoar.getDomain(), &rb);
    Geometry geomFine(m_eblgFine.getDomain(), &rb);
    
    //from multifab_util--this averages over all cells as if EB were not there.
    average_down(*regFineData, *regCoarData, isrc, inco, m_refRat);

    //now replace that answer at the irregular cells.
    EBCellFactory fact(m_eblgCoFi.getEBISL());
    int nvar = inco+idst;
    FabArray<EBCellFAB> dataCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_ghost, MFInfo(), fact);
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      //the false is to turn off incrementOnly
      m_stencil[mfi]->apply(dataCoFi[mfi], a_dataFine[mfi], isrc, idst, inco, false);
    }

    //now copy into the data holder on the designated layout
    //do not need to copy ghost data
    a_dataCoar.copy(a_dataCoar, idst, idst, inco, 0, 0);

  }

  /************************************/
}
