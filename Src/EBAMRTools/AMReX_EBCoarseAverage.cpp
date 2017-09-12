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
#include "AMReX_EBFortND_F.H"


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
                  const bool        & a_useKappaWeightingInStencil,
                  bool a_enableFaceAveraging)
  {
    define(a_eblgFine, a_eblgCoar, a_nref, a_ghost, a_useKappaWeightingInStencil, a_enableFaceAveraging);
  }
  /************************************/
  void
  EBCoarseAverage::
  define(const EBLevelGrid & a_eblgFine,
         const EBLevelGrid & a_eblgCoar,
         const int         & a_nref,
         const int         & a_ghost,
         const bool        & a_useKappaWeightingInStencil,
         bool a_enableFaceAveraging)
  {
    BL_PROFILE("EBCoarseAverage::define(EBLG)");
    BL_ASSERT(a_nref > 0);
    BL_ASSERT(a_eblgFine.getEBISL().getGhost() >= 2);
    BL_ASSERT(a_eblgCoar.getEBISL().getGhost() >= 2);
    BL_ASSERT(a_eblgFine.coarsenable(a_nref));
    
    m_isDefined = true;
    m_enableFaceAveraging = a_enableFaceAveraging;
    m_refRat    = a_nref;
    m_eblgCoar  = a_eblgCoar;
    m_eblgFine  = a_eblgFine;
    m_ghost    = a_ghost;
    m_useKappaWeightingInStencil = a_useKappaWeightingInStencil;
    coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
    m_eblgCoFi.setMaxRefinementRatio(a_nref);
    defineStencils();
    if(m_enableFaceAveraging)
    {
      defineFaceStencils();
    }
  }
  /************************************/
  void
  EBCoarseAverage::
  defineStencils()
  {
    BL_PROFILE("ebcoarseaverage::defineStencils");
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
      const Array<VolIndex>& vofvec = vofit.getVector();
      // cast from VolIndex to BaseIndex
      Array<std::shared_ptr<BaseIndex> >    dstVoF(vofvec.size());
      Array<std::shared_ptr<BaseStencil> > stencil(vofvec.size());
      Array<VoFStencil>  allsten(vofvec.size());
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
  defineFaceStencils()
  {
    BL_PROFILE("ebcoarseaverage::defineFaceStencils");
    BL_ASSERT(m_enableFaceAveraging);

    EBFluxFactory factFine(m_eblgFine.getEBISL());
    EBFluxFactory factCoFi(m_eblgCoFi.getEBISL());
    int nvar = 1; //not important
    FabArray<EBFluxFAB> dummyFine(m_eblgFine.getDBL(),m_eblgFine.getDM(), nvar, m_ghost, MFInfo(), factFine);
    FabArray<EBFluxFAB> dummyCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_ghost, MFInfo(), factCoFi);

    //first the stencils for the coordinate faces
    for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      m_faceStencil[faceDir].define(m_eblgFine.getDBL(), m_eblgFine.getDM());
      for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
      {
        BL_PROFILE("face stencil definition");
        const     Box& gridCoFi = m_eblgCoFi.getDBL()  [mfi];
        const EBISBox& ebisCoFi = m_eblgCoFi.getEBISL()[mfi];
 
        IntVectSet ivs = ebisCoFi.getIrregIVS(gridCoFi);
        FaceIterator faceit(ivs, ebisCoFi.getEBGraph(), faceDir, FaceStop::SurroundingWithBoundary);
        const Array<FaceIndex>& facevec = faceit.getVector();
        // cast from FaceIndex to BaseIndex
        Array<std::shared_ptr<BaseIndex> >   dstFace(facevec.size());
        Array<std::shared_ptr<BaseStencil> > stencil(facevec.size());
        Array<FaceStencil>  allsten(facevec.size());
        // fill stencils for the vofs
        for(int ivec = 0; ivec < facevec.size(); ivec++)
        {
          definePointStencilFace(allsten[ivec], facevec[ivec], mfi);

          dstFace[ivec] = std::shared_ptr<BaseIndex  >((BaseIndex*)(&facevec[ivec]), &null_deleter_ebc_ind);
          stencil[ivec] = std::shared_ptr<BaseStencil>(             &allsten[ivec] , &null_deleter_ebc_sten);
        }
        //the fine  is the source, the coarsened fine is the destination
        m_faceStencil[faceDir][mfi] = std::shared_ptr<AggStencil<EBFaceFAB, EBFaceFAB > >
          (new AggStencil<EBFaceFAB, EBFaceFAB >(dstFace, stencil, dummyFine[mfi][faceDir], dummyCoFi[mfi][faceDir]));

      }
    
    }
    //now for the cut face stencils (irregular faces)
    m_irregStencil.define(m_eblgFine.getDBL(), m_eblgFine.getDM());
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      BL_PROFILE("irreg stencil definition");
      const     Box& gridCoFi = m_eblgCoFi.getDBL()  [mfi];
      const     Box& gridFine = m_eblgFine.getDBL()  [mfi];
      const EBISBox& ebisCoFi = m_eblgCoFi.getEBISL()[mfi];
      const EBISBox& ebisFine = m_eblgFine.getEBISL()[mfi];

      IntVectSet ivsFine = ebisFine.getIrregIVS(gridFine);
      IntVectSet ivs = ebisCoFi.getIrregIVS(gridCoFi);
      VoFIterator vofit(ivs, ebisCoFi.getEBGraph());
      const Array<VolIndex>& vofvec = vofit.getVector();

      // cast from VolIndex to BaseIndex
      Array<std::shared_ptr<BaseIndex> >    dstVoF(vofvec.size());
      Array<std::shared_ptr<BaseStencil> > stencil(vofvec.size());
      Array<VoFStencil>  allsten(vofvec.size());
      // fill stencils for the vofs
      for(int ivec = 0; ivec < vofvec.size(); ivec++)
      {
          definePointStencilIrreg(allsten[ivec], vofvec[ivec], mfi, ivsFine);

        // another cast from VolIndex to BaseIndex
        dstVoF[ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&vofvec[ivec]), &null_deleter_ebc_ind);
        stencil[ivec] = std::shared_ptr<BaseStencil>(            &allsten[ivec] , &null_deleter_ebc_sten);
      }
      //the fine  is the source, the coarsened fine is the destination
      m_irregStencil[mfi] = std::shared_ptr<AggStencil<IrregFAB, IrregFAB > >
        (new AggStencil<IrregFAB, IrregFAB>(dstVoF, stencil, dummyFine[mfi].getEBFlux(), dummyCoFi[mfi].getEBFlux()));

    }
  }
  /************************************/
  void
  EBCoarseAverage::
  definePointStencil(VoFStencil& a_sten, const VolIndex& a_vofCoFi, const MFIter& a_mfi)
  {
    Array<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(a_vofCoFi, m_refRat, a_mfi);
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
        a_sten *= 0.0;
      }
    }
    else
    {
      int numFineCellsPerCoarse = AMREX_D_TERM(m_refRat, *m_refRat, *m_refRat);
      a_sten *= (1.0/Real(numFineCellsPerCoarse));
    }
  }
  /************************************/
  void
  EBCoarseAverage::
  definePointStencilIrreg(VoFStencil& a_sten, 
                          const VolIndex& a_vofCoFi, 
                          const MFIter& a_mfi,
                          const IntVectSet& a_ivsIrregFine)
  {
    BL_ASSERT(m_enableFaceAveraging);
    Array<VolIndex> fineVoFs = m_eblgCoFi.getEBISL().refine(a_vofCoFi, m_refRat, a_mfi);
    a_sten.clear();
    Real sumBndryArea = 0;
    for(int ivof = 0; ivof < fineVoFs.size(); ivof++)
    {
      const VolIndex & fineVoF = fineVoFs[ivof];
      if(a_ivsIrregFine.contains(fineVoF.gridIndex())) //only want cells that are actually cut
      {
        Real area = m_eblgFine.getEBISL()[a_mfi].bndryArea(fineVoF);
        Real weight = area;
        sumBndryArea += area;

        a_sten.add(fineVoF, weight);
      }
    }

    if(sumBndryArea > 1.0e-10)
    {
      a_sten *= (1.0/sumBndryArea);
    }
    else
    {
      a_sten *= 0.0;
    }
  }
  /************************************/
  void
  EBCoarseAverage::
  definePointStencilFace(FaceStencil& a_sten, const FaceIndex& a_faceCoFi, const MFIter& a_mfi)
  {
    BL_ASSERT(m_enableFaceAveraging);
    Array<FaceIndex> fineFaces = m_eblgCoFi.getEBISL().refine(a_faceCoFi, m_refRat, a_mfi);
    a_sten.clear();
    Real sumArea = 0;
    for(int iface = 0; iface < fineFaces.size(); iface++)
    {
      const FaceIndex & fineFace = fineFaces[iface];
      Real areaFrac = m_eblgFine.getEBISL()[a_mfi].areaFrac(fineFace);
      Real weight = areaFrac;
      sumArea += areaFrac;
      a_sten.add(fineFace, weight);
    }

    if(sumArea > 1.0e-10)
    {
      a_sten *= (1.0/sumArea);
    }
    else
    {
      a_sten *= 0.0;
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

    
    EBCellFactory fact(m_eblgCoFi.getEBISL());
    int nvar = inco+idst;
    FabArray<EBCellFAB> dataCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_ghost, MFInfo(), fact);
    Box refbox(IntVect::Zero, IntVect::Zero);
    refbox.refine(m_refRat);
    int ibox = 0;
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      BaseFab<Real>&       regCoar =   dataCoFi[mfi].getSingleValuedFAB();
      const BaseFab<Real>& regFine = a_dataFine[mfi].getSingleValuedFAB();
      Box coarbx = m_eblgCoFi.getDBL()[mfi];

      ebfnd_average(BL_TO_FORTRAN_FAB(regCoar),
                    BL_TO_FORTRAN_FAB(regFine),
                    BL_TO_FORTRAN_BOX(coarbx),
                    BL_TO_FORTRAN_BOX(refbox),
                    &m_refRat, &isrc, &idst, &inco);

      //now do the irregular bit
      //(replace that answer at the irregular cells)
      //the false is to turn off incrementOnly
      m_stencil[mfi]->apply(dataCoFi[mfi], a_dataFine[mfi], isrc, idst, inco, false);
      ibox++;
    }

    //now copy into the data holder on the designated layout
    //do not need to copy ghost data
    a_dataCoar.copy(dataCoFi, idst, idst, inco, 0, 0);

  }
  /************************************/
  void
  EBCoarseAverage::
  average(FabArray<EBFluxFAB>       & a_dataCoar,
          const FabArray<EBFluxFAB> & a_dataFine,
          int isrc, int idst, int inco)
  {
    BL_PROFILE("EBCoarseAverage::average(LD<EBFluxFAB>)");

    BL_ASSERT(m_enableFaceAveraging);
    
    EBFluxFactory fact(m_eblgCoFi.getEBISL());
    int nvar = inco+idst;
    FabArray<EBFluxFAB> dataCoFi(m_eblgCoFi.getDBL(),m_eblgCoFi.getDM(), nvar, m_ghost, MFInfo(), fact);
    Box refbox(IntVect::Zero, IntVect::Zero);
    refbox.refine(m_refRat);
    int ibox = 0;
    for(MFIter mfi(m_eblgFine.getDBL(), m_eblgFine.getDM()); mfi.isValid(); ++mfi)
    {
      //first the coordinate faces
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        EBFaceFAB&       coarface =   dataCoFi[mfi][faceDir];
        const EBFaceFAB& fineface = a_dataFine[mfi][faceDir];
        BaseFab<Real>&       regCoar =   coarface.getSingleValuedFAB();
        const BaseFab<Real>& regFine =   fineface.getSingleValuedFAB();
        Box gridbx = m_eblgCoFi.getDBL()[mfi];
        
        Box coarbx =  surroundingNodes(gridbx, faceDir);
        Box refFace = refbox;
        refFace.growHi(faceDir, -(m_refRat-1));
        
        ebfnd_average_face(BL_TO_FORTRAN_FAB(regCoar),
                           BL_TO_FORTRAN_FAB(regFine),
                           BL_TO_FORTRAN_BOX(coarbx),
                           BL_TO_FORTRAN_BOX(refFace),
                           &faceDir, &m_refRat, &isrc, &idst, &inco);

        //now do the irregular bit
        //(replace that answer at the irregular cells)
        //the false is to turn off incrementOnly
        m_faceStencil[faceDir][mfi]->apply(dataCoFi[mfi][faceDir], a_dataFine[mfi][faceDir], isrc, idst, inco, false);
        ibox++;
      }

      //now do cut face data
      m_irregStencil[mfi]->apply(dataCoFi[mfi].getEBFlux(), a_dataFine[mfi].getEBFlux(), isrc, idst, inco, false);
    }
    //now copy into the data holder on the designated layout
    //do not need to copy ghost data
    a_dataCoar.copy(dataCoFi, idst, idst, inco, 0, 0);

  }

  /************************************/
}
