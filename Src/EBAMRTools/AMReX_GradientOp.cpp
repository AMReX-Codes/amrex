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


#include "AMReX_GradientOp.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBLoHiCenter.H"
#include "AMReX_EBFluxFactory.H"
#include "AMReX_IrregFABFactory.H"
#include "AMReX_EBFortND_F.H"


namespace amrex
{
  void null_deleter_grads_sten(BaseStencil * a_sten)
  {}
  void null_deleter_grads_ind(BaseIndex* a_sten)
  {}
  /************************************/
  void
  GradientOp::
  define(const EBLevelGrid   & a_eblg,
         const Real          & a_dx,
         const int           & a_nComp,
         const int           & a_ghostCellsInData,
         bool a_useLimiting,
         bool a_slowMode)
  {
    m_isDefined     = true;
    m_slowMode      = a_slowMode;
    m_eblg          = a_eblg;
    m_dx            = a_dx;
    m_nComp         = a_nComp;
    m_dataGhost     = a_ghostCellsInData;
    m_useLimiting   = a_useLimiting;
    defineInternals();
  }
  /************************************/
  void
  GradientOp::
  defineInternals()
  {
    BL_PROFILE("NWOEBCFI::defineInternals");

    //variable number does not matter here.
    EBCellFactory  ebcellfact(m_eblg.getEBISL());
    IrregFABFactory irregfact(m_eblg.getEBISL());
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_loSlope[idir].define(m_eblg.getDBL(), m_eblg.getDM(), m_nComp, m_dataGhost, MFInfo(), irregfact);
      m_hiSlope[idir].define(m_eblg.getDBL(), m_eblg.getDM(), m_nComp, m_dataGhost, MFInfo(), irregfact);
    }
    FabArray<EBCellFAB> cellProxy(m_eblg.getDBL(), m_eblg.getDM(), m_nComp, m_dataGhost, MFInfo(), ebcellfact);

    Box domain = m_eblg.getDomain();

    if(m_slowMode)
    {
      m_irregVoFsSlow.define(m_eblg.getDBL(), m_eblg.getDM());
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        m_loStencilSlow[idir].define(m_eblg.getDBL(), m_eblg.getDM());
        m_hiStencilSlow[idir].define(m_eblg.getDBL(), m_eblg.getDM());
      }
    }
    else
    {
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        m_loStencil[idir].define(m_eblg.getDBL(), m_eblg.getDM());
        m_hiStencil[idir].define(m_eblg.getDBL(), m_eblg.getDM());
      }
    }

    m_vofit.define(m_eblg.getDBL(), m_eblg.getDM());
    
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const EBISBox  & ebis =    m_eblg.getEBISL()[mfi];
      Box              grid =    m_eblg.getDBL()  [mfi];
      IntVectSet ivsIrreg = ebis.getIrregIVS(grid);
      VoFIterator & vofit = m_vofit[mfi];
      vofit.define(ivsIrreg, ebis.getEBGraph());
      const Vector<VolIndex>& volvec = vofit.getVector();

      //destination vofs 
      Vector< std::shared_ptr<BaseIndex  > > baseDstVoFs(volvec.size());
      for(int ivec = 0; ivec < volvec.size(); ivec++)
      {
        baseDstVoFs [ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&volvec[ivec]), &null_deleter_grads_ind);
      }

      if(m_slowMode)
      {
        m_irregVoFsSlow[mfi] = volvec;
      }

      // stencils for slopes in irregular cells
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        Vector< std::shared_ptr<BaseStencil> > baseStenLo(volvec.size());
        Vector< std::shared_ptr<BaseStencil> > baseStenHi(volvec.size());
        Vector<VoFStencil> allVoFStenLo(volvec.size());
        Vector<VoFStencil> allVoFStenHi(volvec.size());

        for(int ivec = 0; ivec < volvec.size(); ivec++)
        {
          
          const VolIndex& vof = volvec[ivec];
          getStencils(allVoFStenLo[ivec], allVoFStenHi[ivec], vof, idir, ebis);

          baseStenLo[ivec] = std::shared_ptr<BaseStencil>(&allVoFStenLo[ivec] , &null_deleter_grads_sten);
          baseStenHi[ivec] = std::shared_ptr<BaseStencil>(&allVoFStenHi[ivec] , &null_deleter_grads_sten);
        }
        if(m_slowMode)
        {
          m_loStencilSlow[idir][mfi] = allVoFStenLo;
          m_hiStencilSlow[idir][mfi] = allVoFStenHi;
        }
        else
        {
          m_loStencil[idir][mfi] = std::shared_ptr<AggStencil <EBCellFAB, IrregFAB>  >
            (new AggStencil<EBCellFAB, IrregFAB >(baseDstVoFs, baseStenLo, cellProxy[mfi], m_loSlope[idir][mfi]));
          m_hiStencil[idir][mfi] = std::shared_ptr<AggStencil <EBCellFAB, IrregFAB>  >
            (new AggStencil<EBCellFAB, IrregFAB >(baseDstVoFs, baseStenHi, cellProxy[mfi], m_hiSlope[idir][mfi]));
        }
      }
    }
  }
  
  /*******/
  void
  GradientOp::
  getStencils(VoFStencil         &   a_loSten,
              VoFStencil         &   a_hiSten,
              const VolIndex     &   a_vof,
              const int          &   a_idir,
              const EBISBox      &   a_ebis)
  {
              
    Vector<FaceIndex> facesLo = a_ebis.getFaces(a_vof, a_idir, Side::Lo);
    Vector<FaceIndex> facesHi = a_ebis.getFaces(a_vof, a_idir, Side::Hi);

    bool hasFacesLo = ((facesLo.size() == 1) && (!facesLo[0].isBoundary()));
    bool hasFacesHi = ((facesHi.size() == 1) && (!facesHi[0].isBoundary()));
    if(hasFacesLo && hasFacesHi)
    {
      VolIndex loVoF = facesLo[0].getVoF(Side::Lo);
      VolIndex hiVoF = facesHi[0].getVoF(Side::Hi);
      a_loSten.add(a_vof,  1.0/m_dx);
      a_loSten.add(loVoF, -1.0/m_dx);
      a_hiSten.add(hiVoF,  1.0/m_dx);
      a_hiSten.add(a_vof, -1.0/m_dx);
    }
    else if(hasFacesLo)
    {
      const VolIndex& loVoF = facesLo[0].getVoF(Side::Lo);
      Vector<FaceIndex> facesLower = a_ebis.getFaces(loVoF, a_idir, Side::Lo);
      bool hasLower = ((facesLower.size() == 1) && (!facesLower[0].isBoundary()));
      if(hasLower)
      {
        VolIndex lowerVoF = facesLower[0].getVoF(Side::Lo);
        a_hiSten.add(a_vof   ,  1.0/m_dx);
        a_hiSten.add(loVoF   , -1.0/m_dx);
        a_loSten.add(a_vof   ,  1.0/(2.*m_dx));
        a_loSten.add(lowerVoF, -1.0/(2.*m_dx));
      }
      else
      {
        a_hiSten.add(a_vof   ,  1.0/m_dx);
        a_hiSten.add(loVoF   , -1.0/m_dx);
        a_loSten.add(a_vof   ,  1.0/m_dx);
        a_loSten.add(loVoF   , -1.0/m_dx);
      }
    }
    else if(hasFacesHi)
    {
      const VolIndex& hiVoF = facesHi[0].getVoF(Side::Hi);
      Vector<FaceIndex> facesHigher = a_ebis.getFaces(hiVoF, a_idir, Side::Hi);
      bool hasHigher = ((facesHigher.size() == 1) && (!facesHigher[0].isBoundary()));
      if(hasHigher)
      {
        VolIndex higherVoF = facesHigher[0].getVoF(Side::Hi);
        a_hiSten.add(higherVoF,  1.0/(2.*m_dx));
        a_hiSten.add(    a_vof, -1.0/(2.*m_dx));
        a_loSten.add(    hiVoF,  1.0/m_dx);
        a_loSten.add(    a_vof, -1.0/m_dx);
      }
      else
      {
        a_hiSten.add(    hiVoF,  1.0/m_dx);
        a_hiSten.add(    a_vof, -1.0/m_dx);
        a_loSten.add(    hiVoF,  1.0/m_dx);
        a_loSten.add(    a_vof, -1.0/m_dx);
      }
    }
    else
    {
      //cannot reliably reach in either direction.  just set slopes to zero
      a_loSten.clear();
      a_hiSten.clear();
    }
  }
  /************************************/
  void
  GradientOp::
  gradient(FabArray<EBCellFAB>      & a_gph,
           const FabArray<EBCellFAB>& a_phi)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_phi.nGrow() == m_dataGhost);
    BL_ASSERT(a_gph.nComp() == m_nComp*SpaceDim);
    BL_ASSERT(a_phi.nComp() == m_nComp);

    FabArray<EBCellFAB>& castphi = const_cast<FabArray<EBCellFAB> &>(a_phi);
    castphi.FillBoundary();
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      EBCellFAB       & gph = a_gph[mfi];
      const EBCellFAB & phi = a_phi[mfi];

      BaseFab<Real>       &  regGph = gph.getSingleValuedFAB();
      const BaseFab<Real> &  regPhi = phi.getSingleValuedFAB();;

      const Box& grid = m_eblg.getDBL()[mfi];
      
      int doLimiting = 0;
      if(m_useLimiting) doLimiting = 1;
      //first do everything as if it has no eb.
      ebfnd_gradlim(BL_TO_FORTRAN_FAB(regGph),
                    BL_TO_FORTRAN_FAB(regPhi),
                    BL_TO_FORTRAN_BOX(grid),
                    &m_nComp, &m_dx, &doLimiting);


      //fill our slope holders with low and high slopes
      if(m_slowMode)
      {
        for(int ivar = 0; ivar < m_nComp; ivar++)
        {
          for(int vecdir = 0; vecdir < SpaceDim; vecdir++)
          {
            const Vector<VolIndex>& vofvec = m_irregVoFsSlow[mfi];
            for(int ivof = 0; ivof < vofvec.size(); ivof++)
            {
              const VoFStencil& vofstenLo = m_loStencilSlow[vecdir][mfi][ivof];
              const VoFStencil& vofstenHi = m_hiStencilSlow[vecdir][mfi][ivof];
              Real loSlope = applyVoFStencil(vofstenLo, phi);
              Real hiSlope = applyVoFStencil(vofstenHi, phi);
              m_loSlope[vecdir][mfi](vofvec[ivof], ivar) = loSlope;
              m_hiSlope[vecdir][mfi](vofvec[ivof], ivar) = hiSlope;
            }
          }
        } 
      }
      else
      {
        for(int ivar = 0; ivar < m_nComp; ivar++)
        {
          for(int vecdir = 0; vecdir < SpaceDim; vecdir++)
          {
            int isrc = ivar;
            int idst = ivar;
            int inco = 1;
            bool incrementOnly = false;
            m_loStencil[vecdir][mfi]->apply(m_loSlope[vecdir][mfi], phi, isrc, idst, inco, incrementOnly);
            m_hiStencil[vecdir][mfi]->apply(m_hiSlope[vecdir][mfi], phi, isrc, idst, inco, incrementOnly);
          }
        } 
      }


      VoFIterator vofit = m_vofit[mfi];
      for(vofit.reset(); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit(); 
        for(int ivar = 0; ivar < m_nComp; ivar++)
        {
          for(int vecdir = 0; vecdir < SpaceDim; vecdir++)
          {
            int idst = gradIndex(ivar, vecdir);
            Real deltahi = m_hiSlope[vecdir][mfi](vof, ivar);
            Real deltalo = m_loSlope[vecdir][mfi](vof, ivar);

            Real mono = deltahi*deltalo;
            Real deltaminmod;
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
            a_gph[mfi](vof, idst) = deltaminmod;
          }
        }
      }
    }
  }
  /************************************/
}

