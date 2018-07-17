#include "AMReX_DivergenceOp.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBCellFactory.H"
#include "AMReX_EBLoHiCenter.H"
#include "AMReX_EBFluxFactory.H"
#include "AMReX_IrregFABFactory.H"
#include "AMReX_EBFortND_F.H"


namespace amrex
{
  void null_deleter_divs_sten(BaseStencil * a_sten)
  {}
  void null_deleter_divs_ind(BaseIndex* a_sten)
  {}
  /************************************/
  void
  DivergenceOp::
  define(const EBLevelGrid   & a_eblg,
         const Real          & a_dx,
         const int           & a_nComp,
         const int           & a_ghostCellsInData,
         bool a_multiplyFluxByArea,
         int a_redistRad)
  {
    m_isDefined = true;
    m_eblg          = a_eblg;
    m_multiplyFluxByArea = a_multiplyFluxByArea;
    m_dx            = a_dx;
    m_nComp         = a_nComp;
    m_dataGhost     = a_ghostCellsInData;
    m_redistRad     = a_redistRad;
    defineInternals();
  }
  /************************************/
  void
  DivergenceOp::
  defineInternals()
  {
    BL_PROFILE("NWOEBCFI::defineInternals");


    m_eblevelRedist.define(m_eblg, m_nComp, m_redistRad);
    m_normalizor.define(m_eblg, m_dataGhost);

    //variable number does not matter here.
    EBCellFactory  ebcellfact(m_eblg.getEBISL());
    EBFluxFactory  ebfluxfact(m_eblg.getEBISL());
    IrregFABFactory irregfact(m_eblg.getEBISL());
    m_kappaDivergence.define(     m_eblg.getDBL(), m_eblg.getDM(), m_nComp, m_dataGhost, MFInfo(), ebcellfact);
    m_massDiff.define(            m_eblg.getDBL(), m_eblg.getDM(), m_nComp, m_dataGhost, MFInfo(),  irregfact);
    FabArray<EBFluxFAB> fluxProxy(m_eblg.getDBL(), m_eblg.getDM(), m_nComp, m_dataGhost, MFInfo(), ebfluxfact);
    //have to define kappaDivergence anyway so use it as a proxy
    FabArray<EBCellFAB>& cellProxy = m_kappaDivergence;

    Box domain = m_eblg.getDomain();

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_openStencil[idir].define(m_eblg.getDBL(), m_eblg.getDM());
    }
    m_bdryStencil.define(m_eblg.getDBL(), m_eblg.getDM());
    m_vofit.define(m_eblg.getDBL(), m_eblg.getDM());
    
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      const EBISBox  & ebis =    m_eblg.getEBISL()[mfi];
      Box              grid =    m_eblg.getDBL()  [mfi];
      IntVectSet ivsIrreg = ebis.getIrregIVS(grid);
      VoFIterator & vofit = m_vofit[mfi];
      vofit.define(ivsIrreg, ebis.getEBGraph());
      const Vector<VolIndex>& volvec = vofit.getVector();

      //destination vofs are the same for both open and boundary faces
      Vector< std::shared_ptr<BaseIndex  > > baseDstVoFs(volvec.size());
      for(int ivec = 0; ivec < volvec.size(); ivec++)
      {
        baseDstVoFs [ivec]  = std::shared_ptr<BaseIndex  >((BaseIndex*)(&volvec[ivec]), &null_deleter_divs_ind);
      }

      //THIS IS NOT MEANT TO WORK IN CASES WHERE EMBEDDED BOUNDARIES CROSS COARSE FINE BOUNDARIES
      //first let us dal with the boundary flux stencil
      {
        Vector< std::shared_ptr<BaseStencil> > baseSten(volvec.size());
        Vector<VoFStencil> allvofsten(volvec.size());
        for(int ivec = 0; ivec < volvec.size(); ivec++)
        {
          Real bndryArea = ebis.bndryArea(volvec[ivec]);
          VoFStencil bndryStencil;
          const VolIndex& vof = volvec[ivec];
          Real weight = bndryArea/m_dx;
          allvofsten[ivec].clear();
          allvofsten[ivec].add(vof, weight);
          baseSten[ivec]  = std::shared_ptr<BaseStencil>(&allvofsten[ivec] , &null_deleter_divs_sten);
        }
        m_bdryStencil[mfi] = std::shared_ptr<AggStencil <IrregFAB, EBCellFAB>  >
          (new AggStencil<IrregFAB, EBCellFAB >(baseDstVoFs, baseSten, fluxProxy[mfi].getEBFlux(), cellProxy[mfi]));
      }

      Real full_vol = 1;
      for (int idir = 0; idir < SpaceDim; ++idir)
      {
        full_vol *= m_dx;
      }
      Real inv_vol = 1/full_vol;
      // now do open flux stencils for each face direction
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        Vector< std::shared_ptr<BaseStencil> > baseSten(volvec.size());
        Vector<FaceStencil> allFaceSten(volvec.size());
        for(int ivec = 0; ivec < volvec.size(); ivec++)
        {
          
          const VolIndex& vof = volvec[ivec];
          FaceStencil& dirStencil = allFaceSten[ivec];
          dirStencil.clear();
          for(SideIterator sit; sit.ok(); ++sit)
          {
            Vector<FaceIndex> faces = ebis.getFaces(vof, idir, sit());
            int isign = sign(sit());
            for(int iface = 0; iface < faces.size(); iface++)
            {
              IntVectSet cfivs;//empty--see comment above
              FaceStencil interpSten = EBArith::getInterpStencil(faces[iface], cfivs, ebis, domain);
              Real areaFrac = ebis.areaFrac(faces[iface]);
              if(m_multiplyFluxByArea)
              {
                interpSten *= (isign*areaFrac/m_dx);
              }
              else
              {
                interpSten *= isign * areaFrac * inv_vol;
              }
              dirStencil += interpSten;
            }
          }
          baseSten[ivec]  = std::shared_ptr<BaseStencil>(&allFaceSten[ivec] , &null_deleter_divs_sten);
        }
        m_openStencil[idir][mfi] = std::shared_ptr<AggStencil <EBFaceFAB, EBCellFAB>  >
          (new AggStencil<EBFaceFAB, EBCellFAB >(baseDstVoFs, baseSten, fluxProxy[mfi][idir], cellProxy[mfi]));
      }
    }

  }
  /************************************/
  void
  DivergenceOp::
  hybridDivergence(FabArray<EBCellFAB>      & a_divF,
                   const FabArray<EBFluxFAB>& a_flux,
                   int isrc, int idst, int inco,
                   bool a_trustRegDivF)
  {
/**/
    BL_ASSERT(isDefined());
    BL_ASSERT(a_flux.nGrow() == m_dataGhost);
    BL_ASSERT(a_divF.nGrow() == m_dataGhost);

    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      EBCellFAB       & divF = m_kappaDivergence[mfi];
      const EBFluxFAB & flux = a_flux[mfi];

      if(a_trustRegDivF)
      {
        divF.copy(a_divF[mfi]);
      }
      else
      {
        BaseFab<Real>       &  regDivF = divF.getSingleValuedFAB();
        Vector<const BaseFab<Real>*> regFlux(3, &(flux[0].getSingleValuedFAB()));
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          regFlux[idir] = &(flux[idir].getSingleValuedFAB());
        }
        const Box& grid = m_eblg.getDBL()[mfi];
      
        int multiplyByArea = 0;
        if(m_multiplyFluxByArea) multiplyByArea = 1;

        //first do everything as if it has no eb.
        ebfnd_divflux(BL_TO_FORTRAN_FAB(regDivF),
                      BL_TO_FORTRAN_FAB((*regFlux[0])),
                      BL_TO_FORTRAN_FAB((*regFlux[1])),
                      BL_TO_FORTRAN_FAB((*regFlux[2])),
                      BL_TO_FORTRAN_BOX(grid),
                      &multiplyByArea,
                      &m_dx, &isrc, &idst, &inco);
      }
      //turn off increment only for bndry flux.  This sets the initial divergence at 
      //cut cells to zero
      bool incrementOnly = false;
      m_bdryStencil[mfi]->apply(divF, flux.getEBFlux(), isrc, idst, inco, incrementOnly);
      //turn on incr only so that now the other terms will just get added in
      incrementOnly = true;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        m_openStencil[idir][mfi]->apply(divF, flux[idir], isrc, idst, inco, incrementOnly);
      } 
    }

    //at this point kappaDivergence holds, well, kappa*divF.   Now normalize it 
    //so it will now hold DivFNC for now.
    m_normalizor.normalize(a_divF, m_kappaDivergence, idst, inco);


    //now behold the joys of redistribution.
    m_eblevelRedist.setToZero();
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      VoFIterator vofit = m_vofit[mfi];
      for(vofit.reset(); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        int varmax = idst + inco;
        for(int ivar = idst; ivar < varmax; ivar++)
        {
          Real kappa = m_eblg.getEBISL()[mfi].volFrac(vof);
          const  Real& kappaDivC = m_kappaDivergence[mfi](vof, ivar);
          const  Real& divFNC    =            a_divF[mfi](vof, ivar);
          m_massDiff[mfi](vof, ivar) = (1. - kappa)*(kappaDivC-kappa*divFNC);
          a_divF[mfi](vof, ivar) = kappaDivC + (1.-kappa)*divFNC;
        }
      }
      m_eblevelRedist.increment(m_massDiff[mfi], mfi, idst, inco);
    }
    m_eblevelRedist.redistribute(a_divF, idst, inco);
/**/
  }
  /************************************/
}

