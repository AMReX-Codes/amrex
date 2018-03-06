#include "AMReX_EBNormalizeByVolumeFraction.H"
#include "AMReX_EBArith.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFactory.H"

namespace amrex
{
  void null_deleter_vsten(BaseStencil * a_sten)
 {}
  void null_deleter_vof(BaseIndex* a_sten)
 {}
//----------------------------------------------------------------------------
  void 
  EBNormalizeByVolumeFraction::
  getLocalStencil(VoFStencil      & a_stencil, 
                  const VolIndex  & a_vof, 
                  const MFIter    & a_dit)
  {
    Vector<VolIndex> neighbors;
    EBArith::getAllVoFsWithinRadius(neighbors, a_vof, m_eblg.getEBISL()[a_dit], m_radius);
    Real sumkap = 0;
    a_stencil.clear();
    for(int ivof = 0; ivof < neighbors.size(); ivof++)
    {
      Real kappa = m_eblg.getEBISL()[a_dit].volFrac(neighbors[ivof]);
      //data is assumed to already be kappa weighted
      a_stencil.add(neighbors[ivof], 1.0);
      sumkap += kappa;
    }
    if(sumkap > 1.0e-9)
    {
      a_stencil *= (1.0/sumkap);
    }
  }
//----------------------------------------------------------------------------
  void 
  EBNormalizeByVolumeFraction::
  define(const EBLevelGrid          & a_eblg,
         const int                  & a_ghost,
         const int                    a_radius)
  {
    BL_PROFILE("EBNormalizeByVolFrac::define");
    BL_ASSERT(a_radius > 0);
    m_eblg      = a_eblg;
    m_radius    = a_radius;
    m_ghost     = a_ghost;
    m_isDefined = true;

    int nvar = 1; //this object can be used for any number of variables.  This is just for offsets
    m_stencil.define(m_eblg.getDBL(), m_eblg.getDM());
    EBCellFactory fact(m_eblg.getEBISL());
    FabArray<EBCellFAB> dummy(m_eblg.getDBL(),m_eblg.getDM(), nvar, m_ghost, MFInfo(), fact);
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      BL_PROFILE("vof stencil definition");
      const     Box& grid = m_eblg.getDBL()  [mfi];
      const EBISBox& ebis = m_eblg.getEBISL()[mfi];
 
      IntVectSet ivs = ebis.getIrregIVS(grid);
      VoFIterator vofit(ivs, ebis.getEBGraph());
      const Vector<VolIndex>& vofvec = vofit.getVector();
      // cast from VolIndex to BaseIndex
      Vector<std::shared_ptr<BaseIndex> >    dstVoF(vofvec.size());
      Vector<std::shared_ptr<BaseStencil> > stencil(vofvec.size());
      Vector<VoFStencil> allsten(vofvec.size());
      // fill stencils for the vofs
      for(int ivec = 0; ivec < vofvec.size(); ivec++)
      {
        getLocalStencil(allsten[ivec], vofvec[ivec], mfi);
        // another cast from VolIndex to BaseIndex
        BaseIndex* baseind    = (BaseIndex*)  (& vofvec[ivec]);
        BaseStencil* basesten = (BaseStencil*)(&allsten[ivec]);
        dstVoF[ivec]  = std::shared_ptr<BaseIndex  >(baseind , &null_deleter_vof);
        stencil[ivec] = std::shared_ptr<BaseStencil>(basesten, &null_deleter_vsten);
      }
      m_stencil[mfi] = std::shared_ptr<AggStencil<EBCellFAB, EBCellFAB > >
        (new AggStencil<EBCellFAB, EBCellFAB >(dstVoF, stencil, dummy[mfi], dummy[mfi]));

    }
  }
//----------------------------------------------------------------------------
  void
  EBNormalizeByVolumeFraction::
  normalize(FabArray<EBCellFAB>      & a_Q,
            const FabArray<EBCellFAB>& a_kappaQ,
            int a_dst, int a_nco)
  {
    BL_ASSERT(a_Q.nGrow() == m_ghost);
    BL_ASSERT(a_kappaQ.nGrow() == m_ghost);
    BL_PROFILE("EBNormalizer::normalize");

    EBCellFactory fact(m_eblg.getEBISL());
    FabArray<EBCellFAB>& castkap = const_cast<FabArray<EBCellFAB>&>(a_kappaQ);
    castkap.FillBoundary();

    int begin  = a_dst;
    int length = a_nco;

    a_Q.copy(a_kappaQ, a_dst, a_dst, a_nco, 0, 0);
    for(MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      m_stencil[mfi]->apply(a_Q[mfi], a_kappaQ[mfi], begin, begin, length, false);
    }
  }


}

