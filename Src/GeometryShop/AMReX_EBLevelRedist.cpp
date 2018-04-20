#include "AMReX_EBLevelRedist.H"
#include "AMReX_BaseIVFactory.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_EBIndexSpace.H"

namespace amrex
{
  /***********************/
  EBLevelRedist::EBLevelRedist()
  {
    m_isDefined = false;
  }
  /***********************/
  EBLevelRedist::~EBLevelRedist()
  {
  }
  /***********************/
  EBLevelRedist::
  EBLevelRedist(const EBLevelGrid & a_eblg,
                const int         & a_ncomp,
                int                 a_redistRad)
  {
    define(a_eblg, a_ncomp, a_redistRad);
  }
  /***********************/
  /***********************/
  void
  EBLevelRedist::
  resetWeights(const FabArray<EBCellFAB>& a_modifier,
               const int& a_ivar)
  {
    m_stencil.resetWeights(a_modifier, a_ivar);
  }
  /***********************/
  /***********************/
  void
  EBLevelRedist::
  define(const EBLevelGrid & a_eblg,
         const int         & a_ncomp,
         int                 a_redistRad)
  {
    BL_PROFILE("EBLevelRedist::define");
    m_isDefined = true;
    m_eblg  = a_eblg;
    m_ncomp     = a_ncomp;
    m_redistRad = a_redistRad;
               
    //this sets the stencil to volume-weighted.
    //use resetWeights to set to mass weighted or whatever
    m_stencil.define(m_eblg, m_redistRad);
               
    //define the sets for iterating over
    //to be the irregular cells over the
    //grid grown by the redistribution radius
    m_sets   = shared_ptr<LayoutData<IntVectSet> >
      (new LayoutData<IntVectSet>(m_eblg.getDBL(), m_eblg.getDM()));
    m_vofit  = shared_ptr<LayoutData<VoFIterator> >
      (new LayoutData<VoFIterator>(m_eblg.getDBL(), m_eblg.getDM()));

    for (MFIter mfi(m_eblg.getDBL(), m_eblg.getDM()); mfi.isValid(); ++mfi)
    {
      Box thisBox = m_eblg.getDBL()[mfi];
      thisBox.grow(m_redistRad);
      thisBox &= m_eblg.getDomain();
      (*m_sets) [mfi] = m_eblg.getEBISL()[mfi].getIrregIVS(thisBox);
      (*m_vofit)[mfi].define((*m_sets)[mfi], m_eblg.getEBISL()[mfi].getEBGraph());
    }
    BaseIVFactory<Real> factory(m_eblg.getEBISL(), m_sets);


    m_buffer.define(m_eblg.getDBL(), m_eblg.getDM(), m_ncomp, m_redistRad,  MFInfo(), factory);
    setToZero();
  }
  /***********************/
  /***********************/
  void
  EBLevelRedist::setToZero()
  {
    BL_PROFILE("EBLevelRedist::setToZero");
    BL_ASSERT(isDefined());
    for (MFIter mfi(m_buffer); mfi.isValid(); ++mfi)
    {
      m_buffer[mfi].setVal(0.0);
    }
  }
  /***********************/
  /***********************/
  void
  EBLevelRedist::
  increment(const BaseIVFAB<Real>& a_massDiff,
            const MFIter         & a_datInd,
            int idst, int inco)
  {
    BL_PROFILE("EBLevelRedist::increment");
    BL_ASSERT(isDefined());
    BaseIVFAB<Real>& bufFAB = m_buffer[a_datInd];
    const IntVectSet& fabIVS = a_massDiff.getIVS();
    const IntVectSet& bufIVS = (*m_sets)[a_datInd];
               
    IntVectSet ivs = m_eblg.getEBISL()[a_datInd].getIrregIVS(m_eblg.getDBL()[a_datInd]);;
    BL_ASSERT(fabIVS.contains(ivs));
    BL_ASSERT(bufIVS.contains(ivs));
    VoFIterator& vofit = (*m_vofit)[a_datInd];
    for (vofit.reset(); vofit.ok(); ++vofit)
    {
      for (int ivar = idst; ivar <  (idst + inco); ivar++)
      {
        bufFAB(vofit(), ivar) += a_massDiff(vofit(), ivar);
      }
    }
  }
  /***********************/
  void
  EBLevelRedist::
  redistribute(FabArray<EBCellFAB>& a_solution,
               int idst, int inco)
  {
    redistribute(a_solution, idst, idst, inco);
  }
               
  /***********************/
  void
  EBLevelRedist::
  redistribute(FabArray<EBCellFAB>& a_solution,
               int a_isrc, int a_idst, int a_inco)

  {
    BL_PROFILE("EBLevelRedist::redistribute");
    BL_ASSERT( a_isrc >= 0);
    BL_ASSERT( a_idst >= 0);
    BL_ASSERT((a_isrc + a_inco -1) <  m_ncomp);
    BL_ASSERT((a_idst + a_inco -1) <  a_solution.nComp());
               
    BL_ASSERT(isDefined());
               
    //exchange ghost cell information of the buffer.
    //this way the redistribution from ghost cells will
    //account for fine-fine interfaces.
    m_buffer.FillBoundary();

    //loop over grids.
    for(MFIter mfi(m_buffer); mfi.isValid(); ++mfi)
    {
      const BaseIVFAB<Real>      & bufFAB  =  m_buffer[mfi];
      const BaseIVFAB<VoFStencil>& stenFAB = m_stencil[mfi];
      const Box& grid = m_eblg.getDBL()[mfi];
      EBCellFAB& solFAB = a_solution[mfi];
               
      //loop over the irregular vofs in the grown box.
      EBGraph graph = m_eblg.getEBISL()[mfi].getEBGraph();
      VoFIterator& vofit = (*m_vofit)[mfi];
      for (vofit.reset(); vofit.ok(); ++vofit)
      {
        const VolIndex& srcVoF = vofit();
        const VoFStencil& vofsten = stenFAB(srcVoF, 0);
        for (int isten = 0; isten < vofsten.size(); isten++)
        {
          const Real    & weight = vofsten.weight(isten);
          const VolIndex& dstVoF = vofsten.vof(isten);
          //only add the mass to the solution if it is in the
          //solutions valid region --- because the buffer is over
          //the grown box, this prevents redistribution off
          //the end of the world.
          if (grid.contains(dstVoF.gridIndex()))
          {
            for (int ivar = 0; ivar < a_inco; ivar++)
            {
              int isrc = ivar + a_isrc;
              int idst = ivar + a_idst;
               
              const Real& mass = bufFAB(srcVoF, isrc);
              const Real& solu = solFAB(dstVoF, idst);
              solFAB(dstVoF, idst) = mass*weight + solu;
            }
          } //end check if adding to valid soltuion
        } //end loop over vof stencil
      } //end loop over irregular vofs in grown box
    } //end loop over grids.
  } //if only I got paid by the curly brace
  /***********************/
  bool
  EBLevelRedist::isDefined() const
  {
    return m_isDefined;
  }
  /***********************/
}
