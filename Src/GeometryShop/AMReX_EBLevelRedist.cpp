#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves mon oct 15 2001
#include "EBLevelRedist.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBIndexSpace.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
/***********************/
/***********************/
EBLevelRedist::EBLevelRedist()
{
  m_isDefined = false;
}
/***********************/
/***********************/
EBLevelRedist::~EBLevelRedist()
{
}
/***********************/
/***********************/
EBLevelRedist::
EBLevelRedist(const DisjointBoxLayout& a_dbl,
              const EBISLayout& a_ebisl,
              const ProblemDomain& a_domain,
              const int& a_ncomp,
              int  a_redistRad,
              bool a_do2DStencil)
{
  define(a_dbl, a_ebisl, a_domain, a_ncomp, a_redistRad, a_do2DStencil);
}
/***********************/
/***********************/
void
EBLevelRedist::
resetWeights(const LevelData<EBCellFAB>& a_modifier,
             const int& a_ivar)
{
  m_stencil.resetWeights(a_modifier, a_ivar);
}
/***********************/
/***********************/
void
EBLevelRedist::
define(const DisjointBoxLayout& a_dbl,
       const EBISLayout& a_ebisl,
       const ProblemDomain& a_domain,
       const int& a_ncomp,
       int a_redistRad,
       bool a_do2DStencil)
{
  CH_TIME("EBLevelRedist::define");
  m_isDefined = true;

  m_grids     = a_dbl;
  m_ebisl     = a_ebisl;
  m_domain    = a_domain;
  m_ncomp     = a_ncomp;
  m_redistRad = a_redistRad;

  //this sets the stencil to volume-weighted.
  //use resetWeights to set to mass weighted or whatever
  m_stencil.define(m_grids, m_ebisl, m_domain, m_redistRad, a_do2DStencil);

  //define the sets for iterating over
  //to be the irregular cells over the
  //grid grown by the redistribution radius
  m_sets.define(m_grids);
  int redistRad = m_stencil.getRedistRadius();
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_grids.get(dit());
      thisBox.grow(redistRad);
      thisBox &= m_domain;
      m_sets[dit()] = m_ebisl[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_ebisl, m_sets);
  IntVect ivghost = redistRad*IntVect::Unit;
  m_buffer.define(m_grids, m_ncomp, ivghost, factory);
  setToZero();
}
/***********************/
/***********************/
void
EBLevelRedist::setToZero()
{
  CH_TIME("EBLevelRedist::setToZero");
  CH_assert(isDefined());
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      m_buffer[dit()].setVal(0.0);
    }
}
/***********************/
/***********************/
void
EBLevelRedist::
increment(const BaseIVFAB<Real>& a_massDiff,
          const DataIndex& a_datInd,
          const Interval& a_variables)
{
  CH_TIME("EBLevelRedist::increment");
  CH_assert(isDefined());
  BaseIVFAB<Real>& bufFAB = m_buffer[a_datInd];
  const IntVectSet& fabIVS = a_massDiff.getIVS();
  const IntVectSet& bufIVS = m_sets[a_datInd];

  IntVectSet ivs = m_ebisl[a_datInd].getIrregIVS(m_grids.get(a_datInd));;
  CH_assert(fabIVS.contains(ivs));
  CH_assert(bufIVS.contains(ivs));
  for (VoFIterator vofit(ivs, m_ebisl[a_datInd].getEBGraph());
      vofit.ok(); ++vofit)
    {
      for (int ivar = a_variables.begin();
          ivar <= a_variables.end(); ivar++)
        {
          bufFAB(vofit(), ivar) += a_massDiff(vofit(), ivar);
        }
    }
}
/***********************/
void
EBLevelRedist::
redistribute(LevelData<EBCellFAB>& a_solution,
             const Interval& a_variables)
{
  redistribute(a_solution, a_variables, a_variables);
}

/***********************/
void
EBLevelRedist::
redistribute(LevelData<EBCellFAB>& a_solution,
             const Interval& a_srcVar,
             const Interval& a_dstVar)
{
  CH_TIME("EBLevelRedist::redistribute");
  CH_assert(a_srcVar.size() == a_dstVar.size());
  CH_assert(a_srcVar.begin() >= 0);
  CH_assert(a_dstVar.begin() >= 0);
  CH_assert(a_srcVar.end() <  m_ncomp);
  CH_assert(a_dstVar.end() <  a_solution.nComp());

  CH_assert(isDefined());

  //pout() << "redistribute 0" << endl;
  //exchange ghost cell information of the buffer.
  //this way the redistribution from ghost cells will
  //account for fine-fine interfaces.
  Interval wholeInterv(0, m_ncomp-1);
  m_buffer.exchange(wholeInterv);
  //pout() << "redistribute 1" << endl;
  //loop over grids.
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>& bufFAB = m_buffer[dit()];
      const BaseIVFAB<VoFStencil>& stenFAB = m_stencil[dit()];
      const Box& grid = m_grids.get(dit());
      EBCellFAB& solFAB = a_solution[dit()];

      //pout() << "redistribute 2" << endl;
      //loop over the irregular vofs in the grown box.
      for (VoFIterator vofit(m_sets[dit()], m_ebisl[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& srcVoF = vofit();
          const VoFStencil& vofsten = stenFAB(srcVoF, 0);
          for (int isten = 0; isten < vofsten.size(); isten++)
            {
              const Real& weight = vofsten.weight(isten);
              const VolIndex& dstVoF = vofsten.vof(isten);
              //only add the mass to the solution if it is in the
              //solutions valid region --- because the buffer is over
              //the grown box, this prevents redistribution off
              //the end of the world.
              if (grid.contains(dstVoF.gridIndex()))
                {
                  for (int ivar = 0; ivar < a_srcVar.size(); ivar++)
                    {
                      int isrc = ivar + a_srcVar.begin();
                      int idst = ivar + a_dstVar.begin();

                      const Real& mass = bufFAB(srcVoF, isrc);
                      const Real& solu = solFAB(dstVoF, idst);
                      solFAB(dstVoF, idst) = mass*weight + solu;
                    }
                } //end check if adding to valid soltuion
            } //end loop over vof stencil
        } //end loop over irregular vofs in grown box
    } //end loop over grids.
}
/***********************/
/***********************/
void
EBLevelRedist::
fixExplicitLap(const LevelData<EBCellFAB>& a_kappaLap,
               LevelData<EBCellFAB>& a_solution,
               const Interval& a_variables)
{
  CH_TIME("EBLevelRedist::fixExplicitLap");
  CH_assert(isDefined());

  //loop over grids.
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      const BaseIVFAB<VoFStencil>& stenFAB = m_stencil[dit()];
      const Box& grid = m_grids.get(dit());
      EBCellFAB& solFAB = a_solution[dit()];
      const EBCellFAB& kappaLapFAB = a_kappaLap[dit()];

      //loop over the irregular vofs in the grown box.
      for (VoFIterator vofit(m_sets[dit()], m_ebisl[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const EBISBox& ebisBox = m_ebisl[dit()];
          const VolIndex& solnVoF = vofit();
          const VoFStencil& vofsten = stenFAB(solnVoF, 0);
          const Real& kappa = ebisBox.volFrac(solnVoF);
          // we start the isten loop at 1 since we do not want to include
          // the solution vof
          for (int isten = 1; isten < vofsten.size(); isten++)
            {
              Real weight = vofsten.weight(isten);
              const VolIndex& neighborVoF = vofsten.vof(isten);
              // weight is kappa_i / sum of kappa_i+s, s<=1
              // we want weight to be simply 1 / sum of kappa_i+s, where s=1
              //const Real& kappalocal = ebisBox.volFrac(neighborVoF);
              weight = 1.0/weight;
              weight -= kappa;
              weight = 1.0/weight;
              if (grid.contains(neighborVoF.gridIndex()))
                {
                  for (int ivar = a_variables.begin();
                      ivar <= a_variables.end(); ivar++)
                    {
                      solFAB(solnVoF,ivar) = solFAB(solnVoF,ivar) +
                        (1.0-kappa)*weight*kappaLapFAB(neighborVoF,ivar);
                    }
                }
            } //end loop over vof stencil
        } //end loop over irregular vofs in grown box
    } //end loop over grids.
}
/***********************/
/***********************/
bool
EBLevelRedist::isDefined() const
{
  return m_isDefined;
}
/***********************/
/***********************/
#include "NamespaceFooter.H"
