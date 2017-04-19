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

#ifndef _EBLEVELREDIST_H_
#define _EBLEVELREDIST_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "BaseIVFAB.H"
#include "RedistStencil.H"
#include "NamespaceHeader.H"

///
/**
  */
class EBLevelRedist
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBLevelRedist();

  ///
  ~EBLevelRedist();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().
     ebisl needs ghost cells out to 3*redistradius of the stencil.
     This sets the stencil to volume-weighted.
     Use resetWeights to set to mass weighted or whatever.
  */
  EBLevelRedist(const DisjointBoxLayout& a_dbl,
                const EBISLayout& a_ebisl,
                const ProblemDomain& a_domain,
                const int& a_ncomp,
                int a_redistRad = 1,
                bool a_do2DStencil = false);

  ///
  /**
     Defines this object.  Existing information is overriden.
     ebisl needs ghost cells out to 3*redistradius
  */
  void define(const DisjointBoxLayout& a_dbl,
              const EBISLayout& a_ebisl,
              const ProblemDomain& a_domain,
              const int& a_ncomp,
              int a_redistRad = 1,
              bool a_do2DStencil = false);

  ///
  /**
     Modify the weights in the stencil by multiplying by
     the inputs in a normalized way.
     If you want mass
     weighting, send in the density.
  */
  void resetWeights(const LevelData<EBCellFAB>& modifier,
                    const int& ivar);

  ///
  /**
     Store the input mass difference in the buffers of the class
     by INCREMENTING the buffer. \\
     dm += massDiff
   */
  void increment(const BaseIVFAB<Real>& massDiff,
                 const DataIndex& datInd,
                 const Interval& variables);

  ///
  /**
     Redistribute the data contained in the internal buffers. \\
     solution += weight*dm
  */
  void redistribute(LevelData<EBCellFAB>& solution,
                    const Interval& variables);

  ///
  /**
     version where solution is not necessarily ncomp long.
     intervals have to be the same length.
     srcVar can be up to m_ncomp long.
     dstVar can be up to solution.nComp().
   */
  void redistribute(LevelData<EBCellFAB>& a_solution,
                    const Interval& a_srcVar,
                    const Interval& a_dstVar);

  ///
  /**
     Adjusts an explicit Laplacian calculation so that:
     a_solution_i =
       kappa_i*lap_i + (1-kappa_i)*sum_s(kappa_s*lap_s)/(sum_s kappa_s)
     where s is the set of cells neighboring, but not including i
     NOTE: a_solution must be passed in equal to a_kappaLap
  */
  void fixExplicitLap(const LevelData<EBCellFAB>& a_kappaLap,
                      LevelData<EBCellFAB>& a_solution,
                      const Interval& variables);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Reset the mass buffers of the class.
   */
  void setToZero();

protected:
  int redistRad;
  RedistStencil m_stencil;
  int m_ncomp;
  int m_redistRad;
  DisjointBoxLayout  m_grids;
  ProblemDomain m_domain;
  EBISLayout m_ebisl;
  bool m_isDefined;
  LevelData<BaseIVFAB<Real> > m_buffer;
  LayoutData<IntVectSet> m_sets;
private:

  //forbidden for all the usual reasons
  void operator=(const EBLevelRedist& a_stenin)
  {
    MayDay::Error("Invalid operator");
  }

  EBLevelRedist(const EBLevelRedist& a_stenin)
  {
    MayDay::Error("Invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
