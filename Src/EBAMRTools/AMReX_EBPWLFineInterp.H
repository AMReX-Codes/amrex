#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves thurs, aug 23, 2001

#ifndef _EBPWLFINEINTERP_H_
#define _EBPWLFINEINTERP_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "BaseIVFAB.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBPWLFineInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBPWLFineInterp();

  ///
  ~EBPWLFineInterp();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  EBPWLFineInterp(const DisjointBoxLayout& a_dblFine,
                  const DisjointBoxLayout& a_dblCoar,
                  const EBISLayout& a_ebislFine,
                  const EBISLayout& a_ebislCoar,
                  const ProblemDomain& a_domainCoar,
                  const int& a_nref,
                  const int& a_nvar,
                  const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout& a_ebislFine,
              const EBISLayout& a_ebislCoar,
              const ProblemDomain& a_domainCoar,
              const int& a_nref,
              const int& a_nvar,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Interpolate the fine data from the coarse data
     over the intersection
     of the fine layout with the refined coarse layout. \\
     NOTE: This function sends these arrays
     to Fortran so you must set the covered
     cell values to something.  Unitialized memory
     might cause IEEE problems with some compilers.
     {\bf Arguments:}\\
     coarData: The data over the coarse layout.\\
     fineData: The data over the fine layout.  \\
     Fine and coarse data must
     have the same number of variables.\\
     variables:  The variables to interpolate.  Those not
     in this range will be left alone.   This range of variables
     must be in both the coarse and fine data.
  */
  void
  interpolate(LevelData<EBCellFAB>& a_fineData,
              const LevelData<EBCellFAB>& a_coarseData,
              const Interval& a_variables);

protected:
  void
  interpolateFAB(EBCellFAB& a_fine,
                 const EBCellFAB& a_coar,
                 const DataIndex& a_datInd,
                 const Interval& a_variables) const;

  void makeDerivStencils();

  void setDefaultValues();

  bool m_isDefined;

  DisjointBoxLayout m_coarsenedFineGrids;
  ProblemDomain m_coarDomain;
  EBISLayout m_coarsenedFineEBISL;

  int m_refRat;
  int m_nComp;
  LevelData<EBCellFAB> m_coarsenedFineData;

  //irregular cells
  LayoutData<IntVectSet> m_irregRegions;

  //diffs high and low---need to make the choice
  //about which to use based on the data.
  LayoutData<BaseIVFAB<VoFStencil> >  m_hiStencils[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> >  m_loStencils[SpaceDim];

private:
  //disallowed for all the usual reasons
  EBPWLFineInterp(const EBPWLFineInterp& ebcin)
  {
    MayDay::Error("ebpwl 2 invalid operator");
  }
  void operator=(const EBPWLFineInterp& fabin)
  {
    MayDay::Error("ebpwl 3 invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
