#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri aug 17, 2001

#ifndef _EBCOARSEAVERAGE_H_
#define _EBCOARSEAVERAGE_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBLevelGrid.H"
#include "EBCellFAB.H"
#include "EBFluxFAB.H"
#include "BaseIVFAB.H"
#include "BaseIVFactory.H"
#include "EBFaceFAB.H"
#include "Interval.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

/// replaces coarse level data with an average of fine level data.
/**
   This class replaces data at a coarse level of refinement with an
   average of data at a finer level of refinement, in areas where fine
   data is present.  Coarse level data is not modified where fine
   level data is not present.
  */
class EBCoarseAverage
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBCoarseAverage();

  ///
  ~EBCoarseAverage();

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
  EBCoarseAverage(const DisjointBoxLayout& dblFine,
                  const DisjointBoxLayout& dblCoar,
                  const EBISLayout& ebislFine,
                  const EBISLayout& ebislCoar,
                  const ProblemDomain& domainCoar,
                  const int& nref,
                  const int& nvar,
                  const EBIndexSpace* ebisPtr);

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     eblgFine, eblgCoar: The fine and coarse EBLevelGrids \\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  EBCoarseAverage(const EBLevelGrid& a_eblgFine,
                  const EBLevelGrid& a_eblgCoar,
                  const EBLevelGrid& a_eblgCoFi,
                  const int& nref,
                  const int& nvar);
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
  void define(const DisjointBoxLayout& dblFine,
              const DisjointBoxLayout& dblCoar,
              const EBISLayout& ebislFine,
              const EBISLayout& ebislCoar,
              const ProblemDomain& domainCoar,
              const int& nref,
              const int& nvar,
              const EBIndexSpace* ebisPtr);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     eblgFine, eblgCoar: The fine and coarse EBLevelGrids \\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  void define(const EBLevelGrid& a_eblgFine,
              const EBLevelGrid& a_eblgCoar,
              const EBLevelGrid& a_eblgCoFi,
              const int& nref,
              const int& nvar);
  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Replaces a_coarse_data with the average of a_fine_data within the
     coarsening of the domain of the fine level.  Elsewhere, the
     coarse data is unchanged.  It is an error to call if not
     this->isDefined(). The domain of a_fine_data should be
     the same as the fine domain specified in the most recent call to
     define().  It is expected that the coarse and fine level's
     domains are properly nested.  Both a_coarse_data and a_fine_data
     should have the same number of components specified in the most
     recent call to define().

     {\bf Arguments:}\\
     coarData: The data over the coarse layout.\\
     fineData: The data over the fine layout.  \\
     Fine and coarse data must
     have the same number of variables.\\
     variables:  The variables to average.  Those not
     in this range will be left alone.   This range of variables
     must be in both the coarse and fine data.
  */
  void
  average(LevelData<EBCellFAB>&       a_coarData,
          const LevelData<EBCellFAB>& a_fineData,
          const Interval&             a_variables);

  ///this does area weighted face averaging
  void
  average(LevelData<EBFluxFAB>&       a_coarData,
          const LevelData<EBFluxFAB>& a_fineData,
          const Interval&             a_variables);

  ///this does area-weighted averaging for the irregular area (for data on the EB)
  void
  average(LevelData<BaseIVFAB<Real> >&       a_coarData,
          const LevelData<BaseIVFAB<Real> >& a_fineData,
          const Interval&                    a_variables);

protected:
  void
  averageFAB(EBCellFAB&       a_coar,
             const EBCellFAB& a_fine,
             const DataIndex& a_datInd,
             const Interval&  a_fineInterv,
             const Interval&  a_coarInterv) const;

  void
  averageFAB(EBFaceFAB&       a_coar,
             const EBFaceFAB& a_fine,
             const DataIndex& a_datInd,
             const Interval&  a_fineInterv,
             const Interval&  a_coarInterv,
             const int&       a_dir) const;

  void
  averageFAB(BaseIVFAB<Real>&       a_coar,
             const BaseIVFAB<Real>& a_fine,
             const DataIndex&       a_datInd,
             const Interval&        a_variables) const;


  void
  setDefaultValues();

  bool m_isDefined;
  //true only if the layout is not coarsenable
  //and the grid covers the domain.   This
  //happens in multigrid sometimes
  bool m_useFineBuffer;
  EBLevelGrid m_eblgCoar;
  EBLevelGrid m_eblgFine;
  EBLevelGrid m_eblgCoFi;
  LayoutData<IntVectSet> m_irregSetsCoFi;

  //only defined if using fine buffer
  LayoutData<IntVectSet> m_irregSetsFine;

  int m_refRat;
  int m_nComp;

private:
    //disallowed for all the usual reasons
  EBCoarseAverage(const EBCoarseAverage& ebcin)
  {
    MayDay::Error("ebca 2 invalid operator");
  }
  void operator=(const EBCoarseAverage& fabin)
  {
    MayDay::Error("ebca 3 invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
