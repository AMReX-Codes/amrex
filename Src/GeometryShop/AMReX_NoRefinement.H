
#ifndef AMREX_NOREFINEMENT_H_
#define AMREX_NOREFINEMENT_H_

#include "AMReX_RefinementCriterion.H"



///
/**
   This refinement criterion always false.  Used to stop subdivision.
*/
template <int dim> class NoRefinement: public RefinementCriterion<dim>
{
public:
  /// Null constructor
  NoRefinement();

  /// Destructor
  virtual ~NoRefinement();

  /// Never refine - return false and "a_refineInDir" zero in all components.
  /**
     Never refine - return false and "a_refineInDir" zero in all components.
  */
  virtual bool doRefine(IndexTM<int,dim>          & a_refineInDir,
                        const CutCellMoments<dim> & a_ccm,
                        const int                 & a_numberOfRefinements);

protected:
};



#include "AMReX_NoRefinementImplem.H"

#endif
