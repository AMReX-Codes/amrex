
#ifndef _REFINEMENTCRITERION_H_
#define _REFINEMENTCRITERION_H_

#include "AMReX_IndexTM.H"
#include "AMReX_CutCellMoments.H"


///
/**
   This is the base class for refinement criterion used for the subdivision of
   cells in geometry generation.  If the result if "doRefine()" is true, the
   cell is subdivided in the directions in which "a_refineInDir" is non-zero.
   The results of this subdivision are then combined to get results for the
   original cell.
*/
template <int dim> class RefinementCriterion
{
public:

  // empty constructor
  RefinementCriterion()
    {
      m_constraintsSucceeded       = true;
      m_baseMaxNumberOfRefinements = 1;
    }

  // constructor
  RefinementCriterion(const int& a_baseMaxNumberOfRefinements)
    {
      m_constraintsSucceeded   = true;
      m_baseMaxNumberOfRefinements = a_baseMaxNumberOfRefinements;
    }

  // copy constructor
  RefinementCriterion(const RefinementCriterion<dim>& a_RefinementCriterion)
    :m_constraintsSucceeded  (a_RefinementCriterion.m_constraintsSucceeded),
     m_baseMaxNumberOfRefinements(a_RefinementCriterion.m_baseMaxNumberOfRefinements)
    {
    }

  // Destructor
  ~RefinementCriterion()
    {
    }

  /// Should a cell be subdivided and in which directions
  /**
     This method returns true if the current cell should be subdivided.  The
     subdivsion should occur in all the directions where "a_refineInDir" is
     non-zero.
  */
  virtual bool baseDoRefine(IndexTM<int,dim>          & a_refineInDir,
                            const CutCellMoments<dim> & a_ccm,
                            const int                 & a_numberOfRefinements)
    {
      //false = don't refine
      bool baseRetval = false;

      //refine in no directions
      a_refineInDir = IndexTM<int,dim>::Zero;

      //check whether a priori limit has been reached
      bool exceededMaxNumber = false;
      if (a_numberOfRefinements >= m_baseMaxNumberOfRefinements)
      {
        exceededMaxNumber = true;
      }

      if (!exceededMaxNumber)
      {
        //check whether normal equals the zero vector or whether constraints were active
        if (a_ccm.m_badNormal || !m_constraintsSucceeded)
        {
          baseRetval = true;
          a_refineInDir = IndexTM<int,dim>::Unit;
        }
      }

      bool derivedRetval = doRefine(a_refineInDir,
                                    a_ccm,
                                    a_numberOfRefinements);


      bool retval = baseRetval || derivedRetval;

      //true = refine
      return retval;
    }

  virtual bool doRefine(IndexTM<int,dim>          & a_refineInDir,
                        const CutCellMoments<dim> & a_ccm,
                        const int                 & a_numberOfRefinements)
    {
      //empty implementation, which is useful if the refinement criterion is just the base class
      return false;
    }

  //Records "lsCode" from the least squares calculation
  void setConstrantSuccessStatus(const bool& a_status)
    {
      m_constraintsSucceeded = a_status;
    }

  //Retrieves "lsCode" from the least squares calculation
  bool getConstrantSuccessStatus()
    {
      return m_constraintsSucceeded;
    }

  //set max number of refinements
  void setBaseMaxNumberOfRefinements(const int & a_baseMaxNumberOfRefinements)
    {
      if (a_baseMaxNumberOfRefinements < 0)
      {
        amrex::Abort("FixedRefinement<dim>::setNumberOfRefinements - maxNumberOfRefinements must be >= 0");
      }

      m_baseMaxNumberOfRefinements = a_baseMaxNumberOfRefinements;
    }

  //get max number of refinements
  int getBaseMaxNumberOfRefinements()
    {
      return m_baseMaxNumberOfRefinements;
    }

  void print(ostream& a_out) const
    {
      a_out << "m_constraintsSucceeded  = " << m_constraintsSucceeded << "\n";
    }

  // equals operator
  void operator=(const RefinementCriterion & a_RefinementCriterion)
    {
      m_constraintsSucceeded   = a_RefinementCriterion.m_constraintsSucceeded;
      m_baseMaxNumberOfRefinements = a_RefinementCriterion.m_baseMaxNumberOfRefinements;
    }

protected:
  bool m_constraintsSucceeded;
  int  m_baseMaxNumberOfRefinements;
};

template<int dim> ostream& operator<<(ostream                        & a_out,
                                      const RefinementCriterion<dim> & a_RefinementCriterion)
{
  a_RefinementCriterion.print(a_out);
  return a_out;
}


#endif
