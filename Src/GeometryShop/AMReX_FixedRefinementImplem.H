
#ifndef AMREX_FIXEDREFINEMENTIMPLEM_H_
#define AMREX_FIXEDREFINEMENTIMPLEM_H_

namespace amrex
{

  template <int dim> FixedRefinement<dim>::FixedRefinement()
  {
  }

  template <int dim> FixedRefinement<dim>::FixedRefinement(const int & a_maxNumberOfRefinements)
  {
    setMaxNumberOfRefinements(a_maxNumberOfRefinements);
  }

  template <int dim> FixedRefinement<dim>::~FixedRefinement()
  {
  }

  template <int dim> bool FixedRefinement<dim>::doRefine(IndexTM<int,dim>          & a_refineInDir,
                                                         const CutCellMoments<dim> & a_ccm,
                                                         const int                 & a_numberOfRefinements)
  {
    bool retval = false;

    if (a_numberOfRefinements < m_maxNumberOfRefinements)
    {
      retval = true;

      for (int idir = 0; idir < dim; idir++)
      {
        a_refineInDir[idir] = 1;
      }
    }

    return retval;
  }

  template <int dim> void FixedRefinement<dim>::setMaxNumberOfRefinements(const int & a_maxNumberOfRefinements)
  {
    if (a_maxNumberOfRefinements < 0)
    {
      amrex::Abort("FixedRefinement<dim>::setNumberOfRefinements - maxNumberOfRefinements must be >= 0");
    }

    m_maxNumberOfRefinements = a_maxNumberOfRefinements;
  }

  template <int dim> int FixedRefinement<dim>::getMaxNumberOfRefinements()
  {
    return m_maxNumberOfRefinements;
  }

}

#endif
