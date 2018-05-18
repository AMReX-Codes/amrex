#ifndef AMREX_NOREFINEMENTIMPLEM_H_
#define AMREX_NOREFINEMENTIMPLEM_H_

template <int dim> NoRefinement<dim>::NoRefinement()
{
}

template <int dim> NoRefinement<dim>::~NoRefinement()
{
}

template <int dim> bool NoRefinement<dim>::doRefine(IndexTM<int,dim>          & a_refineInDir,
                                                    const CutCellMoments<dim> & a_ccm,
                                                    const int                 & a_numberOfRefinements)
{
  for (int idir = 0; idir < dim; idir++)
  {
    a_refineInDir[idir] = 0;
  }

  return false;
}


#endif
