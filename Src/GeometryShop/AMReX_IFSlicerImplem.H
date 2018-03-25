
#ifndef _IFSLICERIMPLEM_H_
#define _IFSLICERIMPLEM_H_

#include <iostream>


// empty constructor
template <int dim> IFSlicer<dim>::IFSlicer()
{
  //set to undefined values
  m_IFSlicer     = NULL;
  m_fixedComp  = LARGEINTVAL;
  m_fixedValue = LARGEREALVAL;
}

//copy constructor
template <int dim >IFSlicer<dim>::IFSlicer(const IFSlicer<dim> & a_IFSlicer)
{
  setParams(a_IFSlicer.m_IFSlicer,
            a_IFSlicer.m_fixedComp,
            a_IFSlicer.m_fixedValue);
}

template <int dim> IFSlicer<dim>::IFSlicer(IFSlicer<dim+1> const * a_IFSlicer,
                                           const int             & a_fixedComp,
                                           const Real            & a_fixedValue)
{
  setParams(a_IFSlicer,
            a_fixedComp,
            a_fixedValue);
}

template <int dim> void IFSlicer<dim>::setParams(IFSlicer<dim+1> const * a_IFSlicer,
                                                 const int             & a_fixedComp,
                                                 const Real            & a_fixedValue)
{
  // copy the inputs
  m_IFSlicer   = a_IFSlicer;
  m_fixedComp  = a_fixedComp;
  m_fixedValue = a_fixedValue;
}

// Destructor
template <int dim> IFSlicer<dim>::~IFSlicer()
{
  //lower dimension doesn't own the objects referred to by any pointer.
}

template<int dim >Real IFSlicer<dim>::value(const IndexTM<int, dim> & a_partialDerivative,
                                            const IndexTM<Real,dim> & a_point) const
{
  Real retval = LARGEREALVAL;

  //use fixedComp and fixedValue to pull partialDerivative and point into one higher dimension
  IndexTM<int, dim + 1> partialDerivative;
  IndexTM<Real,dim + 1> point;
  for (int idir = 0; idir < dim +1; ++idir)
  {
    if (idir < m_fixedComp)
    {
      partialDerivative[idir] = a_partialDerivative[idir];
      point            [idir] = a_point            [idir];
    }
    else if (idir > m_fixedComp)
    {
      partialDerivative[idir] = a_partialDerivative[idir - 1];
      point            [idir] = a_point            [idir - 1];
    }
    else
    {
      partialDerivative[m_fixedComp] = 0;
      point            [m_fixedComp] = m_fixedValue;
    }
  }

  //evaluate at one higher dimension
  retval = m_IFSlicer->value(partialDerivative,point);

  return retval;
}

template <int dim> void IFSlicer<dim>::print(ostream& a_out) const
{
  a_out << "fixed component  = " << m_fixedComp  << "\n";
  a_out << "fixed value      = " << m_fixedValue << "\n";
  a_out << "At one dimension higher: " << dim + 1 << "\n";
  this->m_IFSlicer.print(a_out);
}

template <int dim> ostream& operator<<(ostream             & a_out,
                                       const IFSlicer<dim> & a_IFSlicer)
{
  a_IFSlicer.print(a_out);
  return a_out;
}



#endif
