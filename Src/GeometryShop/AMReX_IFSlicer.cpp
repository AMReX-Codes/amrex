
#include "AMReX_IFSlicer.H"


namespace amrex
{

// empty constructor (dim == GLOBALDIM)
  IFSlicer<GLOBALDIM>::IFSlicer()
  {
    m_implicitFunction = NULL;
  }

//copy constructor
  IFSlicer<GLOBALDIM>::IFSlicer(const IFSlicer<GLOBALDIM> & a_IFSlicer)
  {
    m_implicitFunction = a_IFSlicer.m_implicitFunction->newImplicitFunction();
  }

// constructor (dim == GLOBALDIM)
  IFSlicer<GLOBALDIM>::IFSlicer(const BaseIF & a_implicitFunction)
  {
    m_implicitFunction = a_implicitFunction.newImplicitFunction();
  }

// Destructor (dim == GLOBALDIM)
  IFSlicer<GLOBALDIM>::~IFSlicer()
  {
    if (m_implicitFunction != NULL)
    {
      delete m_implicitFunction;
    }
  }

  Real IFSlicer<GLOBALDIM>::value(const IntVect  & a_partialDerivative,
                                  const RealVect & a_point) const
  {
    return m_implicitFunction->value(a_partialDerivative,a_point);
  }

  void IFSlicer<GLOBALDIM>::print(ostream& a_out) const
  {
    amrex::Abort("Not implemented");
    // m_implicitFunction->print(a_out);
  }

}

