
#include "AMReX_IFSlicer.H"


// empty constructor (dim == BL_SPACEDIM)
IFSlicer<BL_SPACEDIM>::IFSlicer()
{
  m_implicitFunction = NULL;
}

//copy constructor
IFSlicer<BL_SPACEDIM>::IFSlicer(const IFSlicer<BL_SPACEDIM> & a_IFSlicer)
{
  m_implicitFunction = a_IFSlicer.m_implicitFunction->newImplicitFunction();
}

// constructor (dim == BL_SPACEDIM)
IFSlicer<BL_SPACEDIM>::IFSlicer(const BaseIF & a_implicitFunction)
{
  m_implicitFunction = a_implicitFunction.newImplicitFunction();
}

// Destructor (dim == BL_SPACEDIM)
IFSlicer<BL_SPACEDIM>::~IFSlicer()
{
  if (m_implicitFunction != NULL)
  {
    delete m_implicitFunction;
  }
}

Real IFSlicer<BL_SPACEDIM>::value(const IntVect  & a_partialDerivative,
                                  const RealVect & a_point) const
{
  return m_implicitFunction->value(a_partialDerivative,a_point);
}

void IFSlicer<BL_SPACEDIM>::print(ostream& a_out) const
{
  amrex::Abort("Not implemented");
  // m_implicitFunction->print(a_out);
}


