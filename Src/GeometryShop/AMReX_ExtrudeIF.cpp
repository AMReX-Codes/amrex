#include "AMReX_ExtrudeIF.H"
#include "AMReX_UnionIF.H"

namespace amrex
{
    ExtrudeIF::ExtrudeIF(const BaseIF& a_impFunc1,
            const bool&   a_inside)
    {
        // Make a copy of the implicit function and note only one was given
        m_impFunc1 = a_impFunc1.newImplicitFunction();

        // Save inside flag
        m_inside = a_inside;
    }

    ExtrudeIF::~ExtrudeIF()
    {
        delete m_impFunc1;
    }

    Real ExtrudeIF::value(const RealVect& a_point) const
    {
        Real retval=0.0;

        Real x = a_point[0];
        Real y = a_point[1];
#if AMREX_SPACEDIM > 2
        RealVect coord3(x,y,0.0);
        retval = m_impFunc1->value(coord3);
#else
        RealVect coord2(x,y);
        retval = m_impFunc1->value(coord2);
#endif
        // Change the sign to change inside to outside
        if (!m_inside)
        {
            retval = -retval;
        }
        return retval;
    }

    BaseIF* ExtrudeIF::newImplicitFunction() const
    {
        ExtrudeIF* extrPtr;

        extrPtr = new ExtrudeIF(*m_impFunc1,m_inside);

        return static_cast<BaseIF*>(extrPtr);
    }
}
