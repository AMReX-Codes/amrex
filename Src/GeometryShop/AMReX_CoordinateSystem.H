#ifndef _COORDINATESYSTEM_H_
#define _COORDINATESYSTEM_H_


#include <map>
using std::map;
#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_IndexTM.H"
#include "AMReX_Notation.H"
#include "AMReX_IFSlicer.H"

using std::endl;
using std::cerr;

namespace amrex
{

  template <int dim> class CoordinateSystem
  {
  public:
    // data types for vector of int and Real in dimension "dim"
    typedef IndexTM<int,dim>  IvDim;
    typedef IndexTM<Real,dim> RvDim;

    // member data
    RvDim m_origin;
    RvDim m_dx;

    // empty constructor
    CoordinateSystem();

    // copy constructor
    CoordinateSystem(const CoordinateSystem<dim>& a_info);

    // makes an CoordinateSystem from a origin and dx
    CoordinateSystem(const RvDim  & a_origin,
                     const RvDim  & a_dx);

    // makes a reduced CoordinateSystem
    CoordinateSystem(const CoordinateSystem<dim+1> & a_hISystem,
                     const int                     & a_idir);

    // destructor
    ~CoordinateSystem();

    // convert "a_point" in coordinate system "a_system" to a point in the
    // current coordinate system.
    RvDim convert(const RvDim                 & a_point,
                  const CoordinateSystem<dim> & a_system) const;

    // convert "a_coord" in direction "a_dir" in coordinate system "a_system"
    // to a value in "a_dir" in the current coordinate system.
    Real convertDir(const Real                  & a_coord,
                    const CoordinateSystem<dim> & a_system,
                    const int                   & a_dir) const;

    // output functions
    void print(ostream& out) const;

    // operators
    void operator=(const CoordinateSystem<dim> & a_coordinateSystem);
  };

}

#include "AMReX_CoordinateSystemImplem.H"

#endif
