
#ifndef AMREX_EB2_IF_SPLINE_H_
#define AMREX_EB2_IF_SPLINE_H_

#include <vector>
#include <AMReX_EB2_IF_Base.H>
#include "AMReX_RealVect.H"
#include "AMReX.H"
#include "AMReX_Vector.H"
#include "AMReX_Array.H"
#include "AMReX_distFcnElement.H"

namespace amrex { namespace EB2 {

/*
 * Implicit function to specify a spline through a set of control points
 */
class SplineIF {
 public:
  //! Constructor
  SplineIF() {}

  //! Destructor
  ~SplineIF() {}

  void addSplineElement(std::vector<amrex::RealVect> pts) {
    SplineDistFcnElement2d * theSpline = new SplineDistFcnElement2d();
    theSpline->set_control_points(pts);
    theSpline->calc_D();
    geomElements.push_back(theSpline);
  }

  void addLineElement(std::vector<amrex::RealVect> pts) {
    LineDistFcnElement2d * theLine = new LineDistFcnElement2d();
    theLine->set_control_points(pts);
    geomElements.push_back(theLine);
  }

  amrex::Real operator() (const amrex::RealArray& p) const {
    amrex::RealVect cp;
    distFcnElement2d * closesetGeomElement;
    amrex::Real dist;
    dist = 1.0e29;
    amrex::RealVect x;
    x[0] = p[0];x[1] = p[1]; x[2] = p[2];
    for (auto * geom : geomElements ) {
      amrex::Real d = geom->cpdist(x, cp);
      if (d < dist) {
        dist = d;
        closesetGeomElement = geom;
      }
    }
    amrex::Real side = closesetGeomElement->cpside(x,cp);
    return dist*side;
  }

  //! private:
  //! The geometry elements used to compute distance function
  amrex::Vector<distFcnElement2d*> geomElements;
};

}}
#endif
