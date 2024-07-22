
#ifndef AMREX_EB2_IF_SPLINE_H_
#define AMREX_EB2_IF_SPLINE_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_RealVect.H>
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_distFcnElement.H>

namespace amrex::EB2 {

/*
 * Implicit function to specify a spline through a set of control points
 */
class SplineIF {
 public:

  void addSplineElement(std::vector<amrex::RealVect> pts) {
    auto* theSpline = new SplineDistFcnElement2d();
    theSpline->set_control_points(std::move(pts));
    theSpline->calc_D();
    geomElements.push_back(theSpline);
  }

  void addLineElement(std::vector<amrex::RealVect> pts) {
    auto* theLine = new LineDistFcnElement2d();
    theLine->set_control_points(std::move(pts));
    geomElements.push_back(theLine);
  }

  [[nodiscard]] amrex::Real operator() (const amrex::RealArray& p) const {
    amrex::RealVect cp;
    distFcnElement2d * closesetGeomElement = nullptr;
    amrex::Real dist;
    dist = 1.0e29_rt;
    amrex::RealVect x(AMREX_D_DECL(p[0], p[1], p[2]));
    for (auto * geom : geomElements ) {
      amrex::Real d = geom->cpdist(x, cp);
      if (d < dist) {
        dist = d;
        closesetGeomElement = geom;
      }
    }
    amrex::Real side = closesetGeomElement ? closesetGeomElement->cpside(x,cp) : Real(0.0);
    return dist*side;
  }

  //! private:
  //! The geometry elements used to compute distance function
  amrex::Vector<distFcnElement2d*> geomElements;
};

}
#endif
