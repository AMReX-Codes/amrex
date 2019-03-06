#ifndef AMREX_DISTFCNELEMENT_H_
#define AMREX_DISTFCNELEMENT_H_

#include <vector>
#include "AMReX_RealVect.H"
#include "AMReX.H"
#include "AMReX_Vector.H"

namespace amrex {

class distFcnElement2d {
 public:
  //! Constructor
  distFcnElement2d() {}
  ~distFcnElement2d() {}

  virtual distFcnElement2d* newDistFcnElement2d() const = 0;

  virtual  amrex::Real cpdist(amrex::RealVect pt, amrex::RealVect& cp) const = 0;
  virtual  amrex::Real cpside(amrex::RealVect pt, amrex::RealVect& cp) const = 0;
  int solve_thomas(std::vector<amrex::Real> diagminus,
                   std::vector<amrex::Real> diag,
                   std::vector<amrex::Real> diagplus,
                   std::vector<amrex::Real> rhs,
                   std::vector<amrex::Real> X);
};


class LineDistFcnElement2d: public distFcnElement2d {
 public:
  LineDistFcnElement2d() {}
  ~LineDistFcnElement2d() {}

  virtual distFcnElement2d* newDistFcnElement2d() const override;

  void set_control_points(std::vector<amrex::RealVect>  pts);

  virtual  amrex::Real cpdist(amrex::RealVect pt, amrex::RealVect& cp) const override;
  virtual  amrex::Real cpside(amrex::RealVect pt, amrex::RealVect& cp) const override;

  void print_control_points();

 protected:
  void single_seg_cpdist(amrex::RealVect pt,
                         amrex::Real x0, amrex::Real x1,
                         amrex::Real y0, amrex::Real y1,
                         amrex::RealVect& cp,
                         amrex::Real& dist) const;

 private:
  std::vector<amrex::Real> control_points_x;
  std::vector<amrex::Real> control_points_y;


};


class SplineDistFcnElement2d: public distFcnElement2d {
 public:
  SplineDistFcnElement2d() {}
  ~SplineDistFcnElement2d() {}

  virtual distFcnElement2d* newDistFcnElement2d() const override;


  void set_control_points(std::vector<amrex::RealVect>  pts);
  void set_bc_points(amrex::RealVect start, amrex::RealVect end);

  void print_control_points() const;
  void print_spline() const;

  void calc_D(bool clamped_bc = true);

  virtual  amrex::Real cpdist(amrex::RealVect pt, amrex::RealVect& cp) const override;
  virtual  amrex::Real cpside(amrex::RealVect pt, amrex::RealVect& cp) const override;

 protected:
  amrex::Real eval(amrex::Real t, amrex::Real y0, amrex::Real y1,
                   amrex::Real D0, amrex::Real D1) const;
  void  dxbydt(amrex::Real t, amrex::Real y0, amrex::Real y1,
               amrex::Real D0, amrex::Real D1, amrex::Real& dyf,
               amrex::Real& d2yf) const;

  void single_spline_cpdist(amrex::RealVect pt,
                            amrex::Real x0, amrex::Real x1,
                            amrex::Real Dx0, amrex::Real Dx1,
                            amrex::Real y0, amrex::Real y1,
                            amrex::Real Dy0, amrex::Real Dy1,
                            amrex::Real& t, amrex::RealVect& cp,
                            amrex::Real& dist) const;

  amrex::Real dist(amrex::RealVect pt,
                   amrex::Real x0, amrex::Real x1,
                   amrex::Real Dx0, amrex::Real Dx1,
                   amrex::Real y0, amrex::Real y1,
                   amrex::Real Dy0, amrex::Real Dy1,
                   amrex::Real& t,
                   amrex::RealVect& spt) const;

 private:
  std::vector<amrex::Real> control_points_x;
  std::vector<amrex::Real> control_points_y;

  amrex::RealVect bc_pt_start;
  amrex::RealVect bc_pt_end;

  std::vector<amrex::Real> Dx;
  std::vector<amrex::Real> Dy;
};

}

#endif
