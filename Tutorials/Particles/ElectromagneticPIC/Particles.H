#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include <CudaManagedAllocator.H>
#include <StructOfArrays.H>

/**

   Note that we use two StructOfArrays to store the particle data,
   due to the current limitation that Thrust's zip_iterators only
   work with 10 or less elements.

 */
struct PIdx
{
    enum { // These are "permanent" attributes that move around with the particles
        x = 0,
        y, z, ux, uy, uz, w,
        nattribs
    };
};

struct TIdx
{
    enum { // These are "temporaries" that are recomputed each time step
        Ex = 0,
        Ey, Ez, Bx, By, Bz, ginv,
        ntemp
    };
};

struct Particles {

    using RealVector = StructOfArrays<PIdx::nattribs, 0>::RealVector;

    StructOfArrays<PIdx::nattribs, 0> attribs;
    StructOfArrays<TIdx::ntemp, 0> temp;
    
    amrex::Real charge;
    amrex::Real mass;

    Particles() {};

    Particles(const Particles& particles) = delete;

    RealVector& x() {return attribs.GetRealData(PIdx::x);}
    RealVector& y() {return attribs.GetRealData(PIdx::y);}
    RealVector& z() {return attribs.GetRealData(PIdx::z);}

    RealVector& ux() {return attribs.GetRealData(PIdx::ux);}
    RealVector& uy() {return attribs.GetRealData(PIdx::uy);}
    RealVector& uz() {return attribs.GetRealData(PIdx::uz);}

    RealVector& w() {return attribs.GetRealData(PIdx::w);}

    RealVector& ex() {return temp.GetRealData(TIdx::Ex);}
    RealVector& ey() {return temp.GetRealData(TIdx::Ey);}
    RealVector& ez() {return temp.GetRealData(TIdx::Ez);}

    RealVector& bx() {return temp.GetRealData(TIdx::Bx);}
    RealVector& by() {return temp.GetRealData(TIdx::By);}
    RealVector& bz() {return temp.GetRealData(TIdx::Bz);}

    RealVector& ginv() {return temp.GetRealData(TIdx::ginv);}    

    const RealVector& x() const {return attribs.GetRealData(PIdx::x);}
    const RealVector& y() const {return attribs.GetRealData(PIdx::y);}
    const RealVector& z() const {return attribs.GetRealData(PIdx::z);}

    const RealVector& ux() const {return attribs.GetRealData(PIdx::ux);}
    const RealVector& uy() const {return attribs.GetRealData(PIdx::uy);}
    const RealVector& uz() const {return attribs.GetRealData(PIdx::uz);}

    const RealVector& w() const {return attribs.GetRealData(PIdx::w);}

    const RealVector& ex() const {return temp.GetRealData(TIdx::Ex);}
    const RealVector& ey() const {return temp.GetRealData(TIdx::Ey);}
    const RealVector& ez() const {return temp.GetRealData(TIdx::Ez);}

    const RealVector& bx() const {return temp.GetRealData(TIdx::Bx);}
    const RealVector& by() const {return temp.GetRealData(TIdx::By);}
    const RealVector& bz() const {return temp.GetRealData(TIdx::Bz);}

    const RealVector& ginv() const {return temp.GetRealData(TIdx::ginv);}

    int size() {return attribs.size();}
};

#endif
