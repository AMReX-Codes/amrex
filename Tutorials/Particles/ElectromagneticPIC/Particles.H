#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

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

    // -------------------------------------------------
    // Data 

    using RealVector = StructOfArrays<PIdx::nattribs, 0>::RealVector;

    StructOfArrays<PIdx::nattribs, 0> attribs;
    StructOfArrays<TIdx::ntemp, 0> temp;
    
    amrex::Real charge;
    amrex::Real mass;

    // -------------------------------------------------
    // Functions

    Particles() {};

    Particles(const Particles& particles) = delete;

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& x() {return attribs.GetRealData(PIdx::x);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& y() {return attribs.GetRealData(PIdx::y);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& z() {return attribs.GetRealData(PIdx::z);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& ux() {return attribs.GetRealData(PIdx::ux);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& uy() {return attribs.GetRealData(PIdx::uy);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& uz() {return attribs.GetRealData(PIdx::uz);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& w() {return attribs.GetRealData(PIdx::w);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& ex() {return temp.GetRealData(TIdx::Ex);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& ey() {return temp.GetRealData(TIdx::Ey);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& ez() {return temp.GetRealData(TIdx::Ez);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& bx() {return temp.GetRealData(TIdx::Bx);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& by() {return temp.GetRealData(TIdx::By);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& bz() {return temp.GetRealData(TIdx::Bz);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    RealVector& ginv() {return temp.GetRealData(TIdx::ginv);}    

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& x() const {return attribs.GetRealData(PIdx::x);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& y() const {return attribs.GetRealData(PIdx::y);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& z() const {return attribs.GetRealData(PIdx::z);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& ux() const {return attribs.GetRealData(PIdx::ux);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& uy() const {return attribs.GetRealData(PIdx::uy);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& uz() const {return attribs.GetRealData(PIdx::uz);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& w() const {return attribs.GetRealData(PIdx::w);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& ex() const {return temp.GetRealData(TIdx::Ex);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& ey() const {return temp.GetRealData(TIdx::Ey);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& ez() const {return temp.GetRealData(TIdx::Ez);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& bx() const {return temp.GetRealData(TIdx::Bx);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& by() const {return temp.GetRealData(TIdx::By);}
    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& bz() const {return temp.GetRealData(TIdx::Bz);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    const RealVector& ginv() const {return temp.GetRealData(TIdx::ginv);}

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    int size() {return attribs.size();}
};

#endif
