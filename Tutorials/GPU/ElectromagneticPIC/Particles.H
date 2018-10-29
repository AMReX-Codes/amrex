#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include <StructOfArrays.H>

/*
    Copy of Particles for GPU use.
*/

    struct ParticlesData
    {
       // Matching access functions
       AMREX_CUDA_DEVICE
       amrex::Real* x() {return m_x;}
       AMREX_CUDA_DEVICE
       amrex::Real* y() {return m_y;}
       AMREX_CUDA_DEVICE
       amrex::Real* z() {return m_z;}

       AMREX_CUDA_DEVICE
       amrex::Real* ux() {return m_ux;}
       AMREX_CUDA_DEVICE
       amrex::Real* uy() {return m_uy;}
       AMREX_CUDA_DEVICE
       amrex::Real* uz() {return m_uz;}

       AMREX_CUDA_DEVICE
       amrex::Real* w() {return m_w;}

       AMREX_CUDA_DEVICE
       amrex::Real* ex() {return m_Ex;}
       AMREX_CUDA_DEVICE
       amrex::Real* ey() {return m_Ey;}
       AMREX_CUDA_DEVICE
       amrex::Real* ez() {return m_Ez;}
 
       AMREX_CUDA_DEVICE
       amrex::Real* bx() {return m_Bx;}
       AMREX_CUDA_DEVICE
       amrex::Real* by() {return m_By;}
       AMREX_CUDA_DEVICE
       amrex::Real* bz() {return m_Bz;}

       AMREX_CUDA_DEVICE
       amrex::Real* ginv() {return m_ginv;}    

       AMREX_CUDA_DEVICE
       int size() {return m_size;}

       public:
       // Data
       int m_size;
       amrex::Real *m_x,  *m_y,  *m_z;  
       amrex::Real *m_ux, *m_uy, *m_uz;
       amrex::Real *m_Ex, *m_Ey, *m_Ez;
       amrex::Real *m_Bx, *m_By, *m_Bz;
       amrex::Real *m_w,  *m_ginv;
    };

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

    ParticlesData data()
    {
       ParticlesData pData;

       pData.m_x     = x().data();
       pData.m_y     = y().data();
       pData.m_z     = z().data();
       pData.m_ux    = ux().data();
       pData.m_uy    = uy().data();
       pData.m_uz    = uz().data();
       pData.m_w     = w().data();
       pData.m_Ex    = ex().data();
       pData.m_Ey    = ey().data();
       pData.m_Ez    = ez().data();
       pData.m_Bx    = bx().data();
       pData.m_By    = by().data();
       pData.m_Bz    = bz().data();
       pData.m_ginv  = ginv().data();
       pData.m_size  = size();

       return pData;
    }

};

#endif
