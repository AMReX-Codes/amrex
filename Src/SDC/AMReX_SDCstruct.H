#ifndef SDCSTRUCT_H_
#define SDCSTRUCT_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

using namespace amrex;

extern "C"  {

  /**
  * \brief External routine for computing quadrature tables
  *
  * \param qtype_in
  * \param nnodes
  * \param nnodes0
  * \param nodes
  * \param nflags
  * \param qmats
  */
  void SDC_quadrature(int* qtype_in, int* nnodes, int* nnodes0,
		      amrex_real* nodes, int* nflags,amrex_real* qmats);
}
class SDCstruct
 {

 public:
   Real* qnodes;      //!< The quadrature nodes
   int*  Nflags;      //!< Flags for quadrature rules

   int Nnodes;         //!< Number of quadrature nodes
   int qtype=1;        //!< Type of quadrature nodes
   int Nsweeps=8;      //!< Number of sweeps per time step
   int Npieces;        //!< Number of terms in RHS


   // The quadrature matrices
   Real* Qall;
   Vector<Vector<Real>> Qgauss;     //!< Gauss collocation
   Vector<Vector<Real>> Qexp;       //!< Explicit
   Vector<Vector<Real>> Qimp;       //!< Implicit
   Vector<Vector<Real>> QLU;        //!< DIRK with LU trick


   // SDC storage
   Vector<MultiFab> sol;           //!< Solution at the nodes
   Vector<Vector<MultiFab> > f;    //!< Functions a nodes access by [npieces][node]
   Vector<MultiFab> res;           //!< Temp storage
   Vector<MultiFab> Ithree;        //!< Integration piece for MISDC


   /**
   * \brief Constructor
   *
   * \param Nnodes_in
   * \param Npieces_in
   * \param sol_in
   */
   SDCstruct(int Nnodes_in,int Npieces_in, MultiFab& sol_in);

   //  Sweeper routines
   void SDC_rhs_integrals(Real dt);

   void SDC_rhs_k_plus_one(MultiFab& rhs, Real dt,int m);

   void SDC_rhs_misdc(MultiFab& rhs, Real dt,int m);

 };

#endif
