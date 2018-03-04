
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _POLYGEOM_H_
#define _POLYGEOM_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_Tuple.H"
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include "AMReX_VolIndex.H"

namespace amrex
{
  class EBISBox;
///
/**
   PolyGeom encapsulates the functions
   to generate centroids, normals, etc given
   the geometric information.
   It is meant to be used simply as a static class.
   Its only data member (a static) is its tolerance.
*/
  class PolyGeom
  {
  public:

    ///
    /**
       return the inverse of the input matrix.   done using cramers rule
    */
    static void invertMatrix(Real      a_AInverse[SpaceDim][SpaceDim],
                             const Real       a_A[SpaceDim][SpaceDim], bool a_test);

    ///
    /**
     // calculate the cofactor of element (row,col)
     //this is the matrix which excludes row and column row, col
     */
    static void getMinor(Real      a_Aminor[SpaceDim-1][SpaceDim-1],
                         const Real     a_A[SpaceDim  ][SpaceDim  ],
                         int a_row, int acol);

    /****/
    static void
    getJacobianAndInverse(Real a_Jacobian[SpaceDim][SpaceDim],
                          Real a_Jinverse[SpaceDim][SpaceDim],
                          RealVect& a_normal,
                          RealVect a_tangents[SpaceDim-1])
      {
        //make the first row = the normal
        //make next row the tangent vectors
        for (int icol = 0; icol < SpaceDim; icol++)
        {
          a_Jacobian[0][icol] = a_normal[icol];
          for (int irow  = 1; irow < SpaceDim; irow++)
          {
            a_Jacobian[irow][icol] = a_tangents[irow-1][icol];
          }
        }
        //compute the inverse of the matrix.   the true is turn on
        //the test function.
        PolyGeom::invertMatrix(a_Jinverse, a_Jacobian, true);
      }

//given a normal vector, get spacedim-1 tangential vectors that are
//unit vectors and
    static void
    getTangentVectors( RealVect              a_tangen[AMREX_SPACEDIM-1],
                       const RealVect&       a_normal)
      {
#if AMREX_SPACEDIM==2
        a_tangen[0][1] =  a_normal[0];
        a_tangen[0][0] = -a_normal[1];
#else
        //get some vector that is not parallel to a_normal
        RealVect notNormal = a_normal + RealVect::Unit;
        Real dotProd = PolyGeom::dot(notNormal, a_normal);
        //there is one bizzare case where I have just scaled up a_normal
        //in this case, add something diff to each component.
        if (std::abs(dotProd) < 1.0e-6)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            notNormal[idir] += Real(idir+1);
          }
          dotProd = PolyGeom::dot(notNormal, a_normal);
          BL_ASSERT(std::abs(dotProd) > 1.0e-6);
        }
        Real sumSquare;
        PolyGeom::unifyVector(notNormal, sumSquare);
        //make the first tangential vector = cross(notnormal, normal);
        a_tangen[0] = PolyGeom::cross(notNormal,   a_normal);
        //make the second one the cross between that and normal
        a_tangen[1] = PolyGeom::cross(a_tangen[0], a_normal);
#endif

        //make them all unit vectors
        for (int itan = 0; itan < AMREX_SPACEDIM-1; itan++)
        {
          Real tmpsumSquare;
          PolyGeom::unifyVector(a_tangen[itan], tmpsumSquare);
        }
      }

    ///
    static Real determinantSD(const Real a_A[SpaceDim][SpaceDim]);

    ///
    static Real determinantSDM1(const Real a_A[SpaceDim-1][SpaceDim-1]);

    ///
    /**
       find equation of line
       between a line formed by a_pointOnLine = p0 and a_direction=v
       l(t) = p0 + v t
       and a point in space a_point.
       returns closest point on line and the normal between the
       closest point and the point in space
    */
    static void
    pointToLine(RealVect& a_closestPt,
                RealVect& a_normal,
                const RealVect& a_point,
                const RealVect& a_pointOnLine,
                const RealVect& a_direction);

    ///return distance from point to a plane
    static Real
    distancePointToPlane(const RealVect& a_point,
                         const RealVect& a_normal,
                         const RealVect& a_pointOnLine);

    ///
    /**
     */
    static void setTolerance(const Real& a_tolerance);

    ///
    /**
     */
    static void setVectDx(const RealVect& a_vectDx);

    ///
    /**
     */
    static void setVolumeTolerance(const Real& a_tolerance);

    ///
    /**
     */
    static void setAreaTolerance(const Real& a_tolerance);

    ///
    /**
     */
    static void setLengthTolerance(const Real& a_tolerance);

    ///
    /**
       compute the cross product between xvec0 and xvec1
       (returns xvec1 x xvec0)
    */
    static RealVect cross(const RealVect& a_xvec1,
                          const RealVect& a_xvec0);

    ///
    /**
       compute the dot product between xvec0 and xvec1
       (returns xvec1 dot xvec0)
    */
    static Real dot(const RealVect& a_xvec1,
                    const RealVect& a_xvec0);

    ///
    /**
     */
    static const Real& getTolerance();

    ///
    /**
     */
    static const RealVect& getVectDx();

    ///
    /**
     */
    static const Real& getVolumeTolerance();

    ///
    /**
     */
    static const Real& getAreaTolerance();

    ///
    /**
     */
    static const Real& getLengthTolerance();

    ///
    /**
       Return the normal to the boundary at the
       input VoF.  If said normal is undefined,
       returns the zero vector.
    */
    static RealVect normal(const VolIndex& a_vof,
                           const EBISBox& a_ebisBox,
                           const Real& a_bndryArea);

    ///
    /**
     */
    static Real bndryArea(const VolIndex& a_vof,
                          const EBISBox& a_ebisBox);

    ///
    /**
       Return the inhomogeneous component of the plane
       eqution ni xi = alpha.
       Volume fraction must be positive and <= 1.0
       and the normals can be in any order.
    */
    static Real computeAlpha(const Real& a_volFrac,
                             const RealVect& a_normal);

    ///
    /**
       Compute the volume fraction in the cell given
       alpha and the normal.
    */
    static Real computeVolume(const Real& a_alpha,
                              const RealVect& a_normal);

    ///
    /**
       Get the normal and alpha given the points in a polygon.
       Solves the linear system \\
       ni xi - alpha = 0\\
       \sum ni*ni = 1\\
       for ni and alpha.
       Always returns nromal[updir] >= 0
    */
    static void computeNormalAndAlpha(Real& a_alpha,
                                      RealVect& a_normal,
                                      const int& a_upDir,
                                      const Tuple<RealVect, AMREX_SPACEDIM>& a_poly);

    // solve ax = b (gives icomp component of x using cramers rule).
    static Real
    matrixSolveComp(const Vector<Vector<Real> >& a_A,
                    const  Vector<Real>& a_rhs,
                    const int& a_icomp);

    //return the determinant of a
    static Real
    determinant(const Vector< Vector<Real> >& a_A);

    /**
       andzvolume computes the volume under
       the curve (covered volume).
       Does so using Scardovelli
       and Zaleski's direct formulae.
       The normal has to be all positive numbers.
       The normals also must be ordered from lowest
       to highest.
    */
    static Real sAndZVolume(const Real& a_alpha,
                            const RealVect& a_normal);

    ///
    PolyGeom()
      {;}
    ///
    ~PolyGeom()
      {;}

    ///
    /**
       Sort vector from low to high.  All comps must be positive.
       Return reverse mapping
    */
    static void sortVector(RealVect& vect, IntVect& ivmap);

    ///make vector all pos and return the signs
    static void posifyVector(RealVect& vect, IntVect& signs);

    ///make vector into a unit vector and return the sum of squares
    static void unifyVector(RealVect& normal, Real& sumSquare);

    ///
    static Tuple<int, AMREX_SPACEDIM-1> computeTanDirs(int upDir);

  protected:

    static Real   twoDFunc(const Real& arg);
    static Real threeDFunc(const Real& arg);

    static Real computeAnyVolume(const Real& a_alpha,
                                 const Real& a_norm0,
                                 const Real& a_norm1,
                                 const Real& a_norm2);

    static RealVect s_vectDx;

    static Real s_tolerance;

    static Real s_lengthtolerance;

    static Real s_areatolerance;

    static RealVect tetCentroid(const RealVect& normal,
                                const Real& alpha);

    static Real tetVolume(const RealVect& normal,
                          const Real& alpha);


  private:
  };

}
#endif
