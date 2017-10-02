
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

#include <cmath>

#include "AMReX_PolyGeom.H"
#include "AMReX_IntVect.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_EBISBox.H"
#include "AMReX_RealVect.H"

namespace amrex
{

  Real PolyGeom::s_tolerance = 1.0e-12;
  Real PolyGeom::s_areatolerance = 1.0e-8;
  Real PolyGeom::s_lengthtolerance = 1.0e-4;

  RealVect PolyGeom::s_vectDx = RealVect::Unit;
/*********/
  Real
  PolyGeom::
  distancePointToPlane(const RealVect& a_point,
                       const RealVect& a_normal,
                       const RealVect& a_pointOnLine)
  {
    RealVect unitNormal = a_normal;
    Real sum;
    PolyGeom::unifyVector(unitNormal, sum);
    RealVect diff = a_point - a_pointOnLine;

    Real retval = PolyGeom::dot(diff, unitNormal);

    return std::abs(retval);
  }



/*********/
  void
  PolyGeom::
  invertMatrix(Real a_AInverse[SpaceDim][SpaceDim], const Real a_A[SpaceDim][SpaceDim], bool a_test)
  {
    // get the determinant of a
    Real detA = 1.0/PolyGeom::determinantSD(a_A);
    if (std::abs(detA) < 1.0e-10)
    {
      amrex::Error("invertMatrix: determinant appears to be zero");
    }

    for (int j=0;j< SpaceDim;j++)
    {
      for (int i=0;i< SpaceDim;i++)
      {
        // get the co-factor (matrix) of A(j,i)
        //this is the same as taking the transpose THEN getting the minor matrix
        Real minor[SpaceDim-1][SpaceDim-1];
        PolyGeom::getMinor(minor, a_A, j, i);
        Real answer = determinantSDM1(minor)/detA;
        if ( (i+j)%2 == 1)
        {
          answer *= -1;
        }
        a_AInverse[i][j] = answer;
      }
    }
    if (a_test)
    {
//      Real AtimesAInverse[SpaceDim][SpaceDim];
//      bool pass = true;
      //multiply a*ainverse and see if correct
      for (int irow = 0; irow < SpaceDim; irow++)
      {
        for (int icol = 0; icol < SpaceDim; icol++)
        {
          //sum a(row, i)*ainverse(i, col).
          //should be 1 if row==col.  0 otherwise.
          int sum = 0;
          for (int ivec = 0; ivec < SpaceDim; ivec++)
          {
            sum += a_A[irow][ivec]*a_AInverse[ivec][icol];
          }
//          AtimesAInverse[irow][icol] = sum;
//          Real rightAns= 0;
//          if (icol == irow)
//          {
//            rightAns= 1;
//          }
//          if (std::abs(sum - rightAns) > 1.0e-8)
//          {
////            pass = false;
//          }
        }
      }
    }

  }

// calculate the cofactor of element (row,col)
//this is the matrix which excludes row and column row, col
  void
  PolyGeom::
  getMinor(Real      a_Aminor[SpaceDim-1][SpaceDim-1],
           const Real     a_A[SpaceDim  ][SpaceDim  ],
           int a_row, int a_col)
  {
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;
    int order = SpaceDim;
    for (int i = 0; i < order; i++ )
    {
      if ( i != a_row )
      {
        colCount = 0;
        for (int j = 0; j < order; j++ )
        {
          // when j is not the element
          if ( j != a_col )
          {
            a_Aminor[rowCount][colCount] = a_A[i][j];
            colCount++;
          }
        }
        rowCount++;
      }
    }
  }


/*********/
  void
  PolyGeom::pointToLine(RealVect& a_closestPt,
                        RealVect& a_normal,
                        const RealVect& a_point,
                        const RealVect& a_pointOnLine,
                        const RealVect& a_direction)
  {
    //get the length of the input direction
    Real  length2 = dot(a_direction, a_direction);
    assert(length2 > 0.0);
    // Get a_point relative to a point on the centerline
    // delta = p - p0
    RealVect delta;

    delta = a_point;
    delta -= a_pointOnLine;

    // Get the dot product of the relative vector (above) and the centerline
    // direction
    //dot = (v dot (p-p0))/v^2 = t0
    Real dotprod = dot(a_direction, delta);
    dotprod /= length2;

    // Find the vector from a_point to the point on the centerline closest to
    // a_point
    //n = p0 + v t0
    a_closestPt = a_direction;
    a_closestPt *= dotprod;
    a_closestPt += a_pointOnLine;

    //normal vector = (n - p)/|n-p|
    a_normal = a_closestPt;
    a_normal -= a_point;
    Real sum;
    unifyVector(a_normal, sum);
  }
/*******************************/
  RealVect
  PolyGeom::normal(const VolIndex& a_vof,
                   const EBISBox& a_ebisBox,
                   const Real& a_bndryArea)
  {
    RealVect retval;
    const IntVect& iv = a_vof.gridIndex();
    if (a_ebisBox.isRegular(iv))
    {
      retval = BASISREALV(0);
    }
    else if (a_ebisBox.isCovered(iv))
    {
      retval =  BASISREALV(0);
    }
    else
    {
      assert(a_ebisBox.isIrregular(iv));
      Real irregArea = a_bndryArea;
      if (irregArea> 0.0)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real hiArea = a_ebisBox.sumArea(a_vof, idir, Side::Hi);
          Real loArea = a_ebisBox.sumArea(a_vof, idir, Side::Lo);
          retval[idir] = ((hiArea - loArea)/irregArea);
        } //end loop over directions
      }
      else
      {
        retval =  RealVect::Zero;
      }

    } //end else (vof is irregular)
    return retval;
  }
/*******************************/
  Real
  PolyGeom::bndryArea(const VolIndex& a_vof,
                      const EBISBox& a_ebisBox)
  {
    Real retval=0.;

    Real irregArea=0.0;
    for (int idir=0; idir < SpaceDim; idir++)
    {
      Real hiArea = a_ebisBox.sumArea(a_vof, idir, Side::Hi);
      Real loArea = a_ebisBox.sumArea(a_vof, idir, Side::Lo);
      irregArea += (hiArea-loArea)*(hiArea-loArea);
    }
    retval = sqrt(irregArea);
    return retval;;
  }
/**********************************************/
  Tuple<int, SpaceDim-1>
    PolyGeom::computeTanDirs(int upDir)
  {
    Tuple<int, SpaceDim-1> tanDirs;
    int itan = 0;
    for (int idir =0; idir < SpaceDim; idir++)
    {
      if (idir != upDir)
      {
        tanDirs[itan] = idir;
        itan++;
      }
    }
    return tanDirs;
  }
/**************/
/**************/
//solve for component of x in A x = rhs
//using Cramer's rule.
/**************/
  Real
  PolyGeom::matrixSolveComp(const Array<Array<Real> >& a_A,
                            const Array<Real>& a_rhs,
                            const int& a_icomp)
  {
    int nVar = a_A.size();
    assert(a_rhs.size() == nVar);

    //set topmat ij to be == Aij where j != a_icomp,
    //set topmat ij to be == rhsj where j == a_icomp,
    assert(a_icomp >= 0);
    assert(a_icomp< nVar);

    Array<Array<Real> > topMat(a_A);
    for (int irow = 0; irow< nVar; irow++)
    {
      topMat[irow][a_icomp] = a_rhs[irow];
    }
    //answer is det(topmat)/det(A)
    Real numer = determinant(topMat);
    Real denom = determinant(a_A);
    Real eps = 1.0e-10;
    if (std::abs(denom) < eps)
      amrex::Error("MatrixSolveComp has encountered unsolvable matrix");
    Real result = numer/denom;
    return result;
  }
/***************/
/***************/
  Real
  PolyGeom::determinant(const Array<Array< Real> >& a_A)
  {
    int nVar = a_A.size();
    assert(nVar >= 2);
    assert(nVar <= 4);
    Real det = 0.0;
    if (nVar == 2)
    {
      det = a_A[0][0]*a_A[1][1] -  a_A[0][1]*a_A[1][0];
    }
    else if (nVar == 3)
    {
      det
        =a_A[0][0]*(a_A[1][1]*a_A[2][2]-a_A[2][1]*a_A[1][2])
        -a_A[0][1]*(a_A[1][0]*a_A[2][2]-a_A[2][0]*a_A[1][2])
        +a_A[0][2]*(a_A[1][0]*a_A[2][1]-a_A[2][0]*a_A[1][1]);
    }
    else
    {
      assert(nVar == 4);
      det=a_A[0][0]*(a_A[1][1]*(a_A[2][2]*a_A[3][3]-a_A[3][2]*a_A[2][3])-
                     a_A[1][2]*(a_A[2][1]*a_A[3][3]-a_A[3][1]*a_A[2][3])+
                     a_A[1][3]*(a_A[2][1]*a_A[3][2]-a_A[3][1]*a_A[2][2]))

        -a_A[0][1]*(a_A[1][0]*(a_A[2][2]*a_A[3][3]-a_A[3][2]*a_A[2][3])-
                    a_A[1][2]*(a_A[2][0]*a_A[3][3]-a_A[3][0]*a_A[2][3])+
                    a_A[1][3]*(a_A[2][0]*a_A[3][2]-a_A[3][0]*a_A[2][2]))

        +a_A[0][2]*(a_A[1][0]*(a_A[2][1]*a_A[3][3]-a_A[3][1]*a_A[2][3])-
                    a_A[1][1]*(a_A[2][0]*a_A[3][3]-a_A[3][0]*a_A[2][3])+
                    a_A[1][3]*(a_A[2][0]*a_A[3][1]-a_A[3][0]*a_A[2][1]))

        -a_A[0][3]*(a_A[1][0]*(a_A[2][1]*a_A[3][2]-a_A[3][1]*a_A[2][2])-
                    a_A[1][1]*(a_A[2][0]*a_A[3][2]-a_A[3][0]*a_A[2][2])+
                    a_A[1][2]*(a_A[2][0]*a_A[3][1]-a_A[3][0]*a_A[2][1]));
    }
#if BL_SPACEDIM==2
#elif BL_SPACEDIM==3
#else
    bogus_ch_spacedim_macro();
#endif
    return det;
  }

/***************/
  Real
  PolyGeom::
  determinantSD(const Real a_A[SpaceDim][SpaceDim])
  {
    Real det = 0.0;

    AMREX_D_TERM(
      det = a_A[0][0];,

      det = a_A[0][0]*a_A[1][1]
      - a_A[0][1]*a_A[1][0];,

      det = a_A[0][0] * (a_A[1][1]*a_A[2][2] - a_A[2][1]*a_A[1][2])
      - a_A[0][1] * (a_A[1][0]*a_A[2][2] - a_A[2][0]*a_A[1][2])
      + a_A[0][2] * (a_A[1][0]*a_A[2][1] - a_A[2][0]*a_A[1][1]));

    return det;
  }
/***************/
  Real
  PolyGeom::
  determinantSDM1(const Real a_A[SpaceDim-1][SpaceDim-1])
  {
    Real det = 0.0;

    AMREX_D_TERM(
      det = 0.0;,

      det = a_A[0][0];,

      det = a_A[0][0]*a_A[1][1]
      - a_A[0][1]*a_A[1][0]);

    return det;
  }
/*******************************/
//solves the linear system
// for each point j
// ni xi - alpha = 0
//with the normalizing eqn
// \sum ni = 1
//for ni and alpha
//always return nromal[updir] >= 0
/*******************************/
  void
  PolyGeom::computeNormalAndAlpha(Real& a_alpha,
                                  RealVect& a_normal,
                                  const int& a_upDir,
                                  const Tuple<RealVect, BL_SPACEDIM>& a_poly)
  {
#if BL_SPACEDIM == 2
    RealVect pt0 = a_poly[0];
    RealVect pt1 = a_poly[1];
    //make pt0 the point with smaller x
    if (pt0[0] > pt1[0])
    {
      RealVect temp = pt0;
      pt0 = pt1;
      pt1 = temp;
    }
    Real signx = 1.0;
    if (pt0[1] < pt1[1])
    {
      signx = -1.0;
    }
    if (std::abs(pt0[0] -pt1[0]) < s_tolerance)
    {
      a_normal[1] = 0.0;
      a_normal[0] = signx;
    }
    else if (std::abs(pt0[1] -pt1[1]) < s_tolerance)
    {
      if (pt0[0] < pt1[0])
      {
        signx = -1.0;
      }
      a_normal[1] = signx;
      a_normal[0] = 0.0;
    }
    else
    {
      //case where points differ in both coords.
      Real ratio = -(pt0[0] -pt1[0])/(pt0[1]-pt1[1]);
      Real denom = 1.0+ ratio*ratio;
      a_normal[0] = signx/sqrt(denom);
      a_normal[1] = ratio*a_normal[0];
    }
#elif BL_SPACEDIM == 3
    //compute cross product of two vectors formed by difference
    //between pt0 and pt2 and pt1 and pt2
    RealVect xvec1 = a_poly[0];
    RealVect xvec0 = a_poly[1];
    xvec1 -= a_poly[2];
    xvec0 -= a_poly[2];
    //a_normal = xvec1 x xvec0
    a_normal = cross(xvec1, xvec0);
    //now normalize so that normal is a normal vector
    Real sum = 0;
    unifyVector(a_normal, sum);
#else
    bogus_ch_spacedim();
#endif

    //if normal[upDir] < 0, reverse the signs
    if (a_normal[a_upDir] < 0.0)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        a_normal[idir] *= -1.0;
      }
    }

    //by definition alpha = sum ni xi for each x in the poly
    a_alpha = 0.0;
    RealVect zeroPt = a_poly[0];
    for (int idir = 0; idir< SpaceDim; idir++)
    {
      a_alpha += a_normal[idir]*zeroPt[idir];
    }
  }
//compute the dot product between xvec0 and xvec1
// (returns xvec1 dot xvec0)
  Real
  PolyGeom::dot(const RealVect& a_xvec1, const RealVect& a_xvec0)
  {
    Real retval= 0.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += a_xvec1[idir]*a_xvec0[idir];
    }
    return retval;
  }
//compute the cross product between xvec0 and xvec1
// (returns xvec1 x xvec0)
  RealVect
  PolyGeom::cross(const RealVect& a_xvec1,const RealVect& a_xvec0)
  {
    RealVect retval;
#if BL_SPACEDIM==3
    retval[0] =   a_xvec0[1]*a_xvec1[2] - a_xvec0[2]*a_xvec1[1];
    retval[1] = -(a_xvec0[0]*a_xvec1[2] - a_xvec0[2]*a_xvec1[0]);
    retval[2] =   a_xvec0[0]*a_xvec1[1] - a_xvec0[1]*a_xvec1[0];
#else
    Real scalar = a_xvec0[0]*a_xvec1[1] - a_xvec0[1]*a_xvec1[0];
    retval[0] =   scalar;
    retval[1] =   scalar;
#endif
    return retval;
  }
/*******************************/
///sort vector whose comps are all positive.
/*******************************/
  void
  PolyGeom::sortVector(RealVect& a_vect, IntVect& a_ivmap)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      a_ivmap[idir] = idir;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      assert(a_vect[idir] >= 0.0);
      for (int jdir = 1; jdir < SpaceDim; jdir++)
      {
        //if bigger normal is first, switch them
        if (a_vect[jdir-1]  > a_vect[jdir])
        {
          {
            Real temp = a_vect[jdir-1];
            a_vect[jdir-1] = a_vect[jdir];
            a_vect[jdir]   = temp;
          }
          {
            int itemp = a_ivmap[jdir-1];
            a_ivmap[jdir-1] = a_ivmap[jdir];
            a_ivmap[jdir] = itemp;
          }
        }
      }
    }
  }

/*******************************/
///make vector all pos and return the signs
/*******************************/
  void
  PolyGeom::posifyVector(RealVect& a_vect, IntVect& a_sign)
  {
    a_sign = IntVect::TheUnitVector();
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_vect[idir] < 0.0)
      {
        a_sign[idir] = -1;
        a_vect[idir] = std::abs(a_vect[idir]);
      }
    }
  }

/*******************************/
///make vector into a unit vector
/*******************************/
  void
  PolyGeom::unifyVector(RealVect& a_normal, Real& a_sumSquare)
  {
    Real sum = 0.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      sum += a_normal[idir]*a_normal[idir];
    }
    sum = sqrt(sum);
    //divide both by the sum of the square
    a_sumSquare = sum;
    if (sum > s_tolerance)
    {
      a_normal /= sum;
    }
  }
/*******************************/
/*******************************/
  Real
  PolyGeom::computeVolume(const Real& a_alpha,
                          const RealVect& a_normal)
  {
    Real alphanew = a_alpha;
    RealVect normnew = a_normal;
    Real sumSquare;
    unifyVector(normnew, sumSquare);
    //Transform alpha for normals that did are not unit vectors
    alphanew /= sumSquare;
    //make the norms all positive
    IntVect signnorm;
    posifyVector(normnew, signnorm);
    //need to alter alpha due to n = -n
    //coord transformation
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (signnorm[idir] == -1)
      {
        //add the normal component that used
        //to be negative
        //use normnew because stuff has been
        //rescaled due to the unify thing
        alphanew += normnew[idir];
      }
    }
    IntVect ivmap;
    sortVector(normnew, ivmap);

    //sandzvolume computes the volume under
    //the curve (covered volume)
    Real volFrac = 1.0 - sAndZVolume(alphanew, normnew);

    return volFrac;
  }

/*******************************/
/*******************************/
/*******************************/
//use bisection algorithm to iterate to alph
/*******************************/
  Real
  PolyGeom::computeAlpha(const Real& a_volFrac,
                         const RealVect& a_normal)
  {
    Real retval = -1.0;

    assert(a_volFrac >= 0.0);
    assert(a_volFrac <= 1.0);

    Real alphaMax = 0.0;
    Real alphaMin = 0.0;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_normal[idir] < 0.0)
      {
        alphaMin += a_normal[idir];
      }
      else
      {
        alphaMax += a_normal[idir];
      }
    }

    Real alphaLo = alphaMin;
    Real alphaHi = alphaMax;
    Real alphaMid= 0.5*(alphaLo+alphaHi);
    Real funcLo = computeVolume(alphaLo, a_normal)  - a_volFrac;
    Real funcHi = computeVolume(alphaHi, a_normal)  - a_volFrac;

    if (abs(funcHi) < s_tolerance)
    {
      retval = alphaHi;
    }
    else if (abs(funcLo) < s_tolerance)
    {
      retval = alphaLo;
    }
    else
    {
      if (funcLo*funcHi > 0.0)
      {
        amrex::Error("computeAlpha (bisection iteration): root not bracketed");
      }
      Real err = std::abs(alphaHi-alphaLo);
      Real funcMid = computeVolume(alphaMid, a_normal)  - a_volFrac;
      int itermax = 1000;
      int iteration = 0;
      while ((iteration < itermax) && (err > s_tolerance))
      {
        Real alphaOld = alphaMid;
        funcMid  = computeVolume(alphaMid, a_normal)  - a_volFrac;
        if (funcLo*funcMid > 0.0)
        {
          alphaLo = alphaMid;
          funcLo = funcMid;
        }
        else
        {
          assert(funcHi*funcMid >= 0.0);
          alphaHi = alphaMid;
          funcHi = funcMid;
        }
        iteration++;
        alphaMid = 0.5*(alphaLo+alphaHi);
        err = std::abs(alphaMid - alphaOld);
      }
      if ((alphaMid < alphaMin) || alphaMid > alphaMax)
        amrex::Error("bisectiter: solution out of bounds");
      retval = alphaMid;

      /*
        if bisection did not converge,
        Houston, we have a problem
      */
      if (iteration >= itermax)
        amrex::Error("bisectiter: did not converge!");
    }

    return retval;
  }
/*******************************/
  Real
  PolyGeom::threeDFunc(const Real& a_arg)
  {
    Real retval = 0.0;
    if (a_arg > 0.0)
      retval = a_arg*a_arg*a_arg;
    return retval;
  }
/*******************************/
  Real
  PolyGeom::twoDFunc(const Real& a_arg)
  {
    Real retval = 0.0;
    if (a_arg > 0.0)
      retval = a_arg*a_arg;
    return retval;
  }
/*******************************/
  Real
  PolyGeom::computeAnyVolume(const Real& a_alpha,
                             const Real& a_norm0,
                             const Real& a_norm1,
                             const Real& a_norm2)
  {
    // These result from expanding the S and Z equation 2 and taking pains
    // to write them in a numerically stable form for each distinct range
    // of alpha.  There are no problems with a_norm0 and a_norm1 being small
    // or zero and 2D is handled simply by setting a_norm0 to zero and calling
    // this method.

    Real vf = 0.0;
    Real alpha = a_alpha;
    int flip = 0;

    if (alpha >= 0.5*(a_norm0+a_norm1+a_norm2))
    {
      alpha = a_norm0+a_norm1+a_norm2 - alpha;
      flip = 1;
    }

    if (alpha < 0)
    {
      vf = 0.0;
    }
    else
      if (alpha >= 0 && alpha < a_norm0)
      {
        vf = 1.0/6.0 * (alpha*alpha/a_norm1/a_norm2) * (alpha/a_norm0);
      }
      else
        if (alpha >= a_norm0 && alpha < a_norm1)
        {
          vf = 1.0/6.0 * (a_norm0*a_norm0)     /a_norm1/a_norm2
            + 1.0/2.0 * (alpha*(alpha-a_norm0)/a_norm1/a_norm2);
        }
        else
          if (a_norm0+a_norm1 <= a_norm2)
          {
            if (alpha >= a_norm1 && alpha < a_norm0+a_norm1)
            {
              vf = 1.0/6.0 * (a_norm0*a_norm0)/a_norm1/a_norm2
                + 1.0/2.0 * ((a_norm1-a_norm0)/a_norm2)
                + 1.0/2.0 * (alpha-a_norm1)/a_norm1
                * (alpha+a_norm1-a_norm0)/a_norm2
                - 1.0/6.0 * (alpha-a_norm1)/a_norm0
                * (alpha-a_norm1)/a_norm1
                * (alpha-a_norm1)/a_norm2;
            }
            else
              if (alpha >= a_norm0+a_norm1)
              {
                vf = (alpha - 1.0/2.0*(a_norm0+a_norm1))/a_norm2;
              }
          }
          else
          {
            if (alpha >= a_norm1 && alpha < a_norm2)
            {
              vf = 1.0/6.0 * (a_norm0*a_norm0)/a_norm1/a_norm2
                + 1.0/2.0 * ((a_norm1-a_norm0)/a_norm2)
                + 1.0/2.0 * (alpha-a_norm1)/a_norm1
                * (alpha+a_norm1-a_norm0)/a_norm2
                - 1.0/6.0 * (alpha-a_norm1)/a_norm0
                * (alpha-a_norm1)/a_norm1
                * (alpha-a_norm1)/a_norm2;
            }
            else
              if (alpha >= a_norm2)
              {
                vf = 1.0/6.0 * (a_norm0*a_norm0)/a_norm1/a_norm2
                  + 1.0/2.0 * ((a_norm1-a_norm0)/a_norm2)
                  + 1.0/2.0 * (alpha-a_norm1)/a_norm1
                  * (alpha+a_norm1-a_norm0)/a_norm2
                  - 1.0/6.0 * (alpha-a_norm1)/a_norm0
                  * (alpha-a_norm1)/a_norm1
                  * (alpha-a_norm1)/a_norm2
                  - 1.0/6.0 * (alpha-a_norm2)/a_norm0
                  * (alpha-a_norm2)/a_norm1
                  * (alpha-a_norm2)/a_norm2;
              }
          }

    if (flip == 1)
    {
      vf = 1.0 - vf;
    }

    return(vf);
  }
/*******************************/
/*******************************/
  Real
  PolyGeom::sAndZVolume(const Real& a_alpha,
                        const RealVect& a_normal)
  {
    Real retval = -1;
    //normals have to be positive.
    for (int idir = 0; idir < SpaceDim; idir++)
      assert(a_normal[idir] >= 0.0);
    //normals have to be ordered by size
    for (int idir = 1; idir < SpaceDim; idir++)
      assert(a_normal[idir] >= a_normal[idir-1]);
    Real norm0,norm1,norm2;

    if (SpaceDim == 2)
    {
      norm0 = 0.0;
      norm1 = a_normal[0];
      norm2 = a_normal[1];
    }
    else
    {
      norm0 = a_normal[0];
      norm1 = a_normal[1];
      norm2 = a_normal[2];
    }

    retval = computeAnyVolume(a_alpha,norm0,norm1,norm2);

    return retval;
  }
/*******************************/
/*******************************/
  void
  PolyGeom::setTolerance(const Real& a_tolerance)
  {
    s_tolerance = a_tolerance;
  }
/*******************************/
/*******************************/
  void
  PolyGeom::setVectDx(const RealVect& a_vectDx)
  {
    s_vectDx = a_vectDx;
  }
/*******************************/
/*******************************/
  void
  PolyGeom::setVolumeTolerance(const Real& a_tolerance)
  {
    s_tolerance = a_tolerance;
  }
/*******************************/
/*******************************/
  void
  PolyGeom::setAreaTolerance(const Real& a_tolerance)
  {
    s_areatolerance = a_tolerance;
  }
/*******************************/
/*******************************/
  void
  PolyGeom::setLengthTolerance(const Real& a_tolerance)
  {
    s_lengthtolerance = a_tolerance;
  }
/*******************************/
/*******************************/
  const Real&
  PolyGeom::getTolerance()
  {
    return s_tolerance;
  }
  const RealVect&
  PolyGeom::getVectDx()
  {
    return s_vectDx;
  }
  const Real&
  PolyGeom::getVolumeTolerance()
  {
    return s_tolerance;
  }
  const Real&
  PolyGeom::getAreaTolerance()
  {
    return s_areatolerance;
  }
  const Real&
  PolyGeom::getLengthTolerance()
  {
    return s_lengthtolerance;
  }
/*******************************/
/*******************************/
//volume of a tetrahedron.  all normals must be > s_tolerance
//(includes positivity)
/*******************************/
  Real
  PolyGeom::tetVolume(const RealVect& a_normal,
                      const Real& a_alpha)
  {
    Real retval = 1.0/6.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      assert(a_normal[idir] > s_tolerance);
      retval *= a_alpha/a_normal[idir];
    }
    return retval;
  }
/*******************************/
//return the centroid of a tetrahedron with
//the given normal and alpha. if any of the normals
//are zero,  not really a tet so it returns
// \vec -1
/*******************************/
  RealVect
  PolyGeom::tetCentroid(const RealVect& a_normal,
                        const Real& a_alpha)
  {
    RealVect retval = RealVect::Zero;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real denom = 4.0*a_normal[idir];
      Real numer = a_alpha;
      if (denom >  s_tolerance)
      {
        retval[idir] = numer/denom;
      }
      else
      {
        //not really a tet.  return nonsense
        retval[idir] = -1.0;
      }
    }
    return retval;
  }
}

