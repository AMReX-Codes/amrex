#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LSQUARES_H_
#define _LSQUARES_H_

#include "Vector.H"
#include "REAL.H"
#include "RealVect.H"
#include "NamespaceHeader.H"
///
/**
    This class solves least squares problems
 */
class LSquares
{
public:
  void LeastSquares(Real** A,
                    Vector<Real>&x,
                    const Vector<Real>&rhs);

  /// gaussian elimination with partial pivoting
  int gaussElim(Real**A,
                Vector<Real>& rhs);

  void swapRows(Real** A,
                const int& rowi,
                const int& rowj,
                const int&numberOfCols);

  void swapRows(Vector<Real>& rhs,
                const int& currRow,
                const int& pivot);

  int findPivot(Real** A,
                const int& currCol,
                const int& currRow,
                const int& numRows,
                int& pivot);

  void addRows(Vector<Real>& rhs,
               const int& rowi,
               const Real& alpha,
               const int& rowj);

  void addRows(Real** A,
               const int& rowi,
               const Real& alpha,
               const int& rowj,
               const int& numberOfCols);

  void timesBeta(Vector<Real>&rhs,
                 const int& currRow,
                 const Real& Beta);

  void timesBeta(Real** A,
                 const int& rowi,
                 const Real& Beta,
                 const int& numberOfcols);

  void transpose(Real** a_A,
                 Real ** a_Atrans,
                 const int& a_numRows,
                 const int& a_numCols);

  void matMul(Real** a_A,
              Real** a_B,
              Real** a_C,
              const int& a_numRowsA,
              const int& a_numColsA,
              const int& a_numColsB);

  void backSolve(Real** a_A,
                 const Vector<Real>& a_rhs,
                 const int& a_numArows,
                 Vector<Real>& a_x);

  void AtimesX(Real** A,
               const Vector<Real>&x,
               const int& numRowsA,
               Vector<Real>& Ax);

  void allocArray(const int& rows,
                  const int& cols,
                  Real**& A);

  void freeArray(const int& rows,
                 const int& cols,
                 Real**& A);

  /// outputs a matrix A
  void output(const int& rows,
              const int& cols,
              Real**& A,
              char* name);
};

#include "NamespaceFooter.H"
#endif
