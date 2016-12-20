#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "GeometryService.H"
#include "GeometryShop.H"
#include "Moments.H"
#include "LSquares.H"
#include "PolyGeom.H"
#include "RealVect.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"


void LSquares::LeastSquares(Real** A,
                            Vector<Real>&x,
                            const Vector<Real>&rhs)

{
  CH_TIMELEAF("LSquares::leastSquares");
  int numColsA = x.size();
  int numRowsA = rhs.size();

  Real** Atrans;
  allocArray(numColsA,numRowsA,Atrans);

  Real** LS; // LS = least squares matrix; Ax = rhs is an overdetermined system
  allocArray(numColsA,numColsA,LS);

  transpose(A,Atrans,numRowsA,numColsA);

  matMul(Atrans,A,LS,numColsA,numRowsA,numColsA);

  // pattern is (A,B,A*B,numrowsA,numcolsA,numcolsB)
  Vector<Real> ATrhs(numColsA);
  AtimesX(Atrans,rhs,numColsA,ATrhs); // numCols = number of rows of Atranspose*A

  gaussElim(LS,ATrhs);
  // char* bug = "LS";
  // output(numColsA,numColsA,LS,bug);

  backSolve(LS,ATrhs,numColsA,x); // x is the answer

  freeArray(numColsA,numRowsA,Atrans);
  freeArray(numColsA,numColsA,LS);
}

void LSquares::AtimesX(Real** A,
                       const Vector<Real>& x,
                       const int& numRowsA,
                       Vector<Real>& Ax)
{
  for (int i = 0; i < numRowsA; i++)
    {
      Real* scanA = A[i];

      Real sum = 0.0;
      int xSize = x.size();

      for (int j = 0; j < xSize; j++)
        {
          sum += *(scanA++) * x[j];
        }

      Ax[i] = sum;
    }
}

// this is written for square A
int  LSquares::gaussElim(Real** A,
                         Vector<Real>& rhs)
{
  //   char* name = "A";
  int currRow = 0;
  int numRows = rhs.size();
  int numCols = rhs.size();

  for (int currCol = 0; currCol < numCols; ++currCol)
    {
      int pivot;

      findPivot(A,currCol, currRow,numRows,pivot);
      // output(numRows,numCols,A,name);
      //
      if (Abs(A[pivot][currCol]) > 1.0e-15)
      {
        swapRows(A,currRow,pivot,numCols);
        swapRows(rhs,currRow,pivot);

        Real Beta = A[currRow][currCol];
        Beta = 1.0/Beta;

        timesBeta(A,currRow,Beta, numCols);
        timesBeta(rhs,currRow,Beta);

        // output(numRows,numCols,A,name);
        for (int rows = currRow; rows < numRows-1; ++rows)
        {
          Real alpha = -A[rows+1][currCol];

          addRows(A,rows+1,alpha,currRow,numCols);
          addRows(rhs,rows+1,alpha,currRow);
          // output(numRows,numCols,A,name);
        }
        currRow += 1;
      }
      else
      {
        // MayDay::Warning("small pivot in gaussElim");
      }
    }

  return 0;
}

void LSquares::swapRows(Vector<Real>& rhs,
                        const int& currRow,
                        const int& pivot)
{
  Real    temp = rhs[currRow];
  rhs[currRow] = rhs[pivot];
  rhs[pivot]   = temp;
}

void LSquares::swapRows(Real** A,
                        const int& rowi,
                        const int& rowj,
                        const int& numCols)
{
  Real *scani = A[rowi];
  Real *scanj = A[rowj];

  for (int count = 0; count < numCols; ++count)
    {
      Real temp;

      temp     = (*scani);
      (*scani) = (*scanj);
      (*scanj) = temp;

      scani++;
      scanj++;
    }
}

int LSquares::findPivot(Real** A,
                        const int& currCol,
                        const int& currRow,
                        const int& numRows,
                        int& pivot)
{
  Real max = 0;
  pivot = currRow;

  for (int count = currRow; count < numRows; count++)
    {
      if (Abs(A[count][currCol]) > max)
        {
          max = Abs(A[count][currCol]);
          pivot = count;
        }
   }

  return 0;
}

void LSquares::addRows(Vector<Real>& rhs,
                       const int& rowi,
                       const Real& alpha,
                       const int& rowj)
{
  rhs[rowi] += alpha * rhs[rowj];
}

void LSquares::addRows(Real** A,
                       const int& rowi,
                       const Real& alpha,
                       const int& rowj,
                       const int& numCols) // rowi += alpha * rowj
{
  Real *scani = A[rowi];
  Real *scanj = A[rowj];

  for (int count = 0; count < numCols; ++count)
    {
      (*scani) += alpha * (*scanj);
      scani++;
      scanj++;
    }
}

void LSquares::timesBeta(Vector<Real>&rhs,
                         const int& currRow,
                         const Real& Beta)
{
  rhs[currRow] *= Beta;
}

void  LSquares::timesBeta(Real** A,
                          const int& rowi,
                          const Real& Beta,
                          const int& numCols)
{
  Real *scanA = A[rowi];

  for (int count = 0; count < numCols; ++count)
    {
      (*scanA) *= Beta;
      scanA++;
    }
}

void LSquares::transpose(Real** a_A,
                         Real** a_Atrans,
                         const int& a_numRowsA,
                         const int& a_numColsA)
{
  for (int irow = 0; irow < a_numColsA; ++irow)
    {
      Real *scanAtrans = a_Atrans[irow];

      for (int icol = 0; icol < a_numRowsA; ++icol)
        {
          (*scanAtrans) = a_A[icol][irow];
          scanAtrans++;
       }
    }
}

void LSquares::matMul(Real** a_A,
                      Real** a_B,
                      Real** a_C,
                      const int& a_numRowsA,
                      const int& a_numColsA,
                      const int& a_numColsB)
{
  for (int i = 0; i < a_numRowsA; ++i)
  {
    Real* scanC = a_C[i];

    for (int j = 0; j < a_numColsB; ++j)
    {
      Real* scanA = a_A[i];
      Real sum = 0.0;

      for (int k = 0; k < a_numColsA; ++k)
      {
        sum += (*scanA) * a_B[k][j];
        scanA++;
      }

      (*scanC) = sum;
      scanC++;
    }
  }
}

void LSquares::backSolve(Real** a_A,
                         const Vector <Real>& a_rhs,
                         const int& a_numArows,
                         Vector<Real>& a_x)
{
  int N = a_numArows;
  for (int n = 0; n < N; ++n)
    {
      a_x[n] = 0.0;
    }

  for (int n = N-1; n >= 0; --n)
    {
      for (int m = 1; m < N-n; ++m)
        {
          a_x[n] -= a_x[n+m] * a_A[n][n+m];
        }
      if (Abs(a_A[n][n]) > 0.0)
        {
          a_x[n] += a_rhs[n] / a_A[n][n]; // this only works for square A
        }
    }
}

void LSquares::allocArray(const int& rows,
                          const int& cols,
                          Real**& A)
{
  A = new Real* [rows];

  for (int i = 0; i < rows;i++)
    {
      A[i] = new Real [cols];
      Real* scanA = A[i];

      for (int j = 0; j < cols; j++)
        {
          *(scanA++) = 0.0;
        }
    }
}

void LSquares::freeArray(const int& rows,
                         const int& cols,
                         Real**& A)
{
  for (int i = 0; i < rows; i++)
    {
      delete[] A[i];
    }

  delete[] A;
}

void LSquares::output(const int& rows,
                      const int& cols,
                      Real**& A,
                      char* name)
{
  // pout() << "outputting " << name << endl;
  for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols; j++)
       {
         pout() << A[i][j] << " ";
       }
      pout() << endl;
    }
}

#include "NamespaceFooter.H"
