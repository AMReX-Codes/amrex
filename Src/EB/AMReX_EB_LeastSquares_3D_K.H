#ifndef AMREX_MLEB_LEASTSQUARES_K_H_
#define AMREX_MLEB_LEASTSQUARES_K_H_
#include <AMReX_Config.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

// This is the 3D version (i.e. for 10x10 (A^T A) matrix)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void decomp_chol_np10(Array2D<Real,0,9,0,9>& aa)
{
    int neq = 10;

    Real p[10];
    Real sum1;
    int ising;

    for (int ii = 0; ii < neq; ii++)
    {
        ising = 0;

        for (int jj = ii; jj < neq; jj++)
        {
            sum1 = aa(ii,jj);

            for (int kk = ii-1; kk >= 0; kk--)
            {
                sum1 = sum1 - aa(ii,kk)*aa(jj,kk);
            }

            if (ii == jj)
            {
                 if (sum1 <= Real(0.0))
                 {
                     p[ii] = Real(0.0);
                     ising = 1;
                 } else {
                     p[ii] = std::sqrt(sum1);
                 }
            } else {
                if (ising == 0)
                   aa(jj,ii) = sum1 / p[ii];
                else
                   aa(jj,ii) = Real(0.0);
            }
        }
    }

    for (int ii = 0; ii < neq; ii++)
    {
        for (int jj = ii+1; jj < neq; jj++)
        {
           aa(ii,jj) = Real(0.0);       // Zero upper triangle
           aa(ii,jj) = aa(jj,ii);  // Zero upper triangle
        }

        aa(ii,ii) = p[ii];
    }
}

// This is the 3D version (i.e. for 10x10 (A^T A) matrix)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void  cholsol_np10(Array2D<Real,0,35,0,9>& Amatrix, Array1D<Real,0,9>& b)
{
    int neq = 10;

    Array2D<Real,0,9,0,9> AtA;

    for (int irow = 0; irow < neq; irow++)
        for (int icol = 0; icol < neq; icol++)
            AtA(irow,icol) =  Real(0.0);

    for (int irow = 0; irow < 36; irow++)
    {
         AtA(0,0) += Amatrix(irow,0)*Amatrix(irow,0); // e^T e
         AtA(0,1) += Amatrix(irow,0)*Amatrix(irow,1); // e^T x
         AtA(0,2) += Amatrix(irow,0)*Amatrix(irow,2); // e^T y
         AtA(0,3) += Amatrix(irow,0)*Amatrix(irow,3); // e^T z
         AtA(0,4) += Amatrix(irow,0)*Amatrix(irow,4); // e^T x^2
         AtA(0,5) += Amatrix(irow,0)*Amatrix(irow,5); // e^T x*y
         AtA(0,6) += Amatrix(irow,0)*Amatrix(irow,6); // e^T y^2
         AtA(0,7) += Amatrix(irow,0)*Amatrix(irow,7); // e^T x*z
         AtA(0,8) += Amatrix(irow,0)*Amatrix(irow,8); // e^T y*z
         AtA(0,9) += Amatrix(irow,0)*Amatrix(irow,9); // e^T z^2

         AtA(1,1) += Amatrix(irow,1)*Amatrix(irow,1); // x^T x
         AtA(1,2) += Amatrix(irow,1)*Amatrix(irow,2); // x^T y
         AtA(1,3) += Amatrix(irow,1)*Amatrix(irow,3); // x^T y
         AtA(1,4) += Amatrix(irow,1)*Amatrix(irow,4); // x^T (x^2)
         AtA(1,5) += Amatrix(irow,1)*Amatrix(irow,5); // x^T (xy)
         AtA(1,6) += Amatrix(irow,1)*Amatrix(irow,6); // x^T (y^2)
         AtA(1,7) += Amatrix(irow,1)*Amatrix(irow,7); // x^T x*z
         AtA(1,8) += Amatrix(irow,1)*Amatrix(irow,8); // x^T y*z
         AtA(1,9) += Amatrix(irow,1)*Amatrix(irow,9); // x^T z^2

         AtA(2,2) += Amatrix(irow,2)*Amatrix(irow,2); // y^T y
         AtA(2,3) += Amatrix(irow,2)*Amatrix(irow,3); // y^T z
         AtA(2,4) += Amatrix(irow,2)*Amatrix(irow,4); // y^T (x^2)
         AtA(2,5) += Amatrix(irow,2)*Amatrix(irow,5); // y^T (xy)
         AtA(2,6) += Amatrix(irow,2)*Amatrix(irow,6); // y^T (y^2)
         AtA(2,7) += Amatrix(irow,2)*Amatrix(irow,7); // y^T x*z
         AtA(2,8) += Amatrix(irow,2)*Amatrix(irow,8); // y^T y*z
         AtA(2,9) += Amatrix(irow,2)*Amatrix(irow,9); // y^T z^2

         AtA(3,3) += Amatrix(irow,3)*Amatrix(irow,3); // z^T z
         AtA(3,4) += Amatrix(irow,3)*Amatrix(irow,4); // z^T (x^2)
         AtA(3,5) += Amatrix(irow,3)*Amatrix(irow,5); // z^T (xy)
         AtA(3,6) += Amatrix(irow,3)*Amatrix(irow,6); // z^T (y^2)
         AtA(3,7) += Amatrix(irow,3)*Amatrix(irow,7); // z^T x*z
         AtA(3,8) += Amatrix(irow,3)*Amatrix(irow,8); // z^T y*z
         AtA(3,9) += Amatrix(irow,3)*Amatrix(irow,9); // z^T z^2

         AtA(4,4) += Amatrix(irow,4)*Amatrix(irow,4); // (x^2)^T (x^2)
         AtA(4,5) += Amatrix(irow,4)*Amatrix(irow,5); // (x^2)^T (xy)
         AtA(4,6) += Amatrix(irow,4)*Amatrix(irow,6); // (x^2)^T (y^2)
         AtA(4,7) += Amatrix(irow,4)*Amatrix(irow,7); // (x^2)^T x*z
         AtA(4,8) += Amatrix(irow,4)*Amatrix(irow,8); // (x^2)^T y*z
         AtA(4,9) += Amatrix(irow,4)*Amatrix(irow,9); // (x^2)^T z^2

         AtA(5,5) += Amatrix(irow,5)*Amatrix(irow,5); // (xy)^T (xy)
         AtA(5,6) += Amatrix(irow,5)*Amatrix(irow,6); // (xy)^T (y^2)
         AtA(5,7) += Amatrix(irow,5)*Amatrix(irow,7); // (xy)^T x*z
         AtA(5,8) += Amatrix(irow,5)*Amatrix(irow,8); // (xy)^T y*z
         AtA(5,9) += Amatrix(irow,5)*Amatrix(irow,9); // (xy)^T z^2

         AtA(6,6) += Amatrix(irow,6)*Amatrix(irow,6); // (y^2)^T (y^2)
         AtA(6,7) += Amatrix(irow,6)*Amatrix(irow,7); // (y^2)^T x*z
         AtA(6,8) += Amatrix(irow,6)*Amatrix(irow,8); // (y^2)^T y*z
         AtA(6,9) += Amatrix(irow,6)*Amatrix(irow,9); // (y^2)^T z^2

         AtA(7,7) += Amatrix(irow,7)*Amatrix(irow,7); // (xz)^T x*z
         AtA(7,8) += Amatrix(irow,7)*Amatrix(irow,8); // (xz)^T y*z
         AtA(7,9) += Amatrix(irow,7)*Amatrix(irow,9); // (xz)^T z^2

         AtA(8,8) += Amatrix(irow,8)*Amatrix(irow,8); // (yz)^T y*z
         AtA(8,9) += Amatrix(irow,8)*Amatrix(irow,9); // (yz)^T z^2

         AtA(9,9) += Amatrix(irow,9)*Amatrix(irow,9); // (z^2)^T z^2
    }

    for (int irow = 0; irow < neq-1; irow++)
        for (int icol = irow+1; icol < neq; icol++)
           AtA(icol,irow) = AtA(irow,icol);

    decomp_chol_np10(AtA);

    if (AtA(0,0) > 0.)
       b(0) = b(0) / AtA(0,0);
    else
       b(0) = 0.;

    for (int ii = 1; ii < neq; ii++)
    {
       if (AtA(ii,ii) > 0.)
       {
          for (int jj = 0; jj < ii; jj++)
            b(ii) = b(ii) - AtA(ii,jj)*b(jj);

          b(ii) = b(ii) / AtA(ii,ii);
       }
       else
       {
          b(ii) = Real(0.0);
       }
    }

    if (AtA(neq-1,neq-1) > 0.)
       b(neq-1) = b(neq-1) / AtA(neq-1,neq-1);
    else
       b(neq-1) = Real(0.0);

    for (int ii = neq-2; ii >= 0; ii--)
    {
       if (AtA(ii,ii) > 0.)
       {
          for (int jj = ii+1; jj < neq; jj++)
            b(ii) = b(ii) - AtA(ii,jj)*b(jj);

          b(ii) = b(ii) / AtA(ii,ii);
       }
       else
       {
          b(ii) = Real(0.0);
       }
    }
}

// This is the 3D version (i.e. for 10x10 (A^T A) matrix) but for A = 54 x 10
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void  cholsol_for_eb(Array2D<Real,0,53,0,9>& Amatrix, Array1D<Real,0,9>& b)
{
    int neq = 10;

    Array2D<Real,0,9,0,9> AtA;

    for (int irow = 0; irow < neq; irow++)
        for (int icol = 0; icol < neq; icol++)
            AtA(irow,icol) =  Real(0.0);

    for (int irow = 0; irow < 54; irow++)
    {
         AtA(0,0) += Amatrix(irow,0)*Amatrix(irow,0); // e^T e
         AtA(0,1) += Amatrix(irow,0)*Amatrix(irow,1); // e^T x
         AtA(0,2) += Amatrix(irow,0)*Amatrix(irow,2); // e^T y
         AtA(0,3) += Amatrix(irow,0)*Amatrix(irow,3); // e^T z
         AtA(0,4) += Amatrix(irow,0)*Amatrix(irow,4); // e^T x^2
         AtA(0,5) += Amatrix(irow,0)*Amatrix(irow,5); // e^T x*y
         AtA(0,6) += Amatrix(irow,0)*Amatrix(irow,6); // e^T y^2
         AtA(0,7) += Amatrix(irow,0)*Amatrix(irow,7); // e^T x*z
         AtA(0,8) += Amatrix(irow,0)*Amatrix(irow,8); // e^T y*z
         AtA(0,9) += Amatrix(irow,0)*Amatrix(irow,9); // e^T z^2

         AtA(1,1) += Amatrix(irow,1)*Amatrix(irow,1); // x^T x
         AtA(1,2) += Amatrix(irow,1)*Amatrix(irow,2); // x^T y
         AtA(1,3) += Amatrix(irow,1)*Amatrix(irow,3); // x^T y
         AtA(1,4) += Amatrix(irow,1)*Amatrix(irow,4); // x^T (x^2)
         AtA(1,5) += Amatrix(irow,1)*Amatrix(irow,5); // x^T (xy)
         AtA(1,6) += Amatrix(irow,1)*Amatrix(irow,6); // x^T (y^2)
         AtA(1,7) += Amatrix(irow,1)*Amatrix(irow,7); // x^T x*z
         AtA(1,8) += Amatrix(irow,1)*Amatrix(irow,8); // x^T y*z
         AtA(1,9) += Amatrix(irow,1)*Amatrix(irow,9); // x^T z^2

         AtA(2,2) += Amatrix(irow,2)*Amatrix(irow,2); // y^T y
         AtA(2,3) += Amatrix(irow,2)*Amatrix(irow,3); // y^T z
         AtA(2,4) += Amatrix(irow,2)*Amatrix(irow,4); // y^T (x^2)
         AtA(2,5) += Amatrix(irow,2)*Amatrix(irow,5); // y^T (xy)
         AtA(2,6) += Amatrix(irow,2)*Amatrix(irow,6); // y^T (y^2)
         AtA(2,7) += Amatrix(irow,2)*Amatrix(irow,7); // y^T x*z
         AtA(2,8) += Amatrix(irow,2)*Amatrix(irow,8); // y^T y*z
         AtA(2,9) += Amatrix(irow,2)*Amatrix(irow,9); // y^T z^2

         AtA(3,3) += Amatrix(irow,3)*Amatrix(irow,3); // z^T z
         AtA(3,4) += Amatrix(irow,3)*Amatrix(irow,4); // z^T (x^2)
         AtA(3,5) += Amatrix(irow,3)*Amatrix(irow,5); // z^T (xy)
         AtA(3,6) += Amatrix(irow,3)*Amatrix(irow,6); // z^T (y^2)
         AtA(3,7) += Amatrix(irow,3)*Amatrix(irow,7); // z^T x*z
         AtA(3,8) += Amatrix(irow,3)*Amatrix(irow,8); // z^T y*z
         AtA(3,9) += Amatrix(irow,3)*Amatrix(irow,9); // z^T z^2

         AtA(4,4) += Amatrix(irow,4)*Amatrix(irow,4); // (x^2)^T (x^2)
         AtA(4,5) += Amatrix(irow,4)*Amatrix(irow,5); // (x^2)^T (xy)
         AtA(4,6) += Amatrix(irow,4)*Amatrix(irow,6); // (x^2)^T (y^2)
         AtA(4,7) += Amatrix(irow,4)*Amatrix(irow,7); // (x^2)^T x*z
         AtA(4,8) += Amatrix(irow,4)*Amatrix(irow,8); // (x^2)^T y*z
         AtA(4,9) += Amatrix(irow,4)*Amatrix(irow,9); // (x^2)^T z^2

         AtA(5,5) += Amatrix(irow,5)*Amatrix(irow,5); // (xy)^T (xy)
         AtA(5,6) += Amatrix(irow,5)*Amatrix(irow,6); // (xy)^T (y^2)
         AtA(5,7) += Amatrix(irow,5)*Amatrix(irow,7); // (xy)^T x*z
         AtA(5,8) += Amatrix(irow,5)*Amatrix(irow,8); // (xy)^T y*z
         AtA(5,9) += Amatrix(irow,5)*Amatrix(irow,9); // (xy)^T z^2

         AtA(6,6) += Amatrix(irow,6)*Amatrix(irow,6); // (y^2)^T (y^2)
         AtA(6,7) += Amatrix(irow,6)*Amatrix(irow,7); // (y^2)^T x*z
         AtA(6,8) += Amatrix(irow,6)*Amatrix(irow,8); // (y^2)^T y*z
         AtA(6,9) += Amatrix(irow,6)*Amatrix(irow,9); // (y^2)^T z^2

         AtA(7,7) += Amatrix(irow,7)*Amatrix(irow,7); // (xz)^T x*z
         AtA(7,8) += Amatrix(irow,7)*Amatrix(irow,8); // (xz)^T y*z
         AtA(7,9) += Amatrix(irow,7)*Amatrix(irow,9); // (xz)^T z^2

         AtA(8,8) += Amatrix(irow,8)*Amatrix(irow,8); // (yz)^T y*z
         AtA(8,9) += Amatrix(irow,8)*Amatrix(irow,9); // (yz)^T z^2

         AtA(9,9) += Amatrix(irow,9)*Amatrix(irow,9); // (z^2)^T z^2
    }

    for (int irow = 0; irow < neq-1; irow++)
        for (int icol = irow+1; icol < neq; icol++)
           AtA(icol,irow) = AtA(irow,icol);

    decomp_chol_np10(AtA);

    if (AtA(0,0) > 0.)
       b(0) = b(0) / AtA(0,0);
    else
       b(0) = 0.;

    for (int ii = 1; ii < neq; ii++)
    {
       if (AtA(ii,ii) > 0.)
       {
          for (int jj = 0; jj < ii; jj++)
            b(ii) = b(ii) - AtA(ii,jj)*b(jj);

          b(ii) = b(ii) / AtA(ii,ii);
       }
       else
       {
          b(ii) = Real(0.0);
       }
    }

    if (AtA(neq-1,neq-1) > 0.)
       b(neq-1) = b(neq-1) / AtA(neq-1,neq-1);
    else
       b(neq-1) = Real(0.0);

    for (int ii = neq-2; ii >= 0; ii--)
    {
       if (AtA(ii,ii) > 0.)
       {
          for (int jj = ii+1; jj < neq; jj++)
            b(ii) = b(ii) - AtA(ii,jj)*b(jj);

          b(ii) = b(ii) / AtA(ii,ii);
       }
       else
       {
          b(ii) = Real(0.0);
       }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_x_of_phi_on_centroids(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Real& yloc_on_xface, Real& zloc_on_xface,
                                bool is_eb_dirichlet, bool is_eb_inhomog)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,35,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 9 are cell centroids, second 9 are EB centroids

    for (int irow = 0; irow < 36; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    //  Columns 0-3: [e x y z]
    for (int ii = i-1; ii <= i; ii++) // Normal to face
      for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
         for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (jj-(j-1)) + 3*(kk-(k-1)) + 9*(ii-(i-1));

                Real x_off = static_cast<Real>(ii-i) + Real(0.5);
                Real y_off = static_cast<Real>(jj-j);
                Real z_off = static_cast<Real>(kk-k);

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0);
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1) - yloc_on_xface;
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2) - zloc_on_xface;

                if (!flag(ii,jj,kk).isRegular())
                {
                    Amatrix(a_ind+18,0) = Real(1.0);
                    Amatrix(a_ind+18,1) = x_off + bcent(ii,jj,kk,0);
                    Amatrix(a_ind+18,2) = y_off + bcent(ii,jj,kk,1) - yloc_on_xface;
                    Amatrix(a_ind+18,3) = z_off + bcent(ii,jj,kk,2) - zloc_on_xface;
                }
            }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 36; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int ii = i-1; ii <= i; ii++)  // Normal to face
          for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
           for (int jj = j-1; jj <= j+1; jj++) // Tangential to face
               if (!flag(ii,jj,kk).isCovered())
               {
                    int a_ind  = (jj-(j-1)) + 3*(kk-(k-1)) + 9*(ii-(i-1));

                    rhs(irow) += Amatrix(a_ind,irow)*  phi(ii,jj,kk,n);

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog)
                        rhs(irow) += Amatrix(a_ind+18,irow)*phieb(ii,jj,kk,n);
               }
    }

    cholsol_np10(Amatrix, rhs);

    return rhs(1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_y_of_phi_on_centroids(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Real& xloc_on_yface, Real& zloc_on_yface,
                                bool is_eb_dirichlet, bool is_eb_inhomog)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,35,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 9 are cell centroids, second 9 are EB centroids

    for (int irow = 0; irow < 36; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    //  Columns 0-2: [e x y]
    for (int jj = j-1; jj <= j; jj++) // Normal to face
      for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
        for (int ii = i-1; ii <= i+1; ii++)  // Tangential to face
            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(kk-(k-1)) + 9*(jj-(j-1));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j) + Real(0.5);
                Real z_off = static_cast<Real>(kk-k);

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0) - xloc_on_yface;
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1);
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2) - zloc_on_yface;

                if (!flag(ii,jj,kk).isRegular())
                {
                    Amatrix(a_ind+18,0) = Real(1.0);
                    Amatrix(a_ind+18,1) = x_off + bcent(ii,jj,kk,0) - xloc_on_yface;
                    Amatrix(a_ind+18,2) = y_off + bcent(ii,jj,kk,1);
                    Amatrix(a_ind+18,3) = z_off + bcent(ii,jj,kk,2) - zloc_on_yface;
                }
            }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 36; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int jj = j-1; jj <= j; jj++) // Normal to face
           for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
             for (int ii = i-1; ii <= i+1; ii++) // Tangential to face
                if (!flag(ii,jj,kk).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(kk-(k-1)) + 9*(jj-(j-1));

                    rhs(irow) += Amatrix(a_ind,irow)*  phi(ii,jj,kk,n);

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog)
                        rhs(irow) += Amatrix(a_ind+18,irow)*phieb(ii,jj,kk,n);
                }
    }

    cholsol_np10(Amatrix, rhs);

    return rhs(2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_z_of_phi_on_centroids(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Real& xloc_on_zface, Real& yloc_on_zface,
                                bool is_eb_dirichlet, bool is_eb_inhomog)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,35,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 9 are cell centroids, second 9 are EB centroids

    for (int irow = 0; irow < 36; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    //  Columns 0-3: [e x y z]
    for (int kk = k-1; kk <= k; kk++) // Normal to face
    {
        for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
         for (int ii = i-1; ii <= i+1; ii++)  // Tangential to face
         {
            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j);
                Real z_off = static_cast<Real>(kk-k) + Real(0.5);

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0) - xloc_on_zface;
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1) - yloc_on_zface;
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2);

                if (!flag(ii,jj,kk).isRegular())
                {
                    Amatrix(a_ind+18,0) = Real(1.0);
                    Amatrix(a_ind+18,1) = x_off + bcent(ii,jj,kk,0) - xloc_on_zface;
                    Amatrix(a_ind+18,2) = y_off + bcent(ii,jj,kk,1) - yloc_on_zface;
                    Amatrix(a_ind+18,3) = z_off + bcent(ii,jj,kk,2);
                }
            }
         }
    }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 36; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int kk = k-1; kk <= k; kk++) // Normal to face
          for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
             for (int ii = i-1; ii <= i+1; ii++) // Tangential to face
                if (!flag(ii,jj,kk).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

                    rhs(irow) += Amatrix(a_ind,irow)*  phi(ii,jj,kk,n);

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog)
                        rhs(irow) += Amatrix(a_ind+18,irow)*phieb(ii,jj,kk,n);
                }
    }

    cholsol_np10(Amatrix, rhs);

    return rhs(3);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_eb_of_phi_on_centroids(int i,int j,int k,int n,
                                 Array4<Real const> const& phi,
                                 Array4<Real const> const& phieb,
                                 Array4<EBCellFlag const> const& flag,
                                 Array4<Real const> const& ccent,
                                 Array4<Real const> const& bcent,
                                 Real& nrmx, Real& nrmy, Real& nrmz,
                                 bool is_eb_inhomog)
{
    Array2D<Real,0,53,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 27 are cell centroids, second 27 are EB centroids

    for (int irow = 0; irow < 54; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    //  Columns 0-3: [e x y z]
    for (int kk = k-1; kk <= k+1; kk++)
      for (int jj = j-1; jj <= j+1; jj++)
         for (int ii = i-1; ii <= i+1; ii++)
         {
            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

                Real x_off = static_cast<Real>(ii-i) - bcent(i,j,k,0);
                Real y_off = static_cast<Real>(jj-j) - bcent(i,j,k,1);
                Real z_off = static_cast<Real>(kk-k) - bcent(i,j,k,2);

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0);
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1);
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2);

                if (!flag(ii,jj,kk).isRegular())
                {
                    Amatrix(a_ind+27,0) = Real(1.0);
                    Amatrix(a_ind+27,1) = x_off + bcent(ii,jj,kk,0);
                    Amatrix(a_ind+27,2) = y_off + bcent(ii,jj,kk,1);
                    Amatrix(a_ind+27,3) = z_off + bcent(ii,jj,kk,2);
                }
            }
         }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 54; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int kk = k-1; kk <= k+1; kk++)
          for (int jj = j-1; jj <= j+1; jj++)
             for (int ii = i-1; ii <= i+1; ii++)
                if (!flag(ii,jj,kk).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

                    rhs(irow) += Amatrix(a_ind,irow)*  phi(ii,jj,kk,n);

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog)
                        rhs(irow) += Amatrix(a_ind+27,irow)*phieb(ii,jj,kk,n);
                }
    }

    cholsol_for_eb(Amatrix, rhs);

    Real dphidn = rhs(1)*nrmx + rhs(2)*nrmy + rhs(3)*nrmz;
    return dphidn;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_x_of_phi_on_centroids_extdir(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Array4<Real const> const& vfrac,
                                Real& yloc_on_xface, Real& zloc_on_xface,
                                bool is_eb_dirichlet, bool is_eb_inhomog,
                                const bool on_x_face, const int domlo_x, const int domhi_x,
                                const bool on_y_face, const int domlo_y, const int domhi_y,
                                const bool on_z_face, const int domlo_z, const int domhi_z)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);

    Array2D<Real,0,53,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 9 are cell centroids, second 9 are EB centroids

    for (int irow = 0; irow < 54; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    const int im = (i > domhi_x) ? 2 : 1;
    const int ip = 2 - im;

    //  Columns 0-3: [e x y z]
    for (int ii = i-im; ii <= i+ip; ii++) // Normal to face
      for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
         for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
         {
            // Don't include corner cells. Could make this even more strict
            // by removing the on_?_face restrictions.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii < domlo_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk > domhi_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk > domhi_z))) {
              continue;
            }

            if ( !phi.contains(ii,jj,kk) )
               continue;

            if (!flag(ii,jj,kk).isCovered())
            {

                int a_ind  = (jj-(j-1)) + 3*(kk-(k-1)) + 9*(ii-(i-im));

                Real x_off = static_cast<Real>(ii-i) + Real(0.5);
                Real y_off = static_cast<Real>(jj-j);
                Real z_off = static_cast<Real>(kk-k);

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,kk) != Real(1.0) || vfrac(ii+2,jj,kk) != Real(1.0)) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,kk) != Real(1.0) || vfrac(ii-2,jj,kk) != Real(1.0)))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,kk) != Real(1.0) || vfrac(ii,jj+2,kk) != Real(1.0)) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,kk) != Real(1.0) || vfrac(ii,jj-2,kk) != Real(1.0)) )
                        continue;
                }

                if(on_z_face){
                    if (kk < domlo_z && (vfrac(ii,jj,kk+1) != Real(1.0) || vfrac(ii,jj,kk+2) != Real(1.0)) )
                        continue;
                    if (kk > domhi_z && (vfrac(ii,jj,kk-1) != Real(1.0) || vfrac(ii,jj,kk-2) != Real(1.0)) )
                        continue;
                }

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0);
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1) - yloc_on_xface;
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2) - zloc_on_xface;

                // Add in information about the location of the EB. Exclude
                // EBs that are outside the domain.
                if (flag(ii,jj,kk).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x &&
                    domlo_y <= jj && jj <= domhi_y &&
                    domlo_z <= kk && kk <= domhi_z)
                {
                    Amatrix(a_ind+27,0) = Real(1.0);
                    Amatrix(a_ind+27,1) = x_off + bcent(ii,jj,kk,0);
                    Amatrix(a_ind+27,2) = y_off + bcent(ii,jj,kk,1) - yloc_on_xface;
                    Amatrix(a_ind+27,3) = z_off + bcent(ii,jj,kk,2) - zloc_on_xface;
                }
            }
         }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 54; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int ii = i-im; ii <= i+ip; ii++)  // Normal to face
          for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
           for (int jj = j-1; jj <= j+1; jj++) { // Tangential to face
               if ( !phi.contains(ii,jj,kk) )
                  continue;

               if (!flag(ii,jj,kk).isCovered())
               {
                    int a_ind  = (jj-(j-1)) + 3*(kk-(k-1)) + 9*(ii-(i-im));
                    Real phi_val = Amatrix(a_ind,0)*phi(ii,jj,kk,n);

                    rhs(irow) += Amatrix(a_ind,irow)*  phi_val;

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog && Amatrix(a_ind+27,irow) != Real(0.0))
                        rhs(irow) += Amatrix(a_ind+27,irow)*phieb(ii,jj,kk,n);
               }
            }
    }

    cholsol_for_eb(Amatrix, rhs);

    return rhs(1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_y_of_phi_on_centroids_extdir(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Array4<Real const> const& vfrac,
                                Real& xloc_on_yface, Real& zloc_on_yface,
                                bool is_eb_dirichlet, bool is_eb_inhomog,
                                const bool on_x_face, const int domlo_x, const int domhi_x,
                                const bool on_y_face, const int domlo_y, const int domhi_y,
                                const bool on_z_face, const int domlo_z, const int domhi_z)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,53,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 9 are cell centroids, second 9 are EB centroids

    for (int irow = 0; irow < 54; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    const int jm = (j > domhi_y) ? 2 : 1;
    const int jp = 2 - jm;

    //  Columns 0-2: [e x y]
    for (int jj = j-jm; jj <= j+jp; jj++) // Normal to face
      for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
        for (int ii = i-1; ii <= i+1; ii++)  // Tangential to face
        {
            // Don't include corner cells. Could make this even more strict
            // by removing the on_?_face restrictions.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii < domlo_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk > domhi_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk > domhi_z))) {
              continue;
            }

            if ( !phi.contains(ii,jj,kk) )
               continue;

            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(kk-(k-1)) + 9*(jj-(j-jm));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j) + Real(0.5);
                Real z_off = static_cast<Real>(kk-k);

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,kk) != Real(1.0) || vfrac(ii+2,jj,kk) != Real(1.0)) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,kk) != Real(1.0) || vfrac(ii-2,jj,kk) != Real(1.0)))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,kk) != Real(1.0) || vfrac(ii,jj+2,kk) != Real(1.0)) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,kk) != Real(1.0) || vfrac(ii,jj-2,kk) != Real(1.0)) )
                        continue;
                }

                if(on_z_face){
                    if (kk < domlo_z && (vfrac(ii,jj,kk+1) != Real(1.0) || vfrac(ii,jj,kk+2) != Real(1.0)) )
                        continue;
                    if (kk > domhi_z && (vfrac(ii,jj,kk-1) != Real(1.0) || vfrac(ii,jj,kk-2) != Real(1.0)) )
                        continue;
                }

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0) - xloc_on_yface;
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1);
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2) - zloc_on_yface;

                // Add in information about the location of the EB. Exclude
                // EBs that are outside the domain.
                if (flag(ii,jj,kk).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x &&
                    domlo_y <= jj && jj <= domhi_y &&
                    domlo_z <= kk && kk <= domhi_z)
                {
                    Amatrix(a_ind+27,0) = Real(1.0);
                    Amatrix(a_ind+27,1) = x_off + bcent(ii,jj,kk,0) - xloc_on_yface;
                    Amatrix(a_ind+27,2) = y_off + bcent(ii,jj,kk,1);
                    Amatrix(a_ind+27,3) = z_off + bcent(ii,jj,kk,2) - zloc_on_yface;
                }
            }
        }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 54; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int jj = j-jm; jj <= j+jp; jj++) // Normal to face
           for (int kk = k-1; kk <= k+1; kk++)  // Tangential to face
             for (int ii = i-1; ii <= i+1; ii++) {// Tangential to face
                if ( !phi.contains(ii,jj,kk) )
                   continue;

                if (!flag(ii,jj,kk).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(kk-(k-1)) + 9*(jj-(j-jm));
                    Real phi_val = Amatrix(a_ind,0) * phi(ii,jj,kk,n);

                    rhs(irow) += Amatrix(a_ind,irow)*  phi_val;

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog && Amatrix(a_ind+27,irow) != Real(0.0))
                        rhs(irow) += Amatrix(a_ind+27,irow)*phieb(ii,jj,kk,n);
                }
             }
    }

    cholsol_for_eb(Amatrix, rhs);

    return rhs(2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_z_of_phi_on_centroids_extdir(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Array4<Real const> const& vfrac,
                                Real& xloc_on_zface, Real& yloc_on_zface,
                                bool is_eb_dirichlet, bool is_eb_inhomog,
                                const bool on_x_face, const int domlo_x, const int domhi_x,
                                const bool on_y_face, const int domlo_y, const int domhi_y,
                                const bool on_z_face, const int domlo_z, const int domhi_z)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,53,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 9 are cell centroids, second 9 are EB centroids

    for (int irow = 0; irow < 54; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    const int km = (k > domhi_z) ? 2 : 1;
    const int kp = 2 - km;

    //  Columns 0-3: [e x y z]
    for (int kk = k-km; kk <= k+kp; kk++) // Normal to face
    {
        for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
         for (int ii = i-1; ii <= i+1; ii++)  // Tangential to face
         {
            // Don't include corner cells. Could make this even more strict
            // by removing the on_?_face restrictions.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii < domlo_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk > domhi_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk > domhi_z))) {
              continue;
            }

            if (!phi.contains(ii,jj,kk))
               continue;

            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-km));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j);
                Real z_off = static_cast<Real>(kk-k) + Real(0.5);

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,kk) != Real(1.0) || vfrac(ii+2,jj,kk) != Real(1.0)) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,kk) != Real(1.0) || vfrac(ii-2,jj,kk) != Real(1.0)))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,kk) != Real(1.0) || vfrac(ii,jj+2,kk) != Real(1.0)) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,kk) != Real(1.0) || vfrac(ii,jj-2,kk) != Real(1.0)) )
                        continue;
                }

                if(on_z_face){
                    if (kk < domlo_z && (vfrac(ii,jj,kk+1) != Real(1.0) || vfrac(ii,jj,kk+2) != Real(1.0)) )
                        continue;
                    if (kk > domhi_z && (vfrac(ii,jj,kk-1) != Real(1.0) || vfrac(ii,jj,kk-2) != Real(1.0)) )
                        continue;
                }

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0) - xloc_on_zface;
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1) - yloc_on_zface ;
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2);

                // Add in information about the location of the EB. Exclude
                // EBs that are outside the domain.
                if (flag(ii,jj,kk).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x &&
                    domlo_y <= jj && jj <= domhi_y &&
                    domlo_z <= kk && kk <= domhi_z)
                {
                    Amatrix(a_ind+27,0) = Real(1.0);
                    Amatrix(a_ind+27,1) = x_off + bcent(ii,jj,kk,0) - xloc_on_zface;
                    Amatrix(a_ind+27,2) = y_off + bcent(ii,jj,kk,1) - yloc_on_zface;
                    Amatrix(a_ind+27,3) = z_off + bcent(ii,jj,kk,2);
                }
            }
         }
    }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 54; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int kk = k-km; kk <= k+kp; kk++) // Normal to face
          for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
             for (int ii = i-1; ii <= i+1; ii++) {// Tangential to face
                if ( !phi.contains(ii,jj,kk) )
                   continue;

                if (!flag(ii,jj,kk).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-km));
                    Real phi_val = Amatrix(a_ind,0) * phi(ii,jj,kk,n);

                    rhs(irow) += Amatrix(a_ind,irow)*  phi_val;

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog && Amatrix(a_ind+27,irow) != Real(0.0))
                        rhs(irow) += Amatrix(a_ind+27,irow)*phieb(ii,jj,kk,n);
                }
             }
    }

    cholsol_for_eb(Amatrix, rhs);

    return rhs(3);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_eb_of_phi_on_centroids_extdir(int i,int j,int k,int n,
                                 Array4<Real const> const& phi,
                                 Array4<Real const> const& phieb,
                                 Array4<EBCellFlag const> const& flag,
                                 Array4<Real const> const& ccent,
                                 Array4<Real const> const& bcent,
                                 Array4<Real const> const& vfrac,
                                 Real& nrmx, Real& nrmy, Real& nrmz,
                                 bool is_eb_inhomog,
                                 const bool on_x_face, const int domlo_x, const int domhi_x,
                                 const bool on_y_face, const int domlo_y, const int domhi_y,
                                 const bool on_z_face, const int domlo_z, const int domhi_z)
{
    Array2D<Real,0,53,0,9> Amatrix;
    Array1D<Real,0, 9    > rhs;

    // Order of column -- first 27 are cell centroids, second 27 are EB centroids

    for (int irow = 0; irow < 54; irow++)
        for (int icol = 0; icol < 10; icol++)
            Amatrix(irow,icol) =  Real(0.0);

    //  Columns 0-3: [e x y z]
    for (int kk = k-1; kk <= k+1; kk++)
      for (int jj = j-1; jj <= j+1; jj++)
         for (int ii = i-1; ii <= i+1; ii++)
         {
            // This is likely overkill for EB grads.
            // Don't include corner cells. Could make this even more strict
            // by removing the on_?_face restrictions.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii < domlo_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj < domlo_y) && (on_z_face && kk > domhi_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk < domlo_z)) ||
                ((on_y_face && jj > domhi_y) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii < domlo_x) && (on_z_face && kk > domhi_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk < domlo_z)) ||
                ((on_x_face && ii > domhi_x) && (on_z_face && kk > domhi_z))) {
              continue;
            }

            if ( !phi.contains(ii,jj,kk) )
               continue;

            if (!flag(ii,jj,kk).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

                Real x_off = static_cast<Real>(ii-i) - bcent(i,j,k,0);
                Real y_off = static_cast<Real>(jj-j) - bcent(i,j,k,1);
                Real z_off = static_cast<Real>(kk-k) - bcent(i,j,k,2);

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,kk) != Real(1.0) || vfrac(ii+2,jj,kk) != Real(1.0)) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,kk) != Real(1.0) || vfrac(ii-2,jj,kk) != Real(1.0)))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,kk) != Real(1.0) || vfrac(ii,jj+2,kk) != Real(1.0)) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,kk) != Real(1.0) || vfrac(ii,jj-2,kk) != Real(1.0)) )
                        continue;
                }

                if(on_z_face){
                    if (kk < domlo_z && (vfrac(ii,jj,kk+1) != Real(1.0) || vfrac(ii,jj,kk+2) != Real(1.0)) )
                        continue;
                    if (kk > domhi_z && (vfrac(ii,jj,kk-1) != Real(1.0) || vfrac(ii,jj,kk-2) != Real(1.0)) )
                        continue;
                }

                Amatrix(a_ind,0) = Real(1.0);
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,kk,0);
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,kk,1);
                Amatrix(a_ind,3) = z_off + ccent(ii,jj,kk,2);

                if (flag(ii,jj,kk).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x &&
                    domlo_y <= jj && jj <= domhi_y &&
                    domlo_z <= kk && kk <= domhi_z)
                {
                    Amatrix(a_ind+27,0) = Real(1.0);
                    Amatrix(a_ind+27,1) = x_off + bcent(ii,jj,kk,0);
                    Amatrix(a_ind+27,2) = y_off + bcent(ii,jj,kk,1);
                    Amatrix(a_ind+27,3) = z_off + bcent(ii,jj,kk,2);
                }
            }
         }

    // Columns 4-9 : [x*x  x*y  y*y  x*z  y*z  z*z]
    for (int irow = 0; irow < 54; irow++)
    {
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,5) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,6) =  Amatrix(irow,2) * Amatrix(irow,2);
       Amatrix(irow,7) =  Amatrix(irow,1) * Amatrix(irow,3);
       Amatrix(irow,8) =  Amatrix(irow,2) * Amatrix(irow,3);
       Amatrix(irow,9) =  Amatrix(irow,3) * Amatrix(irow,3);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 10; irow++)
    {
        rhs(irow) = 0.;

        for (int kk = k-1; kk <= k+1; kk++)
          for (int jj = j-1; jj <= j+1; jj++)
             for (int ii = i-1; ii <= i+1; ii++) {
                if ( !phi.contains(ii,jj,kk) )
                  continue;

                if (!flag(ii,jj,kk).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(jj-(j-1)) + 9*(kk-(k-1));

                    Real phi_val = Amatrix(a_ind,0) * phi(ii,jj,kk,n);

                    rhs(irow) += Amatrix(a_ind,irow)*  phi_val;

                    if (flag(ii,jj,kk).isSingleValued() && is_eb_inhomog && Amatrix(a_ind+27,irow) != Real(0.0))
                        rhs(irow) += Amatrix(a_ind+27,irow)*phieb(ii,jj,kk,n);
                }
             }
    }

    cholsol_for_eb(Amatrix, rhs);

    Real dphidn = rhs(1)*nrmx + rhs(2)*nrmy + rhs(3)*nrmz;
    return dphidn;
}

}

#endif
