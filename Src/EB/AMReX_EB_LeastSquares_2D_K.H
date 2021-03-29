#ifndef AMREX_MLEB_LEASTSQUARES_K_H_
#define AMREX_MLEB_LEASTSQUARES_K_H_
#include <AMReX_Config.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

// This is the 2D version (i.e. for 6x6 (A^T A) matrix)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void decomp_chol_np6(Array2D<Real,0,5,0,5>& aa)
{
    int neq = 6;

    Real p[6];
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
                 if (sum1 <= 0.)
                 {
                     p[ii] = 0.0;
                     ising = 1;
                 } else {
                     p[ii] = std::sqrt(sum1);
                 }
            } else {
                if (ising == 0)
                   aa(jj,ii) = sum1 / p[ii];
                else
                   aa(jj,ii) = 0.0;
            }
        }
    }

    for (int ii = 0; ii < neq; ii++)
    {
        for (int jj = ii+1; jj < neq; jj++)
        {
           aa(ii,jj) = 0.0;       // Zero upper triangle
           aa(ii,jj) = aa(jj,ii);  // Zero upper triangle
        }

        aa(ii,ii) = p[ii];
    }
}

// This is the 2D version (i.e. for 6x6 (A^T A) matrix)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void  cholsol_np6(Array2D<Real,0,11,0,5>& Amatrix, Array1D<Real,0,5>& b)
{
    int neq = 6;

    Array2D<Real,0,5,0,5> AtA;

    for (int irow = 0; irow < neq; irow++)
        for (int icol = 0; icol < neq; icol++)
            AtA(irow,icol) =  0.0;

    for (int irow = 0; irow < 12; irow++)
    {
         AtA(0,0) += Amatrix(irow,0)*Amatrix(irow,0); // e^T e
         AtA(0,1) += Amatrix(irow,0)*Amatrix(irow,1); // e^T x
         AtA(0,2) += Amatrix(irow,0)*Amatrix(irow,2); // e^T y
         AtA(0,3) += Amatrix(irow,0)*Amatrix(irow,3); // e^T x^2
         AtA(0,4) += Amatrix(irow,0)*Amatrix(irow,4); // e^T x*y
         AtA(0,5) += Amatrix(irow,0)*Amatrix(irow,5); // e^T y^2
         AtA(1,1) += Amatrix(irow,1)*Amatrix(irow,1); // x^T x
         AtA(1,2) += Amatrix(irow,1)*Amatrix(irow,2); // x^T y
         AtA(1,3) += Amatrix(irow,1)*Amatrix(irow,3); // x^T (x^2)
         AtA(1,4) += Amatrix(irow,1)*Amatrix(irow,4); // x^T (xy)
         AtA(1,5) += Amatrix(irow,1)*Amatrix(irow,5); // x^T (y^2)
         AtA(2,2) += Amatrix(irow,2)*Amatrix(irow,2); // y^T y
         AtA(2,3) += Amatrix(irow,2)*Amatrix(irow,3); // y^T (x^2)
         AtA(2,4) += Amatrix(irow,2)*Amatrix(irow,4); // y^T (xy)
         AtA(2,5) += Amatrix(irow,2)*Amatrix(irow,5); // y^T (y^2)
         AtA(3,3) += Amatrix(irow,3)*Amatrix(irow,3); // (x^2)^T (x^2)
         AtA(3,4) += Amatrix(irow,3)*Amatrix(irow,4); // (x^2)^T (x*y)
         AtA(3,5) += Amatrix(irow,3)*Amatrix(irow,5); // (x^2)^T (y^2)
         AtA(4,4) += Amatrix(irow,4)*Amatrix(irow,4); // (x*y)^T (x*y)
         AtA(4,5) += Amatrix(irow,4)*Amatrix(irow,5); // (x*y)^T (y^2)
         AtA(5,5) += Amatrix(irow,5)*Amatrix(irow,5); // (y^2)^T (y^2)
    }

    for (int irow = 0; irow < neq-1; irow++)
        for (int icol = irow+1; icol < neq; icol++)
           AtA(icol,irow) = AtA(irow,icol);

    decomp_chol_np6(AtA);

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
          b(ii) = 0.0;
       }
    }

    if (AtA(neq-1,neq-1) > 0.)
       b(neq-1) = b(neq-1) / AtA(neq-1,neq-1);
    else
       b(neq-1) = 0.0;

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
          b(ii) = 0.0;
       }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void  cholsol_for_eb(Array2D<Real,0,17,0,5>& Amatrix, Array1D<Real,0,5>& b)
{
    int neq = 6;

    Array2D<Real,0,5,0,5> AtA;

    for (int irow = 0; irow < neq; irow++)
        for (int icol = 0; icol < neq; icol++)
            AtA(irow,icol) =  0.0;

    for (int irow = 0; irow < 18; irow++)
    {
         AtA(0,0) += Amatrix(irow,0)*Amatrix(irow,0); // e^T e
         AtA(0,1) += Amatrix(irow,0)*Amatrix(irow,1); // e^T x
         AtA(0,2) += Amatrix(irow,0)*Amatrix(irow,2); // e^T y
         AtA(0,3) += Amatrix(irow,0)*Amatrix(irow,3); // e^T x^2
         AtA(0,4) += Amatrix(irow,0)*Amatrix(irow,4); // e^T x*y
         AtA(0,5) += Amatrix(irow,0)*Amatrix(irow,5); // e^T y^2
         AtA(1,1) += Amatrix(irow,1)*Amatrix(irow,1); // x^T x
         AtA(1,2) += Amatrix(irow,1)*Amatrix(irow,2); // x^T y
         AtA(1,3) += Amatrix(irow,1)*Amatrix(irow,3); // x^T (x^2)
         AtA(1,4) += Amatrix(irow,1)*Amatrix(irow,4); // x^T (xy)
         AtA(1,5) += Amatrix(irow,1)*Amatrix(irow,5); // x^T (y^2)
         AtA(2,2) += Amatrix(irow,2)*Amatrix(irow,2); // y^T y
         AtA(2,3) += Amatrix(irow,2)*Amatrix(irow,3); // y^T (x^2)
         AtA(2,4) += Amatrix(irow,2)*Amatrix(irow,4); // y^T (xy)
         AtA(2,5) += Amatrix(irow,2)*Amatrix(irow,5); // y^T (y^2)
         AtA(3,3) += Amatrix(irow,3)*Amatrix(irow,3); // (x^2)^T (x^2)
         AtA(3,4) += Amatrix(irow,3)*Amatrix(irow,4); // (x^2)^T (x*y)
         AtA(3,5) += Amatrix(irow,3)*Amatrix(irow,5); // (x^2)^T (y^2)
         AtA(4,4) += Amatrix(irow,4)*Amatrix(irow,4); // (x*y)^T (x*y)
         AtA(4,5) += Amatrix(irow,4)*Amatrix(irow,5); // (x*y)^T (y^2)
         AtA(5,5) += Amatrix(irow,5)*Amatrix(irow,5); // (y^2)^T (y^2)
    }

    for (int irow = 0; irow < neq-1; irow++)
        for (int icol = irow+1; icol < neq; icol++)
           AtA(icol,irow) = AtA(irow,icol);

    decomp_chol_np6(AtA);

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
          b(ii) = 0.0;
       }
    }

    if (AtA(neq-1,neq-1) > 0.)
       b(neq-1) = b(neq-1) / AtA(neq-1,neq-1);
    else
       b(neq-1) = 0.0;

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
          b(ii) = 0.0;
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
                                Real& yloc_on_xface,
                                bool is_eb_dirichlet, bool is_eb_inhomog)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,11,0,5> Amatrix;
    Array1D<Real,0, 5    > rhs;

    // Order of column -- first  six are cell centroids: (i,j-1) (i,j) (i,j+1) (i-1,j-1) (i-1,j) (i-1,j+1)
    // Order of column -- second six are   EB centroids: (i,j-1) (i,j) (i,j+1) (i-1,j-1) (i-1,j) (i-1,j+1)

    for (int irow = 0; irow < 12; irow++)
        for (int icol = 0; icol < 6; icol++)
            Amatrix(irow,icol) =  0.0;

    // Columns: [e x y x*x x*y y*y]
    for (int ii = i-1; ii <= i; ii++) // Normal to face
        for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
        {
            if (!flag(ii,jj,k).isCovered())
            {
                int a_ind  = (jj-(j-1)) + 3*(ii-(i-1));

                // This shifts the x-offset by a normalized half-dx so that
                // x=0 is the x-coordinate of the FACE centroid.
                Real x_off = static_cast<Real>(ii-i) + 0.5;
                Real y_off = static_cast<Real>(jj-j);

                Amatrix(a_ind,0) = 1.0;
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,k,0);
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,k,1) - yloc_on_xface;
                Amatrix(a_ind,3) = Amatrix(a_ind,1) * Amatrix(a_ind,1);
                Amatrix(a_ind,4) = Amatrix(a_ind,1) * Amatrix(a_ind,2);
                Amatrix(a_ind,5) = Amatrix(a_ind,2) * Amatrix(a_ind,2);

                if (!flag(ii,jj,k).isRegular())
                {
                    int a_ind_eb  = a_ind + 6;

                    Amatrix(a_ind_eb,0) = 1.0;
                    Amatrix(a_ind_eb,1) = x_off + bcent(ii,jj,k,0);
                    Amatrix(a_ind_eb,2) = y_off + bcent(ii,jj,k,1) - yloc_on_xface;
                    Amatrix(a_ind_eb,3) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,1);
                    Amatrix(a_ind_eb,4) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,2);
                    Amatrix(a_ind_eb,5) = Amatrix(a_ind_eb,2) * Amatrix(a_ind_eb,2);
                }
            }
        }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 6; irow++)
    {
        rhs(irow) = 0.; // Only non-zero when inhomogeneous Dirichlet

        for (int ii = i-1; ii <= i; ii++)  // Normal to face
            for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
            {
                if (!flag(ii,jj,k).isCovered())
                {
                    int a_ind  = (jj-(j-1)) + 3*(ii-(i-1));
                    rhs(irow) += Amatrix(a_ind  ,irow)*  phi(ii,jj,k,n);

                    if (flag(ii,jj,k).isSingleValued() &&
                        is_eb_dirichlet && is_eb_inhomog)
                            rhs(irow) += Amatrix(a_ind+6,irow)*phieb(ii,jj,k,n);
                }
            }
    }

    cholsol_np6(Amatrix, rhs);

    return rhs(1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_y_of_phi_on_centroids(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Real& xloc_on_yface,
                                bool is_eb_dirichlet, bool is_eb_inhomog)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,11,0,5> Amatrix;
    Array1D<Real,0, 5    > rhs;

    // Order of column -- first  six are cell centroids: (i-1,j) (i,j) (i+1,j) (i-1,j-1) (i,j-1) (i+1,j-1)
    // Order of column -- second six are   EB centroids: (i-1,j) (i,j) (i+1,j) (i-1,j-1) (i,j-1) (i+1,j-1)

    for (int irow = 0; irow < 12; irow++)
        for (int icol = 0; icol < 6; icol++)
            Amatrix(irow,icol) =  0.0;

    // Columns: [e x y x*x x*y y*y]
    for (int jj = j-1; jj <= j; jj++) // Normal to face
        for (int ii = i-1; ii <= i+1; ii++)  // Tangential to face
        {
            if (!flag(ii,jj,k).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(jj-(j-1));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j) + 0.5;

                Amatrix(a_ind,0) = 1.0;
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,k,0) - xloc_on_yface;
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,k,1);
                Amatrix(a_ind,3) = Amatrix(a_ind,1) * Amatrix(a_ind,1);
                Amatrix(a_ind,4) = Amatrix(a_ind,1) * Amatrix(a_ind,2);
                Amatrix(a_ind,5) = Amatrix(a_ind,2) * Amatrix(a_ind,2);

                if (!flag(ii,jj,k).isRegular())
                {
                    int a_ind_eb  = a_ind + 6;
                    Amatrix(a_ind_eb,0) = 1.0;
                    Amatrix(a_ind_eb,1) = x_off + bcent(ii,jj,k,0) - xloc_on_yface;
                    Amatrix(a_ind_eb,2) = y_off + bcent(ii,jj,k,1);
                    Amatrix(a_ind_eb,3) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,1);
                    Amatrix(a_ind_eb,4) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,2);
                    Amatrix(a_ind_eb,5) = Amatrix(a_ind_eb,2) * Amatrix(a_ind_eb,2);
                }
            }
        }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 6; irow++)
    {
        rhs(irow) = 0.; // Only non-zero when inhomogeneous Dirichlet

        for (int jj = j-1; jj <= j; jj++) // Normal to face
            for (int ii = i-1; ii <= i+1; ii++) // Tangential to face
                if (!flag(ii,jj,k).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(jj-(j-1));
                    rhs(irow) += Amatrix(a_ind  ,irow)*  phi(ii,jj,k,n);

                    if (flag(ii,jj,k).isSingleValued() &&
                        is_eb_dirichlet && is_eb_inhomog)
                            rhs(irow) += Amatrix(a_ind+6,irow)*phieb(ii,jj,k,n);
                }
    }

    cholsol_np6(Amatrix, rhs);

    return rhs(2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_eb_of_phi_on_centroids(int i,int j,int k,int n,
                                Array4<Real const> const& phi,
                                Array4<Real const> const& phieb,
                                Array4<EBCellFlag const> const& flag,
                                Array4<Real const> const& ccent,
                                Array4<Real const> const& bcent,
                                Real& nrmx, Real& nrmy,
                                bool is_eb_inhomog)
{
    Array2D<Real,0,17,0,5> Amatrix;
    Array1D<Real,0, 5    > rhs;

    // Order of column -- first 9 are cell centroids, next 9 are EB centroids

    for (int irow = 0; irow < 18; irow++)
        for (int icol = 0; icol < 6; icol++)
            Amatrix(irow,icol) =  0.0;

    //  Column 0-2: [e x y]
    for (int ii = i-1; ii <= i+1; ii++)
        for (int jj = j-1; jj <= j+1; jj++)
        {
            if (!flag(ii,jj,k).isCovered())
            {
                int a_ind  = (jj-(j-1)) + 3*(ii-(i-1));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j);

                if (flag(i,j,k).isConnected((ii-i),(jj-j),0))
                {
                   Amatrix(a_ind,0) = 1.0;
                   Amatrix(a_ind,1) = x_off + ccent(ii,jj,k,0) - bcent(i,j,k,0);
                   Amatrix(a_ind,2) = y_off + ccent(ii,jj,k,1) - bcent(i,j,k,1);
                }

                if (flag(i,j,k).isConnected((ii-i),(jj-j),0) && !flag(ii,jj,k).isRegular())
                {
                    Amatrix(a_ind+9,0) = 1.0;
                    Amatrix(a_ind+9,1) = x_off + bcent(ii,jj,k,0) - bcent(i,j,k,0);
                    Amatrix(a_ind+9,2) = y_off + bcent(ii,jj,k,1) - bcent(i,j,k,1);
                }
            }
        }

    // Columns 3 : [x*x  x*y  y*y]

    for (int irow = 0; irow < 18; irow++)
    {
       Amatrix(irow,3) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,5) =  Amatrix(irow,2) * Amatrix(irow,2);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 6; irow++)
    {
        rhs(irow) = 0.;

        for (int ii = i-1; ii <= i+1; ii++)
            for (int jj = j-1; jj <= j+1; jj++)
                if (!flag(ii,jj,k).isCovered())
                {
                    int a_ind  = (jj-(j-1)) + 3*(ii-(i-1));
                    rhs(irow) += Amatrix(a_ind,irow) * phi(ii,jj,k,n);
                    if (flag(ii,jj,k).isSingleValued() && is_eb_inhomog)
                        rhs(irow) += Amatrix(a_ind+9,irow)*phieb(ii,jj,k,n);
                }
    }

    cholsol_for_eb(Amatrix, rhs);

    Real dphidn = rhs(1)*nrmx + rhs(2)*nrmy;
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
                                       Real& yloc_on_xface,
                                       bool is_eb_dirichlet, bool is_eb_inhomog,
                                       const bool on_x_face, const int domlo_x, const int domhi_x,
                                       const bool on_y_face, const int domlo_y, const int domhi_y)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,17,0,5> Amatrix;
    Array1D<Real,0, 5    > rhs;

    // Order of column -- first  six are cell centroids: (i,j-1) (i,j) (i,j+1) (i-1,j-1) (i-1,j) (i-1,j+1)
    // Order of column -- second six are   EB centroids: (i,j-1) (i,j) (i,j+1) (i-1,j-1) (i-1,j) (i-1,j+1)

    for (int irow = 0; irow < 18; irow++)
        for (int icol = 0; icol < 6; icol++)
            Amatrix(irow,icol) =  0.0;

    const int im = (i > domhi_x) ? 2 : 1;
    const int ip = 2 - im;

    // Columns: [e x y x*x x*y y*y]
    for (int ii = i-im; ii <= i+ip; ii++) // Normal to face
        for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
        {

            // Don't include corner cells. Could make this even more strict
            // by removing the on_?_face restrictions.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii < domlo_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj < domlo_y))) {
              continue;
            }

            if ( !phi.contains(ii,jj,k) )
               continue;

            if (!flag(ii,jj,k).isCovered())
            {

                int a_ind  = (jj-(j-1)) + 3*(ii-(i-im));
                AMREX_ASSERT(0<= a_ind && a_ind <= 8);

                // This shifts the x-offset by a normalized half-dx so that
                // x=0 is the x-coordinate of the FACE centroid.
                Real x_off = static_cast<Real>(ii-i) + 0.5;
                Real y_off = static_cast<Real>(jj-j);

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,k) != 1.0 || vfrac(ii+2,jj,k) != 1.0) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,k) != 1.0 || vfrac(ii-2,jj,k) != 1.0))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,k) != 1.0 || vfrac(ii,jj+2,k) != 1.0) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,k) != 1.0 || vfrac(ii,jj-2,k) != 1.0) )
                        continue;
                }

                Amatrix(a_ind,0) = 1.0;
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,k,0);
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,k,1) - yloc_on_xface;
                Amatrix(a_ind,3) = Amatrix(a_ind,1) * Amatrix(a_ind,1);
                Amatrix(a_ind,4) = Amatrix(a_ind,1) * Amatrix(a_ind,2);
                Amatrix(a_ind,5) = Amatrix(a_ind,2) * Amatrix(a_ind,2);

                // Add in information about the location of the EB. Exclude
                // EBs that are outside the domain.
                if (flag(ii,jj,k).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x &&
                    domlo_y <= jj && jj <= domhi_y)
                {
                    int a_ind_eb  = a_ind + 9;
                    AMREX_ASSERT(9<= a_ind_eb && a_ind_eb <= 17);

                    Amatrix(a_ind_eb,0) = 1.0;
                    Amatrix(a_ind_eb,1) = x_off + bcent(ii,jj,k,0);
                    Amatrix(a_ind_eb,2) = y_off + bcent(ii,jj,k,1) - yloc_on_xface;
                    Amatrix(a_ind_eb,3) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,1);
                    Amatrix(a_ind_eb,4) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,2);
                    Amatrix(a_ind_eb,5) = Amatrix(a_ind_eb,2) * Amatrix(a_ind_eb,2);
                }
            }
        }

    for (int irow = 0; irow < 6; irow++)
    {
        rhs(irow) = 0.; // Only non-zero when inhomogeneous Dirichlet

        for (int ii = i-im; ii <= i+ip; ii++)  // Normal to face
            for (int jj = j-1; jj <= j+1; jj++)  // Tangential to face
            {
                if ( !phi.contains(ii,jj,k) )
                   continue;

                if (!flag(ii,jj,k).isCovered())
                {
                    int a_ind  = (jj-(j-1)) + 3*(ii-(i-im));
                    AMREX_ASSERT(0<=a_ind && a_ind <9);

                    Real phi_val = Amatrix(a_ind,0) * phi(ii,jj,k,n);

                    rhs(irow) += Amatrix(a_ind  ,irow)*phi_val;

                    if (flag(ii,jj,k).isSingleValued() &&
                        is_eb_dirichlet && is_eb_inhomog && Amatrix(a_ind+9,irow) != 0.0)
                            rhs(irow) += Amatrix(a_ind+9,irow)*phieb(ii,jj,k,n);
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
                                       Real& xloc_on_yface,
                                       bool is_eb_dirichlet, bool is_eb_inhomog,
                                       const bool on_x_face, const int domlo_x, const int domhi_x,
                                       const bool on_y_face, const int domlo_y, const int domhi_y)
{
    AMREX_ALWAYS_ASSERT(is_eb_dirichlet);
    Array2D<Real,0,17,0,5> Amatrix;
    Array1D<Real,0, 5    > rhs;

    // Order of column -- first  six are cell centroids: (i-1,j) (i,j) (i+1,j) (i-1,j-1) (i,j-1) (i+1,j-1)
    // Order of column -- second six are   EB centroids: (i-1,j) (i,j) (i+1,j) (i-1,j-1) (i,j-1) (i+1,j-1)

    for (int irow = 0; irow < 18; irow++)
        for (int icol = 0; icol < 6; icol++)
            Amatrix(irow,icol) =  0.0;

    const int jm = (j > domhi_y) ? 2 : 1;
    const int jp = 2 - jm;

    // Columns: [e x y x*x x*y y*y]
    for (int jj = j-jm; jj <= j+jp; jj++) // Normal to face
        for (int ii = i-1; ii <= i+1; ii++)  // Tangential to face
        {

            // Don't include corner cells.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii < domlo_x) && (on_y_face && jj > domhi_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj < domlo_y))) {
              continue;
            }

            if ( !phi.contains(ii,jj,k) )
               continue;

            if (!flag(ii,jj,k).isCovered())
            {
                int a_ind  = (ii-(i-1)) + 3*(jj-(j-jm));
                AMREX_ASSERT(0<= a_ind && a_ind <= 8);

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j) + 0.5;

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,k) != 1.0 || vfrac(ii+2,jj,k) != 1.0) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,k) != 1.0 || vfrac(ii-2,jj,k) != 1.0))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,k) != 1.0 || vfrac(ii,jj+2,k) != 1.0) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,k) != 1.0 || vfrac(ii,jj-2,k) != 1.0) )
                        continue;
                }

                Amatrix(a_ind,0) = 1.0;
                Amatrix(a_ind,1) = x_off + ccent(ii,jj,k,0) - xloc_on_yface;
                Amatrix(a_ind,2) = y_off + ccent(ii,jj,k,1);
                Amatrix(a_ind,3) = Amatrix(a_ind,1) * Amatrix(a_ind,1);
                Amatrix(a_ind,4) = Amatrix(a_ind,1) * Amatrix(a_ind,2);
                Amatrix(a_ind,5) = Amatrix(a_ind,2) * Amatrix(a_ind,2);

                if (flag(ii,jj,k).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x &&
                    domlo_y <= jj && jj <= domhi_y)
                {

                    int a_ind_eb  = a_ind + 9;
                    AMREX_ASSERT(9<= a_ind_eb && a_ind_eb <= 17);

                    Amatrix(a_ind_eb,0) = 1.0;
                    Amatrix(a_ind_eb,1) = x_off + bcent(ii,jj,k,0) - xloc_on_yface;
                    Amatrix(a_ind_eb,2) = y_off + bcent(ii,jj,k,1);
                    Amatrix(a_ind_eb,3) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,1);
                    Amatrix(a_ind_eb,4) = Amatrix(a_ind_eb,1) * Amatrix(a_ind_eb,2);
                    Amatrix(a_ind_eb,5) = Amatrix(a_ind_eb,2) * Amatrix(a_ind_eb,2);
                }
            }
        }


    // Make the RHS = A^T v
    for (int irow = 0; irow < 6; irow++)
    {
        rhs(irow) = 0.; // Only non-zero when inhomogeneous Dirichlet

        for (int jj = j-jm; jj <= j+jp; jj++) // Normal to face
            for (int ii = i-1; ii <= i+1; ii++) {// Tangential to face
                if ( !phi.contains(ii,jj,k) )
                   continue;
                if (!flag(ii,jj,k).isCovered())
                {
                    int a_ind  = (ii-(i-1)) + 3*(jj-(j-jm));
                    AMREX_ASSERT(0<=irow && irow < 9);

                    Real phi_val = Amatrix(a_ind,0)*phi(ii,jj,k,n);

                    rhs(irow) += Amatrix(a_ind,irow)*phi_val;

                    if (flag(ii,jj,k).isSingleValued() &&
                        is_eb_dirichlet && is_eb_inhomog && Amatrix(a_ind+9,irow)){
                      rhs(irow) += Amatrix(a_ind+9,irow)*phieb(ii,jj,k,n);
                    }
                }
            }

    }

    cholsol_for_eb(Amatrix, rhs);

    return rhs(2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real grad_eb_of_phi_on_centroids_extdir(int i,int j,int k,int n,
                                        Array4<Real const> const& phi,
                                        Array4<Real const> const& phieb,
                                        Array4<EBCellFlag const> const& flag,
                                        Array4<Real const> const& ccent,
                                        Array4<Real const> const& bcent,
                                        Array4<Real const> const& vfrac,
                                        Real& nrmx, Real& nrmy,
                                        bool is_eb_inhomog,
                                        const bool on_x_face, const int domlo_x, const int domhi_x,
                                        const bool on_y_face, const int domlo_y, const int domhi_y)
{
    Array2D<Real,0,17,0,5> Amatrix;
    Array1D<Real,0, 5    > rhs;

    // Order of column -- first 9 are cell centroids, next 9 are EB centroids

    for (int irow = 0; irow < 18; irow++)
        for (int icol = 0; icol < 6; icol++)
            Amatrix(irow,icol) =  0.0;

    //  Column 0-2: [e x y]
    for (int ii = i-1; ii <= i+1; ii++)
        for (int jj = j-1; jj <= j+1; jj++)
        {

            // This is likely overkill for EB grads.
            // Don't include corner cells. Could make this even more strict
            // by removing the on_?_face restrictions.
            if (((on_x_face && ii < domlo_x) && (on_y_face && jj < domlo_y)) ||
                ((on_x_face && ii > domhi_x) && (on_y_face && jj > domhi_y))){
              continue;
            }


            if (!flag(ii,jj,k).isCovered())
            {
                int a_ind  = (jj-(j-1)) + 3*(ii-(i-1));

                Real x_off = static_cast<Real>(ii-i);
                Real y_off = static_cast<Real>(jj-j);

                if(on_x_face){
                    if (ii < domlo_x && (vfrac(ii+1,jj,k) != 1.0 || vfrac(ii+2,jj,k) != 1.0) )
                        continue;
                    if (ii > domhi_x && (vfrac(ii-1,jj,k) != 1.0 || vfrac(ii-2,jj,k) != 1.0))
                        continue;
                }

                if(on_y_face){
                    if (jj < domlo_y && (vfrac(ii,jj+1,k) != 1.0 || vfrac(ii,jj+2,k) != 1.0) )
                        continue;
                    if (jj > domhi_y && (vfrac(ii,jj-1,k) != 1.0 || vfrac(ii,jj-2,k) != 1.0) )
                        continue;
                }


                if (flag(i,j,k).isConnected((ii-i),(jj-j),0))
                {
                   Amatrix(a_ind,0) = 1.0;
                   Amatrix(a_ind,1) = x_off + ccent(ii,jj,k,0) - bcent(i,j,k,0);
                   Amatrix(a_ind,2) = y_off + ccent(ii,jj,k,1) - bcent(i,j,k,1);
                }

                // Include information from the EB unless it is outside the domain
                if (flag(i,j,k).isConnected((ii-i),(jj-j),0) && flag(ii,jj,k).isSingleValued() &&
                    domlo_x <= ii && ii <= domhi_x && domlo_y <= jj && jj <= domhi_y)
                {
                    Amatrix(a_ind+9,0) = 1.0;
                    Amatrix(a_ind+9,1) = x_off + bcent(ii,jj,k,0) - bcent(i,j,k,0);
                    Amatrix(a_ind+9,2) = y_off + bcent(ii,jj,k,1) - bcent(i,j,k,1);
                }
            }
        }

    // Columns 3 : [x*x  x*y  y*y]

    for (int irow = 0; irow < 18; irow++)
    {
       Amatrix(irow,3) =  Amatrix(irow,1) * Amatrix(irow,1);
       Amatrix(irow,4) =  Amatrix(irow,1) * Amatrix(irow,2);
       Amatrix(irow,5) =  Amatrix(irow,2) * Amatrix(irow,2);
    }

    // Make the RHS = A^T v
    for (int irow = 0; irow < 6; irow++)
    {
        rhs(irow) = 0.;

        for (int ii = i-1; ii <= i+1; ii++)
            for (int jj = j-1; jj <= j+1; jj++) {
                if ( !phi.contains(ii,jj,k) )
                   continue;
                if (!flag(ii,jj,k).isCovered())
                {
                    int a_ind  = (jj-(j-1)) + 3*(ii-(i-1));
                    AMREX_ASSERT(0<=a_ind && a_ind <9);

                    // A(a_ind,0) is 1 or 0. Use it to mask phi.
                    Real phi_val = Amatrix(a_ind,0)*phi(ii,jj,k,n);

                    rhs(irow) += Amatrix(a_ind,irow) * phi_val;

                    if (flag(ii,jj,k).isSingleValued() && is_eb_inhomog && Amatrix(a_ind+9,irow) != 0.0)
                        rhs(irow) += Amatrix(a_ind+9,irow)*phieb(ii,jj,k,n);
                }
            }
    }

    cholsol_for_eb(Amatrix, rhs);

    Real dphidn = rhs(1)*nrmx + rhs(2)*nrmy;
    return dphidn;
}

}
#endif
