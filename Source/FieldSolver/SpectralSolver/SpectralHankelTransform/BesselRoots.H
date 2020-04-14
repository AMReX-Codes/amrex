/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
/* -------------------------------------------------------------------------
! program to calculate the first zeroes (root abscissas) of the first
! kind bessel function of integer order n using the subroutine rootj.
! --------------------------------------------------------------------------
! sample run:
!
! (calculate the first 10 zeroes of 1st kind bessel function of order 2).
!
! zeroes of bessel function of order: 2
!
! number of calculated zeroes: 10
!
! table of root abcissas (5 items per line)
!    5.135622    8.417244   11.619841   14.795952   17.959819
    21.116997   24.270112   27.420574   30.569204   33.716520
!
! table of error codes (5 items per line)
!   0   0   0   0   0
!   0   0   0   0   0
!
! --------------------------------------------------------------------------
! reference: from numath library by tuan dang trong in fortran 77
!            [bibli 18].
!
!                               c++ release 1.0 by j-p moreau, paris
!                                          (www.jpmoreau.fr)
! ------------------------------------------------------------------------ */

using amrex::operator""_rt;

void SecantRootFinder(int n, int nitmx, amrex::Real tol, amrex::Real *zeroj, int *ier);

/*----------------------------------------------------------------------
 *    calculate the first nk zeroes of bessel function j(n, x)
 *    including the trivial root (when n > 0)
 *
 *    inputs:
 *      n    order of function j (integer >= 0)                  i*4
 *      nk   number of first zeroes  (integer > 0)               i*4
 *    outputs:
 *      roots(nk)  table of first zeroes (abcissas)              r*8
 *      ier(nk)    table of error codes (must be zeroes)         i*4
 *
 *    reference :
 *    abramowitz m. & stegun irene a.
 *    handbook of mathematical functions
 */
void GetBesselRoots(int n, int nk, HankelTransform::RealVector& roots, amrex::Vector<int> &ier)  {
    amrex::Real zeroj;
    int ierror, ik, k;

    const amrex::Real tol = 1e-14_rt;
    const amrex::Real nitmx = 10;

    const amrex::Real c1 = 1.8557571_rt;
    const amrex::Real c2 = 1.033150_rt;
    const amrex::Real c3 = 0.00397_rt;
    const amrex::Real c4 = 0.0908_rt;
    const amrex::Real c5 = 0.043_rt;

    const amrex::Real t0 = 4.0_rt*n*n;
    const amrex::Real t1 = t0 - 1.0_rt;
    const amrex::Real t3 = 4.0_rt*t1*(7.0_rt*t0 - 31.0_rt);
    const amrex::Real t5 = 32.0_rt*t1*((83.0_rt*t0 - 982.0_rt)*t0 + 3779.0_rt);
    const amrex::Real t7 = 64.0_rt*t1*(((6949.0_rt*t0 - 153855.0_rt)*t0 + 1585743.0_rt)*t0 - 6277237.0_rt);

    roots.resize(nk);
    ier.resize(nk);

    // first zero
    if (n == 0) {
        zeroj = c1 + c2 - c3 - c4 + c5;
        SecantRootFinder(n, nitmx, tol, &zeroj, &ierror);
        ier[0] = ierror;
        roots[0] = zeroj;
        ik = 1;
    }
    else {
        // Include the trivial root
        ier[0] = 0;
        roots[0] = 0.;
        const amrex::Real f1 = std::pow(n, (1.0_rt/3.0_rt));
        const amrex::Real f2 = f1*f1*n;
        const amrex::Real f3 = f1*n*n;
        zeroj = n + c1*f1 + (c2/f1) - (c3/n) - (c4/f2) + (c5/f3);
        SecantRootFinder(n, nitmx, tol, &zeroj, &ierror);
        ier[1] = ierror;
        roots[1] = zeroj;
        ik = 2;
    }

    // other zeroes
    // k counts the nontrivial roots
    // ik counts all roots
    k = 2;
    while (ik < nk) {

        // mac mahon's series for k >> n
        const amrex::Real b0 = (k + 0.5_rt*n - 0.25_rt)*MathConst::pi;
        const amrex::Real b1 = 8.0_rt*b0;
        const amrex::Real b2 = b1*b1;
        const amrex::Real b3 = 3.0_rt*b1*b2;
        const amrex::Real b5 = 5.0_rt*b3*b2;
        const amrex::Real b7 = 7.0_rt*b5*b2;

        zeroj = b0 - (t1/b1) - (t3/b3) - (t5/b5) - (t7/b7);

        const amrex::Real errj = std::abs(jn(n, zeroj));

        // improve solution using procedure SecantRootFinder
        if (errj > tol) SecantRootFinder(n, nitmx, tol, &zeroj, &ierror);

        roots[ik] = zeroj;
        ier[ik] = ierror;

        k++;
        ik++;
    }
}

void SecantRootFinder(int n, int nitmx, amrex::Real tol, amrex::Real *zeroj, int *ier) {

    amrex::Real p0, p1, q0, q1, dp, p;
    amrex::Real c[2];

    c[0] = 0.95_rt;
    c[1] = 0.999_rt;
    *ier = 0;

    for (int ntry=0 ; ntry <= 1 ; ntry++) {
        p0 = c[ntry]*(*zeroj);

        p1 = *zeroj;
        q0 = jn(n, p0);
        q1 = jn(n, p1);
        for (int it=1; it <= nitmx; it++) {
            if (q1 == q0) break;
            p = p1 - q1*(p1 - p0)/(q1 - q0);
            dp = p - p1;
            if (it > 1 && std::abs(dp) < tol) {
                *zeroj = p;
                return;
            }
            p0 = p1;
            q0 = q1;
            p1 = p;
            q1 = jn(n, p1);
        }
    }
    *ier = 3;
    *zeroj = p;
}
