#include <cmath>
using namespace amrex;

void get_face_velocity(const int* level, const amrex::real* time,
                       AMREX_D_DECL(BL_FORT_FAB_ARG(xvel),
                                    BL_FORT_FAB_ARG(yvel),
                                    BL_FORT_FAB_ARG(zvel)),
                       const amrex::real* dx, const amrex::real* problo)
{
    int plo[2], phi[2];
    Real x, y, sx2, sy2;
    Real** psi;
    const Real pi = 3.141592653589793238462643383279502884197;

    plo[0] = min( vx_l1-1, vy_l1-1 );
    plo[1] = min( vx_l2-1, vy_l2-1 );
    phi[0] = max( vx_h1  , vy_h1+1 );
    phi[1] = max( vx_h2+1, vy_h2   );

    // call bl_allocate( psi, plo(1), phi(1), plo(2), phi(2) )

    // streamfunction psi
    for (int i = plo[0]; i <= phi[0]; i++)
    {
        x = ( (Real)i + 0.5 )*dx[0] + prob_lo[0];
        sx2 = sin(pi*x) * sin(pi*x);
        for (int j = plo[1]; j <= phi[1]; j++)
        {
            y = ( (Real)j + 0.5 )*dx[1] + prob_lo[1];
            sy2 = sin(pi*y) * sin(pi*y);
            psi[i][j] = sx2 * sy2 * cos(pi*time/2.0) * (1.0/pi);
        }
    }

    // x velocity
    for (int i = vx_l1; i <= vx_h1; i++)
    {
        x = (Real)i * dx[0] + prob_lo[0];
        for (int j = vx_l2; j <= vx_h2; j++)
        {
            y = ( (Real)j + 0.5 )*dx[1] + prob_lo[1];
            vx[i][j] = -(  ( psi(i,j+1) + psi(i-1,j+1) )
                         - ( psi(i,j-1) + psi(i-1,j-1) ) ) * (0.25/dx[1]);
        }
    }

    // y velocity
    for (int i = vy_l1; i <= vy_h1; i++)
    {
        x = ( (Real)i + 0.5 )*dx[0] + prob_lo[0];
        for (int j = vx_l2; j <= vx_h2; j++)
        {
            y = (Real)i * dx[1] + prob_lo[1];
            vy[i][j] =  (  ( psi(i+1,j) + psi(i+1,j-1) )
                         - ( psi(i-1,j) + psi(i-1,j-1) ) ) * (0.25/dx[0]);
        }
    }

    // call bl_deallocate(psi)

    return;
}
