#include <AMReX_EB2_C.H>

namespace amrex::EB2 {

void intercept_to_edge_centroid (AMREX_D_DECL(Array4<Real> const& excent,
                                              Array4<Real> const& eycent,
                                              Array4<Real> const& ezcent),
                                 AMREX_D_DECL(Array4<Type_t const> const& fx,
                                              Array4<Type_t const> const& fy,
                                              Array4<Type_t const> const& fz),
                                 Array4<Real const> const& levset,
                                 GpuArray<Real,AMREX_SPACEDIM> const& dx,
                                 GpuArray<Real,AMREX_SPACEDIM> const& problo) noexcept
{
    AMREX_D_TERM(const Real dxinv = Real(1.0)/dx[0];,
                 const Real dyinv = Real(1.0)/dx[1];,
                 const Real dzinv = Real(1.0)/dx[2];)
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM (
        Box(excent), xbx, {
            AMREX_LOOP_3D(xbx, i, j, k, {
                if (fx(i,j,k) == Type::regular) {
                    excent(i,j,k) = Real(1.0);
                } else if (fx(i,j,k) == Type::covered) {
                    excent(i,j,k) = Real(-1.0);
                } else {
                    Real xcut = Real(0.5)*(excent(i,j,k) - (problo[0]+i*dx[0]))*dxinv;
                    if (levset(i,j,k) < levset(i+1,j,k)) { // right side covered
                        xcut -= Real(0.5);
                    }
                    excent(i,j,k) = amrex::min(Real(0.5),amrex::max(Real(-0.5),xcut));
                }
            });
        },
        Box(eycent), ybx, {
            AMREX_LOOP_3D(ybx, i, j, k, {
                if (fy(i,j,k) == Type::regular) {
                    eycent(i,j,k) = Real(1.0);
                } else if (fy(i,j,k) == Type::covered) {
                    eycent(i,j,k) = Real(-1.0);
                } else {
                    Real ycut = Real(0.5)*(eycent(i,j,k) - (problo[1]+j*dx[1]))*dyinv;
                    if (levset(i,j,k) < levset(i,j+1,k)) { // right side covered
                        ycut -= Real(0.5);
                    }
                    eycent(i,j,k) = amrex::min(Real(0.5),amrex::max(Real(-0.5),ycut));
                }
            });
        },
        Box(ezcent), zbx, {
            AMREX_LOOP_3D(zbx, i, j, k, {
                if (fz(i,j,k) == Type::regular) {
                    ezcent(i,j,k) = Real(1.0);
                } else if (fz(i,j,k) == Type::covered) {
                    ezcent(i,j,k) = Real(-1.0);
                } else {
                    Real zcut = Real(0.5)*(ezcent(i,j,k) - (problo[2]+k*dx[2]))*dzinv;
                    if (levset(i,j,k) < levset(i,j,k+1)) { // right side covered
                        zcut -= Real(0.5);
                    }
                    ezcent(i,j,k) = amrex::min(Real(0.5),amrex::max(Real(-0.5),zcut));
                }
            });
        }
    );
}

}
