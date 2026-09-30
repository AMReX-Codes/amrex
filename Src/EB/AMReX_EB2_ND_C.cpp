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
    // Classify by level-set sign, not type: 2D face types may have been
    // promoted to regular by the small-area tolerance in build_faces.
    AMREX_D_TERM(const Real dxinv = Real(1.0)/dx[0];,
                 const Real dyinv = Real(1.0)/dx[1];,
                 const Real dzinv = Real(1.0)/dx[2];)
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM (
        Box(excent), xbx, {
            AMREX_LOOP_3D(xbx, i, j, k, {
                bool const cut = (levset(i,j,k) < Real(0.0)) != (levset(i+1,j,k) < Real(0.0));
                if (fx(i,j,k) == Type::regular && !cut) {
                    excent(i,j,k) = Real(1.0);
                } else if (fx(i,j,k) == Type::covered && !cut) {
                    excent(i,j,k) = Real(-1.0);
                } else {
                    Real xcut = Real(0.5)*(excent(i,j,k) - (problo[0]+Real(i)*dx[0]))*dxinv;
                    if (levset(i,j,k) < levset(i+1,j,k)) { // right side covered
                        xcut -= Real(0.5);
                    }
                    excent(i,j,k) = amrex::min(Real(0.5),amrex::max(Real(-0.5),xcut));
                }
            });
        },
        Box(eycent), ybx, {
            AMREX_LOOP_3D(ybx, i, j, k, {
                bool const cut = (levset(i,j,k) < Real(0.0)) != (levset(i,j+1,k) < Real(0.0));
                if (fy(i,j,k) == Type::regular && !cut) {
                    eycent(i,j,k) = Real(1.0);
                } else if (fy(i,j,k) == Type::covered && !cut) {
                    eycent(i,j,k) = Real(-1.0);
                } else {
                    Real ycut = Real(0.5)*(eycent(i,j,k) - (problo[1]+Real(j)*dx[1]))*dyinv;
                    if (levset(i,j,k) < levset(i,j+1,k)) { // right side covered
                        ycut -= Real(0.5);
                    }
                    eycent(i,j,k) = amrex::min(Real(0.5),amrex::max(Real(-0.5),ycut));
                }
            });
        },
        Box(ezcent), zbx, {
            AMREX_LOOP_3D(zbx, i, j, k, {
                bool const cut = (levset(i,j,k) < Real(0.0)) != (levset(i,j,k+1) < Real(0.0));
                if (fz(i,j,k) == Type::regular && !cut) {
                    ezcent(i,j,k) = Real(1.0);
                } else if (fz(i,j,k) == Type::covered && !cut) {
                    ezcent(i,j,k) = Real(-1.0);
                } else {
                    Real zcut = Real(0.5)*(ezcent(i,j,k) - (problo[2]+Real(k)*dx[2]))*dzinv;
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
