#include <AMReX_EB2_C.H>

namespace amrex { namespace EB2 {

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_eb_data (const int i, const int j,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm)
{
    constexpr Real almostone = 1.0-1.e-15;
    constexpr Real small = 1.e-14;
    constexpr Real tiny = 1.e-15;

    const Real axm = apx(i  ,j  ,0);
    const Real axp = apx(i+1,j  ,0);
    const Real aym = apy(i  ,j  ,0);
    const Real ayp = apy(i  ,j+1,0);
    const Real apnorm = std::hypot(axm-axp,aym-ayp);
    const Real nx = (axm-axp) * (1.0/apnorm);
    const Real ny = (aym-ayp) * (1.0/apnorm);

    const Real nxabs = amrex::Math::abs(nx);
    const Real nyabs = amrex::Math::abs(ny);

    Real x_ym;
    Real x_yp;
    Real signx;
    Real y_xm;
    Real y_xp;
    Real signy;
    if (nx > 0.0) {
        x_ym = -0.5 + aym;
        x_yp = -0.5 + ayp;
        signx = 1.0;
    } else {
        x_ym = 0.5 - aym;
        x_yp = 0.5 - ayp;
        signx = -1.0;
    }

    if (ny > 0.0) {
        y_xm = -0.5 + axm;
        y_xp = -0.5 + axp;
        signy = 1.0;
    } else {
        y_xm = 0.5 - axm;
        y_xp = 0.5 - axp;
        signy = -1.0;
    }

    barea(i,j,0) = nx*(axm-axp) + ny*(aym-ayp);
    bcent(i,j,0,0) = 0.5*(x_ym+x_yp);
    bcent(i,j,0,1) = 0.5*(y_xm+y_xp);
    bnorm(i,j,0,0) = nx;
    bnorm(i,j,0,1) = ny;

    if (nxabs < tiny || nyabs > almostone) {
        barea(i,j,0) = 1.0;
        bcent(i,j,0,0) = 0.0;
        bnorm(i,j,0,0) = 0.0;
        bnorm(i,j,0,1) = signy;
        vfrac(i,j,0) = 0.5*(axm+axp);
        vcent(i,j,0,0) = 0.0;
        vcent(i,j,0,1) = (0.125*(ayp-aym) + ny*0.5*bcent(i,j,0,1)*bcent(i,j,0,1)) / (vfrac(i,j,0) + 1.e-30);
    } else if (nyabs < tiny || nxabs > almostone) {
        barea(i,j,0) = 1.0;
        bcent(i,j,0,1) = 0.0;
        bnorm(i,j,0,0) = signx;
        bnorm(i,j,0,1) = 0.0;
        vfrac(i,j,0) = 0.5*(aym+ayp);
        vcent(i,j,0,0) = (0.125*(axp-axm) + nx*0.5*bcent(i,j,0,0)*bcent(i,j,0,0)) / (vfrac(i,j,0) + 1.e-30);
        vcent(i,j,0,1) = 0.0;
    } else {
        Real aa = nxabs/ny;
        const Real dx = x_ym - x_yp;
        const Real dx2 = dx * (x_ym + x_yp);
        const Real dx3 = dx * (x_ym*x_ym + x_ym*x_yp + x_yp*x_yp);
        const Real af1 = 0.5*(axm+axp) + aa*0.5*dx2;
        vcent(i,j,0,0) = 0.125*(axp-axm) + aa*(1./6.)*dx3;

        aa = nyabs/nx;
        const Real dy = y_xm - y_xp;
        const Real dy2 = dy * (y_xm + y_xp);
        const Real dy3 = dy * (y_xm*y_xm + y_xm*y_xp + y_xp*y_xp);
        const Real af2 = 0.5*(aym+ayp) + aa*0.5*dy2;
        vcent(i,j,0,1) = 0.125*(ayp-aym) + aa*(1./6.)*dy3;

        vfrac(i,j,0) = 0.5*(af1+af2);
        if (vfrac(i,j,0) > 1.0-small) {
            vfrac(i,j,0) = 1.0;
            vcent(i,j,0,0) = 0.0;
            vcent(i,j,0,1) = 0.0;
        } else if (vfrac(i,j,0) < small) {
            vfrac(i,j,0) = 0.0;
            vcent(i,j,0,0) = 0.0;
            vcent(i,j,0,1) = 0.0;
        } else {
            vcent(i,j,0,0) *= (1.0/vfrac(i,j,0));
            vcent(i,j,0,1) *= (1.0/vfrac(i,j,0));
            vcent(i,j,0,0) = amrex::min(amrex::max(vcent(i,j,0,0),-0.5),0.5);
            vcent(i,j,0,1) = amrex::min(amrex::max(vcent(i,j,0,1),-0.5),0.5);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_covered(const int i, const int j,
                 Array4<EBCellFlag> const& cell,
                 Array4<Real> const& vfrac, Array4<Real> const& vcent,
                 Array4<Real> const& barea, Array4<Real> const& bcent,
                 Array4<Real> const& bnorm) 
{
   vfrac(i,j,0) = 0.0;
   vcent(i,j,0,0) = 0.0;
   vcent(i,j,0,1) = 0.0;
   barea(i,j,0) = 0.0;
   bcent(i,j,0,0) = -1.0;
   bcent(i,j,0,1) = -1.0;
   bnorm(i,j,0,0) = 0.0;
   bnorm(i,j,0,1) = 0.0;
   cell(i,j,0).setCovered();
}

}

void build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Real const> const& levset,
                  Array4<Real const> const& interx, Array4<Real const> const& intery,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& fcx, Array4<Real> const& fcy,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo)
{
    constexpr Real small = 1.e-14;
    const Real dxinv = 1.0/dx[0];
    const Real dyinv = 1.0/dx[1];
    const Box& ndbxg1 = amrex::grow(amrex::surroundingNodes(bx),1);
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( ndbxg1, tbx,
    {
        Box lbx = amrex::grow(amrex::surroundingNodes(bx,0),1);
        auto lo = amrex::max_lbound(tbx, lbx);
        auto hi = amrex::min_ubound(tbx, lbx);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i)
        {
            if (fx(i,j,0) == Type::regular) {
                apx(i,j,0) = 1.0;
                fcx(i,j,0) = 0.0;
            } else if (fx(i,j,0) == Type::covered) {
                apx(i,j,0) = 0.0;
                fcx(i,j,0) = 0.0;
            } else {
                if (levset(i,j,0) < 0.0) {
                    apx(i,j,0) = (interx(i,j,0)-(problo[1]+j*dx[1]))*dyinv;
                    fcx(i,j,0) = 0.5*apx(i,j,0) - 0.5;
                } else {
                    apx(i,j,0) = 1.0 - (interx(i,j,0)-(problo[1]+j*dx[1]))*dyinv;
                    fcx(i,j,0) = 0.5 - 0.5*apx(i,j,0);
                }

                if (apx(i,j,0) > 1.0-small) {
                    apx(i,j,0) = 1.0;
                    fcx(i,j,0) = 0.0;
                    fx(i,j,0) = Type::regular;
                } else if (apx(i,j,0) < small) {
                    apx(i,j,0) = 0.0;
                    fcx(i,j,0) = 0.0;
                    fx(i,j,0) = Type::covered;
                }
            }
        }}

        lbx = amrex::grow(amrex::surroundingNodes(bx,1),1);
        lo = amrex::max_lbound(tbx, lbx);
        hi = amrex::min_ubound(tbx, lbx);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i)
        {
            if (fy(i,j,0) == Type::regular) {
                apy(i,j,0) = 1.0;
                fcy(i,j,0) = 0.0;
            } else if (fy(i,j,0) == Type::covered) {
                apy(i,j,0) = 0.0;
                fcy(i,j,0) = 0.0;
            } else {
                if (levset(i,j,0) < 0.0) {
                    apy(i,j,0) = (intery(i,j,0)-(problo[0]+i*dx[0]))*dxinv;
                    fcy(i,j,0) = 0.5*apy(i,j,0) - 0.5;
                } else {
                    apy(i,j,0) = 1.0 - (intery(i,j,0)-(problo[0]+i*dx[0]))*dxinv;
                    fcy(i,j,0) = 0.5 - 0.5*apy(i,j,0);
                }

                if (apy(i,j,0) > 1.0-small) {
                    apy(i,j,0) = 1.0;
                    fcy(i,j,0) = 0.0;
                    fy(i,j,0) = Type::regular;
                } else if (apy(i,j,0) < small) {
                    apy(i,j,0) = 0.0;
                    fcy(i,j,0) = 0.0;
                    fy(i,j,0) = Type::covered;
                }
            }
        }}
    });

    const Box& bxg1 = amrex::grow(bx, 1);
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        if (cell(i,j,0).isSingleValued()) {
            if (fx(i,j,0) == Type::regular && fx(i+1,j,0) == Type::regular &&
                fy(i,j,0) == Type::regular && fy(i,j+1,0) == Type::regular)
            {
                cell(i,j,0).setRegular();
            }
            else if (fx(i,j,0) == Type::covered && fx(i+1,j,0) == Type::covered &&
                     fy(i,j,0) == Type::covered && fy(i,j+1,0) == Type::covered)
            {
                cell(i,j,0).setCovered();
            }
            else
            {
                int ncuts = 0;
                if (fx(i  ,j  ,0) == Type::irregular) ++ncuts;
                if (fx(i+1,j  ,0) == Type::irregular) ++ncuts;
                if (fy(i  ,j  ,0) == Type::irregular) ++ncuts;
                if (fy(i  ,j+1,0) == Type::irregular) ++ncuts;
                if (ncuts > 2) {
                    amrex::Error("amerx::EB2::build_faces: more than 2 cuts not supported");
                }
            }
        }
    });
}

void build_cells (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Real small_volfrac,
                  Geometry const& geom, bool extend_domain_face)
{
    const Box& bxg1 = amrex::grow(bx,1);
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        if (cell(i,j,k).isRegular()) {
            vfrac(i,j,0) = 1.0;
            vcent(i,j,0,0) = 0.0;
            vcent(i,j,0,1) = 0.0;
            barea(i,j,0) = 0.0;
            bcent(i,j,0,0) = -1.0;
            bcent(i,j,0,1) = -1.0;
            bnorm(i,j,0,0) = 0.0;
            bnorm(i,j,0,1) = 0.0;
        } else if (cell(i,j,k).isCovered()) {
            vfrac(i,j,0) = 0.0;
            vcent(i,j,0,0) = 0.0;
            vcent(i,j,0,1) = 0.0;
            barea(i,j,0) = 0.0;
            bcent(i,j,0,0) = -1.0;
            bcent(i,j,0,1) = -1.0;
            bnorm(i,j,0,0) = 0.0;
            bnorm(i,j,0,1) = 0.0;
        } else {

            set_eb_data(i,j,apx,apy,vfrac,vcent,barea,bcent,bnorm);

            // remove small cells
            if (vfrac(i,j,0) < small_volfrac) {
               set_covered(i,j,cell,vfrac,vcent,barea,bcent,bnorm);
            }
        }
    });

    // set cells in the extended region to covered if the
    // corresponding cell on the domain face is covered
    if(extend_domain_face) {

       Box gdomain = geom.Domain();
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
           if (geom.isPeriodic(idim)) {
               gdomain.setSmall(idim, std::min(gdomain.smallEnd(idim), bxg1.smallEnd(idim)));
               gdomain.setBig(idim, std::max(gdomain.bigEnd(idim), bxg1.bigEnd(idim)));
           }
       }

       if (! gdomain.contains(bxg1)) {
       AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
       {
              const auto & dlo = gdomain.loVect();
              const auto & dhi = gdomain.hiVect();

              // find the cell(ii,jj,kk) on the corr. domain face
              // this would have already been set to correct value
              bool in_extended_domain = false;
              int ii = i;
              int jj = j;
              int kk = k;
              if(i < dlo[0]) {
                  in_extended_domain = true;
                  ii = dlo[0];
              }
              else if(i > dhi[0]) {
                  in_extended_domain = true;
                  ii = dhi[0];
              }

              if(j < dlo[1]) {
                  in_extended_domain = true;
                  jj = dlo[1];
              }
              else if(j > dhi[1]) {
                  in_extended_domain = true;
                  jj = dhi[1];
              }

              // set cell in extendable region to covered if necessary
              if( in_extended_domain && (! cell(i,j,k).isCovered()) 
                  && cell(ii,jj,kk).isCovered() ) 
              {
                  set_covered(i,j,cell,vfrac,vcent,barea,bcent,bnorm);
              }
          });
       }
    }


    // fix face for small cells whose vfrac has been set to zero
    const Box xbx = Box(bx).surroundingNodes(0).grow(1,1);
    AMREX_HOST_DEVICE_FOR_3D ( xbx, i, j, k,
    {
        if (vfrac(i-1,j,0) == 0._rt || vfrac(i,j,0) == 0._rt) {
            fx(i,j,0) = Type::covered;
            apx(i,j,0) = 0.0;
            // race conditions do not happeen because multiple cuts are not allowed
            if (! cell(i,j,0).isCovered())
            {
                cell(i,j,0).setSingleValued();
                set_eb_data(i,j,apx,apy,vfrac,vcent,barea,bcent,bnorm);
            }
            if (! cell(i-1,j,0).isCovered())
            {
                cell(i-1,j,0).setSingleValued();
                set_eb_data(i-1,j,apx,apy,vfrac,vcent,barea,bcent,bnorm);
            }
        }
    });
    //
    const Box ybx = Box(bx).surroundingNodes(1).grow(0,1);
    AMREX_HOST_DEVICE_FOR_3D ( ybx, i, j, k,
    {
        if (vfrac(i,j-1,0) == 0._rt || vfrac(i,j,0) == 0._rt) {
            fy(i,j,0) = Type::covered;
            apy(i,j,0) = 0.0;
            if (! cell(i,j,0).isCovered())
            {
                cell(i,j,0).setSingleValued();
                set_eb_data(i,j,apx,apy,vfrac,vcent,barea,bcent,bnorm);
            }
            if (! cell(i,j-1,0).isCovered())
            {
                cell(i,j-1,0).setSingleValued();
                set_eb_data(i,j-1,apx,apy,vfrac,vcent,barea,bcent,bnorm);
            }
        }
    });

    // Build neighbors.  By default, all neighbors are already set.
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        auto flg = cell(i,j,0);

        if (fx(i  ,j  ,0) == Type::covered) flg.setDisconnected(IntVect(-1, 0));
        if (fx(i+1,j  ,0) == Type::covered) flg.setDisconnected(IntVect( 1, 0));
        if (fy(i  ,j  ,0) == Type::covered) flg.setDisconnected(IntVect( 0,-1));
        if (fy(i  ,j+1,0) == Type::covered) flg.setDisconnected(IntVect( 0, 1));

        if (((fx(i,j,0) == Type::covered) || fy(i-1,j,0) == Type::covered) &&
            ((fx(i,j-1,0) == Type::covered) || fy(i,j,0) == Type::covered))
        {
            flg.setDisconnected(IntVect(-1,-1));
        }

        if (((fx(i+1,j,0) == Type::covered) || fy(i+1,j,0) == Type::covered) &&
            ((fx(i+1,j-1,0) == Type::covered) || fy(i,j,0) == Type::covered))
        {
            flg.setDisconnected(IntVect(1,-1));
        }

        if (((fx(i,j,0) == Type::covered) || fy(i-1,j+1,0) == Type::covered) &&
            ((fx(i,j+1,0) == Type::covered) || fy(i,j+1,0) == Type::covered))
        {
            flg.setDisconnected(IntVect(-1,1));
        }

        if (((fx(i+1,j,0) == Type::covered) || fy(i+1,j+1,0) == Type::covered) &&
            ((fx(i+1,j+1,0) == Type::covered) || fy(i,j+1,0) == Type::covered))
        {
            flg.setDisconnected(IntVect(1,1));
        }

        cell(i,j,0) = flg;
    });
}

}}
