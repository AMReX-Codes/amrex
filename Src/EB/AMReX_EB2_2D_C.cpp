#include <AMReX_EB2_C.H>

namespace amrex::EB2 {

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_eb_data (const int i, const int j,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm) noexcept
{
#ifdef AMREX_USE_FLOAT
    constexpr Real almostone = 1.0_rt-1.e-6_rt;
    constexpr Real small = 1.e-5_rt;
    constexpr Real tiny = 1.e-6_rt;
#else
    constexpr Real almostone = 1.0-1.e-15;
    constexpr Real small = 1.e-14;
    constexpr Real tiny = 1.e-15;
#endif

    const Real axm = apx(i  ,j  ,0)*dx[1];
    const Real axp = apx(i+1,j  ,0)*dx[1];
    const Real aym = apy(i  ,j  ,0)*dx[0];
    const Real ayp = apy(i  ,j+1,0)*dx[0];
    const Real daxp = (axm-axp);
    const Real dayp = (aym-ayp);
    const Real apnorm = std::hypot(daxp,dayp) + 1.e-30_rt*std::sqrt(dx[0]*dx[1]);
    const Real nx = daxp * (1.0_rt/apnorm);
    const Real ny = dayp * (1.0_rt/apnorm);
    const Real bareascaling = std::sqrt( (nx*dx[0])*(nx*dx[0]) +
            (ny*dx[1])*(ny*dx[1]) );

    const Real nxabs = std::abs(nx);
    const Real nyabs = std::abs(ny);

    Real x_ym, x_yp, y_xm, y_xp;
    if (nx > 0.0_rt) {
        x_ym = -0.5_rt*dx[0] + aym;
        x_yp = -0.5_rt*dx[0] + ayp;
    } else {
        x_ym = 0.5_rt*dx[0] - aym;
        x_yp = 0.5_rt*dx[0] - ayp;
    }

    if (ny > 0.0_rt) {
        y_xm = -0.5_rt*dx[1] + axm;
        y_xp = -0.5_rt*dx[1] + axp;
    } else {
        y_xm = 0.5_rt*dx[1] - axm;
        y_xp = 0.5_rt*dx[1] - axp;
    }

    barea(i,j,0) = (nx*daxp + ny*dayp)/bareascaling;
    bcent(i,j,0,0) = 0.5_rt*(x_ym+x_yp);
    bcent(i,j,0,1) = 0.5_rt*(y_xm+y_xp);
    bnorm(i,j,0,0) = nx;
    bnorm(i,j,0,1) = ny;

    if (nxabs < tiny || nyabs > almostone) {
        vfrac(i,j,0) = 0.5_rt*(axm+axp)/dx[1];
        vcent(i,j,0,0) = 0.0_rt;
        if (vfrac(i,j,0) > almostone) {
            vcent(i,j,0,1) = 0.0_rt;
        } else {
            vcent(i,j,0,1) = (-0.125_rt*dayp*dx[1]*dx[1] + ny*dx[0]*0.5_rt*bcent(i,j,0,1)*bcent(i,j,0,1)) / ((vfrac(i,j,0) + 1.e-30_rt) * (dx[0]*dx[1]*dx[1]));
        }
    } else if (nyabs < tiny || nxabs > almostone) {
        vfrac(i,j,0) = 0.5_rt*(aym+ayp)/dx[0];
        if (vfrac(i,j,0) > almostone) {
            vcent(i,j,0,0) = 0.0_rt;
        } else {
            vcent(i,j,0,0) = (-0.125_rt*daxp*dx[0]*dx[0] + nx*dx[1]*0.5_rt*bcent(i,j,0,0)*bcent(i,j,0,0)) / ((vfrac(i,j,0) + 1.e-30_rt) * (dx[0]*dx[0]*dx[1]));
        }
        vcent(i,j,0,1) = 0.0_rt;
    } else {
        Real aa = nxabs/ny*dx[0]/dx[1];
        const Real dxx = x_ym - x_yp;
        const Real dx2 = dxx * (x_ym + x_yp);
        const Real dx3 = dxx * (x_ym*x_ym + x_ym*x_yp + x_yp*x_yp);
        const Real af1 = 0.5_rt*(axm+axp)*dx[0] + aa*0.5_rt*dx2;
        vcent(i,j,0,0) = -0.125_rt*daxp*dx[0]*dx[0] + aa*(1._rt/6._rt)*dx3;

        aa = nyabs/nx*dx[1]/dx[0];
        const Real dy = y_xm - y_xp;
        const Real dy2 = dy * (y_xm + y_xp);
        const Real dy3 = dy * (y_xm*y_xm + y_xm*y_xp + y_xp*y_xp);
        const Real af2 = 0.5_rt*(aym+ayp)*dx[1] + aa*0.5_rt*dy2;
        vcent(i,j,0,1) = -0.125_rt*dayp*dx[1]*dx[1] + aa*(1._rt/6._rt)*dy3;

        vfrac(i,j,0) = 0.5_rt*(af1+af2)/(dx[0]*dx[1]);

        if (vfrac(i,j,0) > 1.0_rt-small) {
            vfrac(i,j,0) = 1.0_rt;
            vcent(i,j,0,0) = 0.0_rt;
            vcent(i,j,0,1) = 0.0_rt;
        } else if (vfrac(i,j,0) < small) {
            vfrac(i,j,0) = 0.0_rt;
            vcent(i,j,0,0) = 0.0_rt;
            vcent(i,j,0,1) = 0.0_rt;
        } else {
            vcent(i,j,0,0) *= (1.0_rt/(vfrac(i,j,0)*dx[0]*dx[0]*dx[1]));
            vcent(i,j,0,1) *= (1.0_rt/(vfrac(i,j,0)*dx[0]*dx[1]*dx[1]));
            vcent(i,j,0,0) = amrex::min(amrex::max(vcent(i,j,0,0),Real(-0.5)),Real(0.5));
            vcent(i,j,0,1) = amrex::min(amrex::max(vcent(i,j,0,1),Real(-0.5)),Real(0.5));
        }
    }
    bcent(i,j,0,0) /= dx[0];
    bcent(i,j,0,1) /= dx[1];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_covered(const int i, const int j,
                 Array4<EBCellFlag> const& cell,
                 Array4<Real> const& vfrac, Array4<Real> const& vcent,
                 Array4<Real> const& barea, Array4<Real> const& bcent,
                 Array4<Real> const& bnorm) noexcept
{
   vfrac(i,j,0) = 0.0_rt;
   vcent(i,j,0,0) = 0.0_rt;
   vcent(i,j,0,1) = 0.0_rt;
   barea(i,j,0) = 0.0_rt;
   bcent(i,j,0,0) = -1.0_rt;
   bcent(i,j,0,1) = -1.0_rt;
   bnorm(i,j,0,0) = 0.0_rt;
   bnorm(i,j,0,1) = 0.0_rt;
   cell(i,j,0).setCovered();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool set_eb_cell (int i, int j, Array4<EBCellFlag> const& cell,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Real small_volfrac) noexcept
{
    bool is_small_cell = false;
    if (cell(i,j,0).isRegular()) {
        vfrac(i,j,0) = 1.0_rt;
        vcent(i,j,0,0) = 0.0_rt;
        vcent(i,j,0,1) = 0.0_rt;
        barea(i,j,0) = 0.0_rt;
        bcent(i,j,0,0) = -1.0_rt;
        bcent(i,j,0,1) = -1.0_rt;
        bnorm(i,j,0,0) = 0.0_rt;
        bnorm(i,j,0,1) = 0.0_rt;
    } else if (cell(i,j,0).isCovered()) {
        vfrac(i,j,0) = 0.0_rt;
        vcent(i,j,0,0) = 0.0_rt;
        vcent(i,j,0,1) = 0.0_rt;
        barea(i,j,0) = 0.0_rt;
        bcent(i,j,0,0) = -1.0_rt;
        bcent(i,j,0,1) = -1.0_rt;
        bnorm(i,j,0,0) = 0.0_rt;
        bnorm(i,j,0,1) = 0.0_rt;
    } else {
        set_eb_data(i,j,apx,apy,dx,vfrac,vcent,barea,bcent,bnorm);
        // remove small cells
        if (vfrac(i,j,0) < small_volfrac) {
            set_covered(i,j,cell,vfrac,vcent,barea,bcent,bnorm);
            is_small_cell = true;
        }
    }
    return is_small_cell;
}

}

int build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                 Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                 Array4<Real const> const& levset,
                 Array4<Real const> const& interx, Array4<Real const> const& intery,
                 Array4<Real> const& apx, Array4<Real> const& apy,
                 Array4<Real> const& fcx, Array4<Real> const& fcy,
                 GpuArray<Real,AMREX_SPACEDIM> const& dx,
                 GpuArray<Real,AMREX_SPACEDIM> const& problo,
                 bool cover_multiple_cuts) noexcept
{
#ifdef AMREX_USE_FLOAT
    constexpr Real small = 1.e-5_rt;
#else
    constexpr Real small = 1.e-14;
#endif
    const Real dxinv = 1.0_rt/dx[0];
    const Real dyinv = 1.0_rt/dx[1];
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
                apx(i,j,0) = 1.0_rt;
                fcx(i,j,0) = 0.0_rt;
            } else if (fx(i,j,0) == Type::covered) {
                apx(i,j,0) = 0.0_rt;
                fcx(i,j,0) = 0.0_rt;
            } else {
                if (levset(i,j,0) < 0.0_rt) {
                    apx(i,j,0) = (intery(i,j,0)-(problo[1]+j*dx[1]))*dyinv;
                    fcx(i,j,0) = 0.5_rt*apx(i,j,0) - 0.5_rt;
                } else {
                    apx(i,j,0) = 1.0_rt - (intery(i,j,0)-(problo[1]+j*dx[1]))*dyinv;
                    fcx(i,j,0) = 0.5_rt - 0.5_rt*apx(i,j,0);
                }

                if (apx(i,j,0) > 1.0_rt-small) {
                    apx(i,j,0) = 1.0_rt;
                    fcx(i,j,0) = 0.0_rt;
                    fx(i,j,0) = Type::regular;
                } else if (apx(i,j,0) < small) {
                    apx(i,j,0) = 0.0_rt;
                    fcx(i,j,0) = 0.0_rt;
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
                apy(i,j,0) = 1.0_rt;
                fcy(i,j,0) = 0.0_rt;
            } else if (fy(i,j,0) == Type::covered) {
                apy(i,j,0) = 0.0_rt;
                fcy(i,j,0) = 0.0_rt;
            } else {
                if (levset(i,j,0) < 0.0_rt) {
                    apy(i,j,0) = (interx(i,j,0)-(problo[0]+i*dx[0]))*dxinv;
                    fcy(i,j,0) = 0.5_rt*apy(i,j,0) - 0.5_rt;
                } else {
                    apy(i,j,0) = 1.0_rt - (interx(i,j,0)-(problo[0]+i*dx[0]))*dxinv;
                    fcy(i,j,0) = 0.5_rt - 0.5_rt*apy(i,j,0);
                }

                if (apy(i,j,0) > 1.0_rt-small) {
                    apy(i,j,0) = 1.0_rt;
                    fcy(i,j,0) = 0.0_rt;
                    fy(i,j,0) = Type::regular;
                } else if (apy(i,j,0) < small) {
                    apy(i,j,0) = 0.0_rt;
                    fcy(i,j,0) = 0.0_rt;
                    fy(i,j,0) = Type::covered;
                }
            }
        }}
    });

    Gpu::Buffer<int> nmulticuts = {0};
    int* hp = nmulticuts.hostData();
    int* dp = nmulticuts.data();

    const Box& bxg1 = amrex::grow(bx, 1);
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        amrex::ignore_unused(k);
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
                    Gpu::Atomic::Add(dp,1);
                }
            }
        }
    });

    nmulticuts.copyToHost();

    if (*hp > 0 && !cover_multiple_cuts) {
        amrex::Abort("amrex::EB2::build_faces: more than 2 cuts not supported");
    }

    return *hp;
}

void build_cells (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Array4<Real> const& levset,
                  Real small_volfrac, Geometry const& geom, bool extend_domain_face,
                  int& nsmallcells, int const nmulticuts) noexcept
{
    Gpu::Buffer<int> smc = {0};
    int* hp = smc.hostData();
    int* dp = smc.data();

    const Box& bxg1 = amrex::grow(bx,1);
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bxg1, i, j, k,
    {
        amrex::ignore_unused(k);
        bool is_small = set_eb_cell(i, j, cell, apx, apy, dx, vfrac, vcent, barea, bcent,
                                    bnorm, small_volfrac);
        if (is_small) {
            Gpu::Atomic::Add(dp, 1);
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
                  Gpu::Atomic::Add(dp, 1);
                  set_covered(i,j,cell,vfrac,vcent,barea,bcent,bnorm);
              }
          });
       }
    }

    smc.copyToHost();
    nsmallcells += *hp;

    if (nsmallcells > 0 || nmulticuts > 0) {
        Box const& nbxg1 = amrex::surroundingNodes(bxg1);
        AMREX_HOST_DEVICE_FOR_3D(nbxg1, i, j, k,
        {
            if (levset(i,j,k) < Real(0.0)) {
                if ((bxg1.contains(i-1,j-1,k)
                     &&       cell(i-1,j-1,k).isCovered()) ||
                    (bxg1.contains(i  ,j-1,k)
                     &&       cell(i  ,j-1,k).isCovered()) ||
                    (bxg1.contains(i-1,j  ,k)
                     &&       cell(i-1,j  ,k).isCovered()) ||
                    (bxg1.contains(i  ,j  ,k)
                     &&       cell(i  ,j  ,k).isCovered()))
                {
                    levset(i,j,k) = Real(0.0);
                }
            }
        });
    }

    set_connection_flags(bxg1, cell, fx, fy);
}

void set_connection_flags (Box const& bxg1,
                           Array4<EBCellFlag> const& cell,
                           Array4<Type_t> const& fx, Array4<Type_t> const& fy) noexcept
{
    // Build neighbors.  By default, all neighbors are already set.
    AMREX_HOST_DEVICE_FOR_3D ( bxg1, i, j, k,
    {
        amrex::ignore_unused(k);

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

}
