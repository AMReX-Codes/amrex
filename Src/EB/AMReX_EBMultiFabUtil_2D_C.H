#ifndef AMREX_EB_MULTIFAB_UTIL_2D_C_H_
#define AMREX_EB_MULTIFAB_UTIL_2D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real v)
{
    if (f(i-1,j-1,k).isCovered() and f(i  ,j-1,k).isCovered() and
        f(i-1,j  ,k).isCovered() and f(i  ,j  ,k).isCovered())
    {
        d(i,j,k,n+icomp) = v;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real const * AMREX_RESTRICT v)
{
    if (f(i-1,j-1,k).isCovered() and f(i  ,j-1,k).isCovered() and
        f(i-1,j  ,k).isCovered() and f(i  ,j  ,k).isCovered())
    {
        d(i,j,k,n+icomp) = v[n];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_with_vol (int i, int j, int k,
                          Array4<Real const> const& fine, int fcomp,
                          Array4<Real> const& crse, int ccomp,
                          Array4<Real const> const& fv, Array4<Real const> const& vfrc,
                          Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real cv = 0.0;
        constexpr int kk = 0;
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = fv(ii,jj,kk)*vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,kk,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown (int i, int j, int k,
                 Array4<Real const> const& fine, int fcomp,
                 Array4<Real> const& crse, int ccomp,
                 Array4<Real const> const& vfrc,
                 Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real cv = 0.0;
        constexpr int kk = 0;
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,kk,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_x (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int ii = i*ratio.x;
    constexpr int kk = 0;
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real a = 0.0;
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }
        if (a > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(ii,j*ratio.y,kk,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_y (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int jj = j*ratio.y;
    constexpr int kk = 0;
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real a = 0.0;
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }
        if (a > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,jj,kk,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_boundaries (int i, int j, int k,
                            Array4<Real const> const& fine, int fcomp,
                            Array4<Real> const& crse, int ccomp,
                            Array4<Real const> const& ba,
                            Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real cv = 0.0;
        constexpr int kk = 0;
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = ba(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = 0.0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_compute_divergence (int i, int j, int k, int n, Array4<Real> const& divu,
                            Array4<Real const> const& u, Array4<Real const> const& v,
                            Array4<int const> const& ccm, Array4<EBCellFlag const> const& flag,
                            Array4<Real const> const& vfrc, Array4<Real const> const& apx,
                            Array4<Real const> const& apy, Array4<Real const> const& fcx,
                            Array4<Real const> const& fcy, GpuArray<Real,2> const& dxinv,
                            bool already_on_centroids)
{
    if (flag(i,j,k).isCovered())
    {
        divu(i,j,k,n) = 0.0;
    }
    else if (flag(i,j,k).isRegular())
    {
        divu(i,j,k,n) = dxinv[0] * (u(i+1,j,k,n)-u(i,j,k,n))
            +           dxinv[1] * (v(i,j+1,k,n)-v(i,j,k,n));
    }
    else if (already_on_centroids)
    {
        divu(i,j,k,n) = (1.0/vfrc(i,j,k)) *
            ( dxinv[0] * (apx(i+1,j,k)*u(i+1,j,k,n)-apx(i,j,k)*u(i,j,k,n))
            + dxinv[1] * (apy(i,j+1,k)*v(i,j+1,k,n)-apy(i,j,k)*v(i,j,k,n)) );
    }
    else
    {
        Real fxm = u(i,j,k,n);
        if (apx(i,j,k) != 0.0 and apx(i,j,k) != 1.0) {
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i,j,k)));
            Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? amrex::Math::abs(fcx(i,j,k)) : 0.0;
            fxm = (1.0-fracy)*fxm + fracy*u(i,jj,k,n);
        }

        Real fxp = u(i+1,j,k,n);
        if (apx(i+1,j,k) != 0.0 and apx(i+1,j,k) != 1.0) {
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i+1,j,k)));
            Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? amrex::Math::abs(fcx(i+1,j,k)) : 0.0;
            fxp = (1.0-fracy)*fxp + fracy*u(i+1,jj,k,n);
        }

        Real fym = v(i,j,k,n);
        if (apy(i,j,k) != 0.0 and apy(i,j,k) != 1.0) {
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k)));
            Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? amrex::Math::abs(fcy(i,j,k)) : 0.0;
            fym = (1.0-fracx)*fym + fracx*v(ii,j,k,n);
        }

        Real fyp = v(i,j+1,k,n);
        if (apy(i,j+1,k) != 0.0 and apy(i,j+1,k) != 1.0) {
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j+1,k)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? amrex::Math::abs(fcy(i,j+1,k)) : 0.0;
            fyp = (1.0-fracx)*fyp + fracx*v(ii,j+1,k,n);
        }

        divu(i,j,k,n) = (1.0/vfrc(i,j,k)) *
            ( dxinv[0] * (apx(i+1,j,k)*fxp-apx(i,j,k)*fxm)
            + dxinv[1] * (apy(i,j+1,k)*fyp-apy(i,j,k)*fym) );
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avg_fc_to_cc (int i, int j, int k, int n, Array4<Real> const& cc,
                      Array4<Real const> const& fx, Array4<Real const> const& fy,
                      Array4<Real const> const& ax, Array4<Real const> const& ay,
                      Array4<EBCellFlag const> const& flag)
{
    if (flag(i,j,k).isCovered()) {
        cc(i,j,k,n+0) = 0.0;
        cc(i,j,k,n+1) = 0.0;
    } else {
        if (ax(i,j,k) == 0.0) {
            cc(i,j,k,n+0) = fx(i+1,j,k);
        } else if (ax(i+1,j,k) == 0.0) {
            cc(i,j,k,n+0) = fx(i,j,k);
        } else {
            cc(i,j,k,n+0) = 0.5 * (fx(i,j,k) + fx(i+1,j,k));
        }

        if (ay(i,j,k) == 0.0) {
            cc(i,j,k,n+1) = fy(i,j+1,k);
        } else if (ay(i,j+1,k) == 0.0) {
            cc(i,j,k,n+1) = fy(i,j,k);
        } else {
            cc(i,j,k,n+1) = 0.5 * (fy(i,j,k) + fy(i,j+1,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2cent (Box const& box,
                        const Array4<Real>& phicent,
                        Array4<Real const> const& phicc,
                        Array4<EBCellFlag const> const& flag,
                        Array4<Real const> const& cent,
                        int ncomp) noexcept
{
  amrex::Loop(box, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (flag(i,j,k).isCovered())
    {
      phicent(i,j,k,n) = phicc(i,j,k,n); //do nothing
    }
    else
    {
      if (flag(i,j,k).isRegular())
      {
        phicent(i,j,k,n) = phicc(i,j,k,n);
      }
      else
      {
        Real gx = cent(i,j,k,0);
        Real gy = cent(i,j,k,1);
        int ii = (gx < 0.0) ? i - 1 : i + 1;
        int jj = (gy < 0.0) ? j - 1 : j + 1;
        gx = amrex::Math::abs(gx);
        gy = amrex::Math::abs(gy);
        Real gxy = gx*gy;

        phicent(i,j,k,n) = ( 1.0 - gx - gy + gxy ) * phicc(i ,j ,k ,n)
          +                (            gy - gxy ) * phicc(i ,jj,k ,n)
          +                (       gx      - gxy ) * phicc(ii,j ,k ,n)
          +                (                 gxy ) * phicc(ii,jj,k ,n);
      }
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2facecent (Box const& box,
                            Array4<Real const> const& phi,
                            Array4<Real const> const& apx, Array4<Real const> const& apy,
                            Array4<Real const> const& fcx,
                            Array4<Real const> const& fcy,
                            Array4<Real> const& edg_x,
                            Array4<Real> const& edg_y,
                            int ncomp,
                            const Box& domain,
                            const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);

  const Box ubx = amrex::surroundingNodes(box,0);
  const Box vbx = amrex::surroundingNodes(box,1);

  //
  // ===================== X =====================
  //
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
      if (apx(i,j,k) == 0)
      {
        edg_x(i,j,k,n) = 1e40;
      }
      else
      {
        if (apx(i,j,k) == 1)
        {
          if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
          {
            edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
          }
          else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
          {
            edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
          }
          else
          {
            edg_x(i,j,k,n) = 0.5 * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
          }
        }
        else
        {
          int ii = i - 1;
          int jj = (fcx(i,j,k) < 0.0) ? j - 1 : j + 1;
          Real gx = 0.5;
          Real gy = amrex::Math::abs(fcx(i,j,k));
          Real gxy = gx*gy;
          edg_x(i,j,k,n) = ( 1.0 - gx - gy + gxy ) * phi(i ,j ,k ,n)
           +               (            gy - gxy ) * phi(i ,jj,k ,n)
           +               (       gx      - gxy ) * phi(ii,j ,k ,n)
           +               (                 gxy ) * phi(ii,jj,k ,n);
         }
       }
    });
  

    //
  // ===================== Y =====================
  //
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
      if (apy(i,j,k) == 0)
      {
        edg_y(i,j,k,n) = 1e40;
      }
      else
      {
        if (apy(i,j,k) == 1)
        {
          if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
          {
            edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
          }
          else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
          {
            edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
          }
          else
          {
            edg_y(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
          }
        }
        else
        {
          int ii = (fcy(i,j,k) < 0.0) ? i - 1 : i + 1;
          int jj = j - 1;
          Real gx = amrex::Math::abs(fcy(i,j,k));
          Real gy = 0.5;
          Real gxy = gx*gy;
          edg_y(i,j,k,n) = ( 1.0 - gx - gy + gxy ) * phi(i ,j ,k ,n)
           +               (            gy - gxy ) * phi(i ,jj,k ,n)
           +               (       gx      - gxy ) * phi(ii,j ,k ,n)
           +               (                 gxy ) * phi(ii,jj,k ,n);
         }
       }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_centroid2facecent (Box const& box,
                                  Array4<Real const> const& phi,
                                  Array4<Real const> const& apx, Array4<Real const> const& apy,
                                  Array4<Real const> const& cvol,
                                  Array4<Real const> const& ccent,
                                  Array4<Real const> const& fcx,
                                  Array4<Real const> const& fcy,
                                  Array4<Real> const& edg_x,
                                  Array4<Real> const& edg_y,
                                  int ncomp,
                                  const Box& domain,
                                  const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);

  const Box ubx = amrex::surroundingNodes(box,0);
  const Box vbx = amrex::surroundingNodes(box,1);

  //
  // ===================== X =====================
  //
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
      if (apx(i,j,k) == 0)
      {
          edg_x(i,j,k,n) = 1e40;
      }
      else if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
      {
          edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
      }
      else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
      {
          edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
      }
      else if (apx(i,j,k) == 1 && cvol(i,j,k) == 1 && cvol(i-1,j,k) == 1)
      {
          edg_x(i,j,k,n) = 0.5 * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
      } 
      else 
      {
          int ii = i - 1;
          int jj = (fcx(i,j,k) < 0.0) ? j - 1 : j + 1;

          Real d0    = 0.5 + ccent( i,j,k,0);                                 // distance in x-dir from centroid(i  ,j) to face(i,j)
          Real d1    = 0.5 - ccent(ii,j,k,0);                                 // distance in x-dir from centroid(i-1,j) to face(i,j)

          Real d0_jj = 0.5 + ccent( i,jj,k,0);                                 // distance in x-dir from centroid(i  ,jj) to face(i,jj)
          Real d1_jj = 0.5 - ccent(ii,jj,k,0);                                 // distance in x-dir from centroid(i-1,jj) to face(i,jj)

          Real a0    = d1 / (d0 + d1);                                         // coefficient of phi(i  ,j,k,n)
          Real a1    = d0 / (d0 + d1);                                         // coefficient of phi(i-1,j,k,n)

          Real a0_jj = d1_jj / (d0_jj + d1_jj);                                // coefficient of phi(i  ,jj,k,n)
          Real a1_jj = d0_jj / (d0_jj + d1_jj);                                // coefficient of phi(i-1,jj,k,n)

          Real phi01    = a0    *   phi(i, j,k,n)  + a1    *   phi(ii, j,k,n); // interpolation of phi   in x-direction
          Real y01      = a0    * ccent(i, j,k,1)  + a1    * ccent(ii, j,k,1); // interpolation of y-loc in x-direction

          Real phi01_jj = a0_jj *   phi(i,jj,k,n)  + a1_jj *   phi(ii,jj,k,n); // interpolation of phi   in x-direction
          Real y01_jj   = a0_jj * ccent(i,jj,k,1)  + a1_jj * ccent(ii,jj,k,1); // interpolation of y-loc in x-direction
          
          edg_x(i,j,k,n) = ( (     fcx(i,j,k) - y01   ) * phi01_jj + 
                             (1. - fcx(i,j,k) + y01_jj) * phi01    ) / (1. - y01 + y01_jj);
      }
    });
  

  //
  // ===================== Y =====================
  //
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == 0)
        {
            edg_y(i,j,k,n) = 1e40;
        }
        else if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
        {
            edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
        }
        else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
        {
            edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
        }
        else if (apy(i,j,k) == 1 && cvol(i,j,k) == 1 && cvol(i,j-1,k) == 1)
        {
            edg_y(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
        }
        else 
        {
          int jj = j - 1;
          int ii = (fcy(i,j,k) < 0.0) ? i - 1 : i + 1;

          Real d0    = 0.5 + ccent(i, j,k,1);                                 // distance in y-dir from centroid(i,j  ) to face(i,j)
          Real d1    = 0.5 - ccent(i,jj,k,1);                                 // distance in y-dir from centroid(i,j-1) to face(i,j)

          Real d0_ii = 0.5 + ccent(ii, j,k,0);                                 // distance from centroid(ii,j  ) to face(ii,j)
          Real d1_ii = 0.5 - ccent(ii,jj,k,0);                                 // distance from centroid(ii,j-1) to face(ii,j)

          Real a0    = d1 / (d0 + d1);                                         // coefficient of phi(i,j  ,k,n)
          Real a1    = d0 / (d0 + d1);                                         // coefficient of phi(i,j-1,k,n)

          Real a0_ii = d1_ii / (d0_ii + d1_ii);                                // coefficient of phi(i  ,jj,k,n)
          Real a1_ii = d0_ii / (d0_ii + d1_ii);                                // coefficient of phi(i-1,jj,k,n)

          Real phi01    = a0    *   phi( i,j,k,n)  + a1    *   phi( i,jj,k,n); // interpolation of phi   in y-direction
          Real x01      = a0    * ccent( i,j,k,0)  + a1    * ccent( i,jj,k,0); // interpolation of x-loc in y-direction

          Real phi01_ii = a0_ii *   phi(ii,j,k,n)  + a1_ii *   phi(ii,jj,k,n); // interpolation of phi   in y-direction
          Real x01_ii   = a0_ii * ccent(ii,j,k,0)  + a1_ii * ccent(ii,jj,k,0); // interpolation of x-loc in y-direction
          
          edg_y(i,j,k,n) = ( (     fcy(i,j,k) - x01   ) * phi01_ii + 
                             (1. - fcy(i,j,k) + x01_ii) * phi01    ) / (1. - x01 + x01_ii);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2face (Box const& box,
                        Array4<Real const> const& phi,
                        Array4<Real> const& edg_x,
                        Array4<Real> const& edg_y,
                        int ncomp,
                        const Box& domain,
                        const BCRec* bc) noexcept
{
  
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);
  
  const Box ubx = amrex::surroundingNodes(box,0);
  const Box vbx = amrex::surroundingNodes(box,1);
  
  //
  // ===================== X =====================
  //
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
        {
          edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
        }
        else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
        {
          edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
        }
        else
        {
          edg_x(i,j,k,n) = 0.5 * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
        } 
    });
  

  //
  // ===================== Y =====================
  //
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
        {
          edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
        }
        else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
        {
          edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
        }
        else
        {
          edg_y(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
        }
    });
}
}

#endif
