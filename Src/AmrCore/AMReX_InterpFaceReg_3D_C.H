#ifndef AMREX_INTERP_FACE_REG_3D_C_H_
#define AMREX_INTERP_FACE_REG_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interp_face_reg (int i, int j, int k, IntVect const& rr, Array4<Real> const& fine, int scomp,
                      Array4<Real const> const& crse, Array4<Real> const& slope, int ncomp,
                      Box const& domface, int idim)
{
    int ic = amrex::coarsen(i,rr[0]);
    int jc = amrex::coarsen(j,rr[1]);
    int kc = amrex::coarsen(k,rr[2]);
    if (idim == 0) {
        if (jc == domface.smallEnd(1) || jc == domface.bigEnd(1)) {
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) = crse(ic,jc,kc,n);
            }
        } else {
            Real sfy = Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
                Real dc = Real(0.5) * (crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n));
                Real df = Real(2.0) * (crse(ic,jc+1,kc,n) - crse(ic,jc  ,kc,n));
                Real db = Real(2.0) * (crse(ic,jc  ,kc,n) - crse(ic,jc-1,kc,n));
                Real sy = (df*db >= Real(0.0)) ?
                    amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
                sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));
                if (dc != Real(0.0)) {
                    sfy = amrex::min(sfy, sy / dc);
                }
                slope(i,j,k,n) = dc;
            }
            Real yoff = (j - jc*rr[1] + Real(0.5)) / Real(rr[1]) - Real(0.5);
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) = crse(ic,jc,kc,n) + yoff * slope(i,j,k,n) * sfy;
            }
        }

        if (kc != domface.smallEnd(2) && kc != domface.bigEnd(2)) {
            Real sfz = Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
                Real dc = Real(0.5) * (crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n));
                Real df = Real(2.0) * (crse(ic,jc,kc+1,n) - crse(ic,jc,kc  ,n));
                Real db = Real(2.0) * (crse(ic,jc,kc  ,n) - crse(ic,jc,kc-1,n));
                Real sz = (df*db >= Real(0.0)) ?
                    amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
                sz = amrex::Math::copysign(Real(1.),dc)*amrex::min(sz,amrex::Math::abs(dc));
                if (dc != Real(0.0)) {
                    sfz = amrex::min(sfz, sz / dc);
                }
                slope(i,j,k,n) = dc;
            }
            Real zoff = (k - kc*rr[2] + Real(0.5)) / Real(rr[2]) - Real(0.5);
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) += zoff * slope(i,j,k,n) * sfz;
            }
        }
    } else if (idim == 1) {
        if (ic == domface.smallEnd(0) || ic == domface.bigEnd(0)) {
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) = crse(ic,jc,kc,n);
            }
        } else {
            Real sfx = Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
                Real dc = Real(0.5) * (crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n));
                Real df = Real(2.0) * (crse(ic+1,jc,kc,n) - crse(ic  ,jc,kc,n));
                Real db = Real(2.0) * (crse(ic  ,jc,kc,n) - crse(ic-1,jc,kc,n));
                Real sx = (df*db >= Real(0.0)) ?
                    amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
                sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));
                if (dc != Real(0.0)) {
                    sfx = amrex::min(sfx, sx / dc);
                }
                slope(i,j,k,n) = dc;
            }
            Real xoff = (i - ic*rr[0] + Real(0.5)) / Real(rr[0]) - Real(0.5);
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) = crse(ic,jc,kc,n) + xoff * slope(i,j,k,n) * sfx;
            }
        }

        if (kc != domface.smallEnd(2) && kc != domface.bigEnd(2)) {
            Real sfz = Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
                Real dc = Real(0.5) * (crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n));
                Real df = Real(2.0) * (crse(ic,jc,kc+1,n) - crse(ic,jc,kc  ,n));
                Real db = Real(2.0) * (crse(ic,jc,kc  ,n) - crse(ic,jc,kc-1,n));
                Real sz = (df*db >= Real(0.0)) ?
                    amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
                sz = amrex::Math::copysign(Real(1.),dc)*amrex::min(sz,amrex::Math::abs(dc));
                if (dc != Real(0.0)) {
                    sfz = amrex::min(sfz, sz / dc);
                }
                slope(i,j,k,n) = dc;
            }
            Real zoff = (k - kc*rr[2] + Real(0.5)) / Real(rr[2]) - Real(0.5);
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) += zoff * slope(i,j,k,n) * sfz;
            }
        }
    } else {
        if (ic == domface.smallEnd(0) || ic == domface.bigEnd(0)) {
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) = crse(ic,jc,kc,n);
            }
        } else {
            Real sfx = Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
                Real dc = Real(0.5) * (crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n));
                Real df = Real(2.0) * (crse(ic+1,jc,kc,n) - crse(ic  ,jc,kc,n));
                Real db = Real(2.0) * (crse(ic  ,jc,kc,n) - crse(ic-1,jc,kc,n));
                Real sx = (df*db >= Real(0.0)) ?
                    amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
                sx = amrex::Math::copysign(Real(1.),dc)*amrex::min(sx,amrex::Math::abs(dc));
                if (dc != Real(0.0)) {
                    sfx = amrex::min(sfx, sx / dc);
                }
                slope(i,j,k,n) = dc;
            }
            Real xoff = (i - ic*rr[0] + Real(0.5)) / Real(rr[0]) - Real(0.5);
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) = crse(ic,jc,kc,n) + xoff * slope(i,j,k,n) * sfx;
            }
        }

        if (jc != domface.smallEnd(1) && jc != domface.bigEnd(1)) {
            Real sfy = Real(1.0);
            for (int n = 0; n < ncomp; ++n) {
                Real dc = Real(0.5) * (crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n));
                Real df = Real(2.0) * (crse(ic,jc+1,kc,n) - crse(ic,jc  ,kc,n));
                Real db = Real(2.0) * (crse(ic,jc  ,kc,n) - crse(ic,jc-1,kc,n));
                Real sy = (df*db >= Real(0.0)) ?
                    amrex::min(amrex::Math::abs(df),amrex::Math::abs(db)) : Real(0.);
                sy = amrex::Math::copysign(Real(1.),dc)*amrex::min(sy,amrex::Math::abs(dc));
                if (dc != Real(0.0)) {
                    sfy = amrex::min(sfy, sy / dc);
                }
                slope(i,j,k,n) = dc;
            }
            Real yoff = (j - jc*rr[1] + Real(0.5)) / Real(rr[1]) - Real(0.5);
            for (int n = 0; n < ncomp; ++n) {
                fine(i,j,k,n+scomp) += yoff * slope(i,j,k,n) * sfy;
            }
        }
    }
}

}

#endif
