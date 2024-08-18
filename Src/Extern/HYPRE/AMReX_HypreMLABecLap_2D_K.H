#ifndef AMREX_HYPRE_ML_ABECLAP_2D_K_H_
#define AMREX_HYPRE_ML_ABECLAP_2D_K_H_

#include <AMReX_Array.H>
#include <AMReX_Orientation.H>

namespace amrex {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_f2c_set_values (IntVect const& cell, Real* values,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dx, Real sb,
                                  GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& b,
                                  GpuArray<Array4<int const>,AMREX_SPACEDIM*2> const& bmask,
                                  IntVect const& refratio, int not_covered)
{
    Array2D<Real,-1,1,-1,1> tmp;
    for (auto& x : tmp) { x = Real(0.0); }

    Array2D<int,-1,1,-1,1> used;
    for (auto& x : used) { x = 0; }

    for (OrientationIter ori; ori; ++ori) {
        auto const face = ori();
        int const idir = face.coordDir();
        int const idir1 = 1-idir; // the transverse direction
        IntVect offset(0);
        offset[idir] = face.isLow() ? -1 : +1;
        IntVect const cell_out = cell + offset;
        auto const& msk = bmask[face];
        if (msk.contains(cell_out) && msk(cell_out) == not_covered) {
            // There is a coarse cell on the other side of the face. There
            // are three cases for the coarse cells involved. (1)
            // Interpolation using 3 coarse cells. (2) Upward biased
            // interpolation using 2 coarse cells. (3) Doward biased
            // interpolation using 2 coarse cells. (Here up and down means
            // the y-dirction, if we assume idir is the x-direction.)
            IntVect offset_t(0);
            offset_t[idir1] = refratio[idir1];
            Real bcoeff = b[idir] ? b[idir](face.isLow() ? cell : cell_out) : Real(1.0); // b is on face
            Real poly_coef[3];
            {
                Real xx[3] = {Real(-0.5)*Real(refratio[idir]), Real(0.5), Real(1.5)};
                poly_interp_coeff<3>(Real(-0.5), xx, poly_coef);
            }
            Real fac = -(sb / (dx[idir]*dx[idir])) * bcoeff * poly_coef[0];
            int const rr1 = refratio[idir1];
            int const i1 = cell[idir1];
            int const i1c = amrex::coarsen(i1, rr1);
            Real xInt = Real(-0.5) + (i1-i1c*rr1+Real(0.5))/Real(rr1);
            Real xc[] = {Real(-1.0), Real(0.0), Real(1.0)};
            Real c[] = {Real(0.0), Real(0.0), Real(0.0)};
            int cc[] = {0, 0, 0};
            if (msk(cell_out-offset_t) == not_covered &&
                msk(cell_out+offset_t) == not_covered)
            {
                poly_interp_coeff<3>(xInt, xc, c);
                cc[0] = cc[1] = cc[2] = 1;
            } else if (msk(cell_out+offset_t) == not_covered) {
                poly_interp_coeff<2>(xInt, &(xc[1]), &(c[1]));
                cc[1] = cc[2] = 1;
            } else {
                poly_interp_coeff<2>(xInt, xc, c);
                cc[0] = cc[1] = 1;
            }
            if (face == Orientation(0, Orientation::low)) {
                for (int m = 0; m < 3; ++m) {
                    tmp(-1,m-1) += c[m] * fac;
                    used(-1,m-1) += cc[m];
                }
            } else if (face == Orientation(0, Orientation::high)) {
                for (int m = 0; m < 3; ++m) {
                    tmp(1,m-1) += c[m] * fac;
                    used(1,m-1) += cc[m];
                }
            } else if (face == Orientation(1, Orientation::low)) {
                for (int m = 0; m < 3; ++m) {
                    tmp(m-1,-1) += c[m] * fac;
                    used(m-1,-1) += cc[m];
                }
            } else if (face == Orientation(1, Orientation::high)) {
                for (int m = 0; m < 3; ++m) {
                    tmp(m-1,1) += c[m] * fac;
                    used(m-1,1) += cc[m];
                }
            }
        }
    }

    auto const* ptmp = tmp.begin();
    auto const* pused = used.begin();
    for (int m = 0; m < 9; ++m) {
        if (pused[m]) {
            (*values) += ptmp[m];
            ++values;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_c2f (int i, int j, int k,
                       Array4<GpuArray<Real,2*AMREX_SPACEDIM+1>> const& stencil,
                       GpuArray<HYPRE_Int,AMREX_SPACEDIM>* civ, HYPRE_Int* nentries,
                       int* entry_offset, Real* entry_values,
                       Array4<int const> const& offset_from,
                       Array4<int const> const& nentries_to,
                       Array4<int const> const& offset_to,
                       GpuArray<Real,AMREX_SPACEDIM> const& dx, Real sb,
                       Array4<int const> const& offset_bx,
                       Array4<int const> const& offset_by,
                       Real const* bx, Real const* by,
                       Array4<int const> const& fine_mask,
                       IntVect const& rr, Array4<int const> const& crse_mask)
{
    if (fine_mask(i,j,k)) {
        // Let's set off-diagonal elements to zero
        for (int m = 1; m < 2*AMREX_SPACEDIM+1; ++m) {
            stencil(i,j,k)[m] = Real(0.0);
        }
    } else if (nentries_to(i,j,k) > 0) {
        int const fromoff = offset_from(i,j,k);
        civ[fromoff][0] = i;
        civ[fromoff][1] = j;
        nentries[fromoff] = nentries_to(i,j,k);
        int foff = offset_to(i,j,k);
        entry_offset[fromoff] = foff;

        // We must iterate the faces in the lexicographical order of fine
        // neighbor cells, because that's the order when non-stencil entries
        // were added to Hypre's graph. Also note that a coarse cell will
        // not have entries to fine cells at both ends of a direction. Thus
        // we do not have to worry about the order between fine cells at the
        // small and big ends of the same direction.

        if (fine_mask(i,j-1,k)) {
            stencil(i,j,k)[0] += stencil(i,j,k)[3];
            stencil(i,j,k)[3] = Real(0.0);
            // Reflux: sb/h^2*by(i,j,k)*(phi(i,j,k)-phi(i,j-1,k)) is replaced by
            // sb/h*sum_{fine faces}(dphi/dy*by)/n_fine_faces
            Real dyf = dx[1] / Real(rr[1]);
            Real dycinv = Real(1.0) / dx[1];
            Real dyfinv = Real(1.0) / dyf;
            Real cc[3];
            Real yy[3] = {dx[1]*Real(-0.5), dyf*Real(0.5), dyf*Real(1.5)};
            poly_interp_coeff<3>(dyf*Real(-0.5), yy, cc);
            for (int irx = 0; irx < rr[0]; ++irx) {
                Real bym = by ? by[offset_by(i,j,k)+irx] : Real(1.0);
                Real fac = sb*dycinv*dyfinv*bym/Real(rr[0]);
                // int ii = i*rr[0] + irx
                // int jj = j*rr[1]
                // dphi/dy = (phi_interp - phi_fine(jj-1)) / dy_fine
                //         = (phi_coarse*cc[0] + phi_fine(jj-1)*(cc[1]-1)
                //                             + phi_fine(jj-2)* cc[2]) / dy_fine
                // So the entry for fine cell (jj-1) is (cc[1]-1)*fac
                //                  fine cell (jj-2) is  cc[2]   *fac
                entry_values[foff+irx      ] += fac* cc[2];
                entry_values[foff+irx+rr[0]] += fac*(cc[1]-Real(1.0));

                // The coarse cell's stencils need updates too.
                Real xInt = Real(-0.5) + (irx+Real(0.5))/Real(rr[0]);
                Real xc[3] = {Real(-1.0), Real(0.0), Real(1.0)};
                Real ct[3] = {Real(0.0), Real(0.0), Real(0.0)};
                if ( fine_mask(i-1,j,k) ||
                    !crse_mask(i-1,j,k))
                {
                    poly_interp_coeff<2>(xInt, &(xc[1]), &(ct[1]));
                } else if ( fine_mask(i+1,j,k) ||
                           !crse_mask(i+1,j,k))
                {
                    poly_interp_coeff<2>(xInt, xc, ct);
                } else {
                    poly_interp_coeff<3>(xInt, xc, ct);
                }
                // phi_coarse = ct[0]*phi(i-1) + ct[1]*phi(i) + ct[2]*phi(i+1)
                stencil(i,j,k)[0] += (fac*cc[0])*ct[1];
                stencil(i,j,k)[1] += (fac*cc[0])*ct[0];
                stencil(i,j,k)[2] += (fac*cc[0])*ct[2];
            }
            foff += 2*rr[0];
        }

        if (fine_mask(i-1,j,k)) {
            stencil(i,j,k)[0] += stencil(i,j,k)[1];
            stencil(i,j,k)[1] = Real(0.0);
            // Reflux: sb/h^2*bx(i,j,k)*(phi(i,j,k)-phi(i-1,j,k)) is replaced by
            // sb/h*sum_{fine faces}(dphi/dx*bx).
            Real dxf = dx[0] / Real(rr[0]);
            Real dxcinv = Real(1.0) / dx[0];
            Real dxfinv = Real(1.0) / dxf;
            Real cc[3];
            Real xx[3] = {dx[0]*Real(-0.5), dxf*Real(0.5), dxf*Real(1.5)};
            poly_interp_coeff<3>(dxf*Real(-0.5), xx, cc);
            for (int iry = 0; iry < rr[1]; ++iry) {
                Real bxm = bx ? bx[offset_bx(i,j,k)+iry] : Real(1.0);
                Real fac = sb*dxcinv*dxfinv*bxm/Real(rr[1]);
                // int ii = i*rr[0]
                // int jj = j*rr[1] + iry
                // dphi/dx = (phi_interp - phi_fine(ii-1)) / dx_fine
                //         = (phi_coarse*cc[0] + phi_fine(ii-1)*(cc[1]-1)
                //                             + phi_fine(ii-2)* cc[2]) / dx_fine
                // So the entry for fine cell(ii-1) is (cc[1]-1)*fac
                //                  fine cell(ii-2) is  cc[2]   *fac
                entry_values[foff++] = fac* cc[2];
                entry_values[foff++] = fac*(cc[1] - Real(1.0));

                // The coarse cell's stencils need updates too.
                Real yInt = Real(-0.5) + (iry+Real(0.5))/Real(rr[1]);
                Real yc[3] = {Real(-1.0), Real(0.0), Real(1.0)};
                Real ct[3] = {Real(0.0), Real(0.0), Real(0.0)};
                if ( fine_mask(i,j-1,k) ||
                    !crse_mask(i,j-1,k))
                {
                    poly_interp_coeff<2>(yInt, &(yc[1]), &(ct[1]));
                } else if ( fine_mask(i,j+1,k) ||
                           !crse_mask(i,j+1,k))
                {
                    poly_interp_coeff<2>(yInt, yc, ct);
                } else {
                    poly_interp_coeff<3>(yInt, yc, ct);
                }
                // phi_coarse = ct[0]*phi(j-1) + ct[1]*phi(j) + ct[2]*phi(j+1)
                stencil(i,j,k)[0] += (fac*cc[0])*ct[1];
                stencil(i,j,k)[3] += (fac*cc[0])*ct[0];
                stencil(i,j,k)[4] += (fac*cc[0])*ct[2];
            }
        }

        if (fine_mask(i+1,j,k)) {
            stencil(i,j,k)[0] += stencil(i,j,k)[2];
            stencil(i,j,k)[2] = Real(0.0);
            // Reflux: sb/h^2*bx(i+1,j,k)*(phi(i,j,k)-phi(i+1,j,k)) is replaced by
            // sb/h*sum_{fine faces}(-dphi/dx*bx).
            Real dxf = dx[0] / Real(rr[0]);
            Real dxcinv = Real(1.0) / dx[0];
            Real dxfinv = Real(1.0) / dxf;
            Real cc[3];
            Real xx[3] = {dx[0]*Real(-0.5), dxf*Real(0.5), dxf*Real(1.5)};
            poly_interp_coeff<3>(dxf*Real(-0.5), xx, cc);
            for (int iry = 0; iry < rr[1]; ++iry) {
                Real bxp = bx ? bx[offset_bx(i+1,j,k)+iry] : Real(1.0);
                Real fac = sb*dxcinv*dxfinv*bxp/Real(rr[1]);
                // int ii = i*rr[0] + (rr[0]-1)
                // int jj = j*rr[1] + iry
                // -dphi/dx = (phi_interp - phi_fine(ii+1)) / dx_fine
                //          = (phi_coarse*cc[0] + phi_fine(ii+1)*(cc[1]-1)
                //                              + phi_fine(ii+2)* cc[2]) / dx_fine
                // So the entry for fine cell(ii+1) is (cc[1]-1)*fac
                //                  fine cell(ii+2) is  cc[2]   *fac
                entry_values[foff++] = fac*(cc[1] - Real(1.0));
                entry_values[foff++] = fac* cc[2];

                // The coarse cell's stencils need updates too.
                Real yInt = Real(-0.5) + (iry+Real(0.5))/Real(rr[1]);
                Real yc[3] = {Real(-1.0), Real(0.0), Real(1.0)};
                Real ct[3] = {Real(0.0), Real(0.0), Real(0.0)};
                if ( fine_mask(i,j-1,k) ||
                    !crse_mask(i,j-1,k))
                {
                    poly_interp_coeff<2>(yInt, &(yc[1]), &(ct[1]));
                } else if ( fine_mask(i,j+1,k) ||
                           !crse_mask(i,j+1,k))
                {
                    poly_interp_coeff<2>(yInt, yc, ct);
                } else {
                    poly_interp_coeff<3>(yInt, yc, ct);
                }
                // phi_coarse = ct[0]*phi(j-1) + ct[1]*phi(j) + ct[2]*phi(j+1)
                stencil(i,j,k)[0] += (fac*cc[0])*ct[1];
                stencil(i,j,k)[3] += (fac*cc[0])*ct[0];
                stencil(i,j,k)[4] += (fac*cc[0])*ct[2];
            }
        }

        if (fine_mask(i,j+1,k)) {
            stencil(i,j,k)[0] += stencil(i,j,k)[4];
            stencil(i,j,k)[4] = Real(0.0);
            // Reflux: sb/h^2*by(i,j+1,k)*(phi(i,j,k)-phi(i,j+1,k)) is replaced by
            // sb/h*sum_{fine faces}(-dphi/dy*by)
            Real dyf = dx[1] / Real(rr[1]);
            Real dycinv = Real(1.0) / dx[1];
            Real dyfinv = Real(1.0) / dyf;
            Real cc[3];
            Real yy[3] = {dx[1]*Real(-0.5), dyf*Real(0.5), dyf*Real(1.5)};
            poly_interp_coeff<3>(dyf*Real(-0.5), yy, cc);
            for (int irx = 0; irx < rr[0]; ++irx) {
                Real byp = by ? by[offset_by(i,j+1,k)+irx] : Real(1.0);
                Real fac = sb*dycinv*dyfinv*byp/Real(rr[0]);
                // int ii = i*rr[0] + irx
                // int jj = j*rr[1] + (rr[1]-1)
                // -dphi/dy = (phi_interp - phi_fine(jj+1)) / dy_fine
                //          = (phi_coarse*cc[0] + phi_fine(jj+1)*(cc[1]-1)
                //                              + phi_fine(jj+2)* cc[2]) / dy_fine
                // So the entry for fine cell (jj+1) is (cc[1]-1)*fac
                //                  fine cell (jj+2) is  cc[2]   *fac
                entry_values[foff+irx      ] += fac*(cc[1]-Real(1.0));
                entry_values[foff+irx+rr[0]] += fac* cc[2];

                // The coarse cell's stencils need updates too.
                Real xInt = Real(-0.5) + (irx+Real(0.5))/Real(rr[0]);
                Real xc[3] = {Real(-1.0), Real(0.0), Real(1.0)};
                Real ct[3] = {Real(0.0), Real(0.0), Real(0.0)};
                if ( fine_mask(i-1,j,k) ||
                    !crse_mask(i-1,j,k))
                {
                    poly_interp_coeff<2>(xInt, &(xc[1]), &(ct[1]));
                } else if ( fine_mask(i+1,j,k) ||
                           !crse_mask(i+1,j,k))
                {
                    poly_interp_coeff<2>(xInt, xc, ct);
                } else {
                    poly_interp_coeff<3>(xInt, xc, ct);
                }
                // phi_coarse = ct[0]*phi(i-1) + ct[1]*phi(i) + ct[2]*phi(i+1)
                stencil(i,j,k)[0] += (fac*cc[0])*ct[1];
                stencil(i,j,k)[1] += (fac*cc[0])*ct[0];
                stencil(i,j,k)[2] += (fac*cc[0])*ct[2];
            }
            // not needed: foff += 2*rr[0];
        }
    }
}

}

#endif
