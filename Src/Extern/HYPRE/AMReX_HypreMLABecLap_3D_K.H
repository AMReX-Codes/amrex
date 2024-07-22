#ifndef AMREX_HYPRE_ML_ABECLAP_3D_K_H_
#define AMREX_HYPRE_ML_ABECLAP_3D_K_H_

namespace amrex {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_f2c_set_values (IntVect const& cell, Real* values,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dx, Real sb,
                                  GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& b,
                                  GpuArray<Array4<int const>,AMREX_SPACEDIM*2> const& bmask,
                                  IntVect const& refratio, int not_covered)
{
    Array3D<Real,-1,1,-1,1,-1,1> tmp;
    for (auto& x : tmp) { x = Real(0.0); }

    Array3D<bool,-1,1,-1,1,-1,1> used;
    for (auto& x : used) { x = false; }

    for (OrientationIter ori; ori; ++ori) {
        auto const face = ori();
        int const idir = face.coordDir();
        int const idir1 = (idir+1 < AMREX_SPACEDIM) ? idir+1 : idir+1-AMREX_SPACEDIM;
        int const idir2 = (idir+2 < AMREX_SPACEDIM) ? idir+2 : idir+2-AMREX_SPACEDIM;
        IntVect offset(0);
        offset[idir] = face.isLow() ? -1 : +1;
        IntVect const cell_out = cell + offset;
        auto const& msk = bmask[face];
        if (msk.contains(cell_out) && msk(cell_out) == not_covered) {
            // There is a coarse cell on the other side of the face.
            int const rr1 = refratio[idir1];
            int const rr2 = refratio[idir2];
            IntVect offset_t1(0);
            IntVect offset_t2(0);
            IntVect offset_tr1(0);
            IntVect offset_tr2(0);
            offset_t1 [idir1] = 1;
            offset_t2 [idir2] = 1;
            offset_tr1[idir1] = rr1;
            offset_tr2[idir2] = rr2;

            Real bcoeff = b[idir] ? b[idir](face.isLow() ? cell : cell_out) : Real(1.0); // b is on face
            Real poly_coef[3];
            {
                Real xx[3] = {Real(-0.5)*Real(refratio[idir]), Real(0.5), Real(1.5)};
                poly_interp_coeff<3>(Real(-0.5), xx, poly_coef);
            }
            Real fac = -(sb / (dx[idir]*dx[idir])) * bcoeff * poly_coef[0];

            used(offset[0],offset[1],offset[2]) = true;
            tmp (offset[0],offset[1],offset[2]) = fac;

            int it  = cell[idir1];
            int itc = amrex::coarsen(it, rr1);
            Real const xt1 = Real(-0.5) + (it-itc*rr1+Real(0.5))/Real(rr1);
            it = cell[idir2];
            itc = amrex::coarsen(it, rr2);
            Real const xt2 = Real(-0.5) + (it-itc*rr2+Real(0.5))/Real(rr2);

            if (msk(cell_out-offset_tr1) == not_covered &&
                msk(cell_out+offset_tr1) == not_covered)
            {
                IntVect iv = offset - offset_t1;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += fac*(Real(-0.5)*xt1 + Real(0.5)*xt1*xt1);

                iv = offset + offset_t1;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += fac*(Real(0.5)*xt1 + Real(0.5)*xt1*xt1);

                tmp(offset[0],offset[1],offset[2]) -= fac*(xt1*xt1);
            }
            else if (msk(cell_out+offset_tr1) == not_covered)
            {
                IntVect iv = offset + offset_t1;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += fac*xt1;

                tmp(offset[0],offset[1],offset[2]) -= fac*xt1;
            }
            else
            {
                IntVect iv = offset - offset_t1;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) -= fac*xt1;

                tmp(offset[0],offset[1],offset[2]) += fac*xt1;
            }

            if (msk(cell_out-offset_tr2) == not_covered &&
                msk(cell_out+offset_tr2) == not_covered)
            {
                IntVect iv = offset - offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += fac*(Real(-0.5)*xt2 + Real(0.5)*xt2*xt2);

                iv = offset + offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += fac*(Real(0.5)*xt2 + Real(0.5)*xt2*xt2);

                tmp(offset[0],offset[1],offset[2]) -= fac*(xt2*xt2);
            }
            else if (msk(cell_out+offset_tr2) == not_covered)
            {
                IntVect iv = offset + offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += fac*xt2;

                tmp(offset[0],offset[1],offset[2]) -= fac*xt2;
            }
            else
            {
                IntVect iv = offset - offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) -= fac*xt2;

                tmp(offset[0],offset[1],offset[2]) += fac*xt2;
            }

            if (msk(cell_out-offset_tr1-offset_tr2) == not_covered &&
                msk(cell_out+offset_tr1-offset_tr2) == not_covered &&
                msk(cell_out-offset_tr1+offset_tr2) == not_covered &&
                msk(cell_out+offset_tr1+offset_tr2) == not_covered)
            {
                Real tmp2 = fac*xt1*xt2*Real(0.25);

                IntVect iv = offset - offset_t1 - offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += tmp2;

                iv = offset + offset_t1 + offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) += tmp2;

                iv = offset - offset_t1 + offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) -= tmp2;

                iv = offset + offset_t1 - offset_t2;
                used(iv[0],iv[1],iv[2]) = true;
                tmp(iv[0],iv[1],iv[2]) -= tmp2;
            }
        }
    }

    auto const* ptmp = tmp.begin();
    auto const* pused = used.begin();
    for (int m = 0; m < 27; ++m) {
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
                       Array4<int const> const& offset_bz,
                       Real const* bx, Real const* by, Real const* bz,
                       Array4<int const> const& fine_mask,
                       IntVect const& rr)
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
        civ[fromoff][2] = k;
        nentries[fromoff] = nentries_to(i,j,k);
        int const tooff = offset_to(i,j,k);
        entry_offset[fromoff] = tooff;

        // Fist, we need to figure out how many corner coarse cells are
        // involved. The coarse cell entries must be stored ahead of fine
        // cell entries because that's how we sorted the entried when
        // building the hypre graph.
        bool corner[3] = {false, false, false};
        if ((fine_mask(i-1,j,k) || fine_mask(i+1,j,k)) &&
            (! fine_mask(i,j-1,k-1)) &&
            (! fine_mask(i,j+1,k-1)) &&
            (! fine_mask(i,j-1,k+1)) &&
            (! fine_mask(i,j+1,k+1)))
        {
            corner[0] = true;
        }
        if((fine_mask(i,j-1,k) || fine_mask(i,j+1,k)) &&
            (! fine_mask(i-1,j,k-1)) &&
            (! fine_mask(i+1,j,k-1)) &&
            (! fine_mask(i-1,j,k+1)) &&
            (! fine_mask(i+1,j,k+1)))
        {
            corner[1] = true;
        }
        if((fine_mask(i,j,k-1) || fine_mask(i,j,k+1)) &&
            (! fine_mask(i-1,j-1,k)) &&
            (! fine_mask(i+1,j-1,k)) &&
            (! fine_mask(i-1,j+1,k)) &&
            (! fine_mask(i+1,j+1,k)))
        {
            corner[2] = true;
        }
        int nentries_c = 4 * (int(corner[0]) + int(corner[1]) + int(corner[2]));
        int foff = tooff + nentries_c;

        // We must iterate the faces in the lexicographical order of fine
        // neighbor cells, because that's the order when non-stencil entries
        // were added to Hypre's graph. Also note that a coarse cell will
        // not have entries to fine cells at both ends of a direction. Thus
        // we do not have to worry about the order between fine cells at the
        // small and big ends of the same direction.

        if (fine_mask(i,j,k-1)) {
            stencil(i,j,k)[0] += stencil(i,j,k)[5];
            stencil(i,j,k)[5] = Real(0.0);
            // Reflux: sb/h^2*bz(i,j,k)*(phi(i,j,k)-phi(i,j,k-1)) is replaced by
            // sb/h*sum_{fine_faces}(dphi/dz*bz)/n_fine_faces
            Real dzf = dx[2] / Real(rr[2]);
            Real dzcinv = Real(1.0) / dx[2];
            Real dzfinv = Real(1.0) / dzf;
            Real cc[3];
            Real zz[3] = {dx[2]*Real(-0.5), dzf*Real(0.5), dzf*Real(1.5)};
            poly_interp_coeff<3>(dzf*Real(-0.5), zz, cc);
            for (int iry = 0; iry < rr[1]; ++iry) {
            for (int irx = 0; irx < rr[0]; ++irx) {
                Real bzm = bz ? bz[offset_bz(i,j,k)+irx+iry*rr[0]] : Real(1.0);
                Real fac = sb*dzcinv*dzfinv*bzm/Real(rr[0]*rr[1]);
                // int ii = i*rr[0] + irx
                // int jj = j*rr[1] + iry
                // int kk = k*rr[2]
                // dphi/dz = (phi_interp - phi_fine(kk-1))/dz_fine
                //         = (phi_coarse*cc[0] + phi_fine(kk-1)*(cc[1]-1)
                //                             + phi_fine(kk-2)* cc[2]) / dz_fine
                // So the entry for fine cell (kk-1) is (cc[1]-1)*fac
                //                  fine cell (kk-2) is  cc[2]   *fac
                entry_values[foff+irx+iry*rr[0]            ] += fac* cc[2];
                entry_values[foff+irx+iry*rr[0]+rr[0]*rr[1]] += fac*(cc[1]-Real(1.0));

                // The stencil and non-stencil coarse cells need updates too.
                Real x = Real(-0.5) + (irx+Real(0.5))/Real(rr[0]);
                Real y = Real(-0.5) + (iry+Real(0.5))/Real(rr[1]);
                Real fac0 = fac*cc[0];
                Real s0 = Real(1.0);

                if (!fine_mask(i-1,j,k) && !fine_mask(i+1,j,k)) {
                    s0 -= x*x;
                    stencil(i,j,k)[1] += fac0*Real(0.5)*x*(x-Real(1.0));
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x*(x+Real(1.0));
                } else if (!fine_mask(i-1,j,k)) {
                    s0 += Real(0.5)*x;
                    stencil(i,j,k)[1] += fac0*Real(-0.5)*x;
                } else {
                    s0 += Real(-0.5)*x;
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x;
                }

                if (!fine_mask(i,j-1,k) && !fine_mask(i,j+1,k)) {
                    s0 -= y*y;
                    stencil(i,j,k)[3] += fac0*Real(0.5)*y*(y-Real(1.0));
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y*(y+Real(1.0));
                } else if (!fine_mask(i,j-1,k)) {
                    s0 += Real(0.5)*y;
                    stencil(i,j,k)[3] += fac0*Real(-0.5)*y;
                } else {
                    s0 += Real(-0.5)*y;
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y;
                }

                stencil(i,j,k)[0] += fac0*s0;

                if (corner[2]) {
                    int offset = tooff + (corner[0] ? 2 : 0) + (corner[1] ? 2 : 0);
                    entry_values[offset++] += fac0*Real( 0.25)*x*y;
                    entry_values[offset++] += fac0*Real(-0.25)*x*y;
                    entry_values[offset++] += fac0*Real(-0.25)*x*y;
                    entry_values[offset  ] += fac0*Real( 0.25)*x*y;
                }
            }}
            foff += 2*rr[0]*rr[1];
        }

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
            for (int irz = 0; irz < rr[2]; ++irz) {
            for (int irx = 0; irx < rr[0]; ++irx) {
                Real bym = by ? by[offset_by(i,j,k)+irx+irz*rr[0]] : Real(1.0);
                Real fac = sb*dycinv*dyfinv*bym/Real(rr[0]*rr[2]);
                // int ii = i*rr[0] + irx
                // int jj = j*rr[1]
                // int kk = k*rr[2] + irz
                // dphi/dy = (phi_interp - phi_fine(jj-1)) / dy_fine
                //         = (phi_coarse*cc[0] + phi_fine(jj-1)*(cc[1]-1)
                //                             + phi_fine(jj-2)* cc[2]) / dy_fine
                // So the entry for fine cell (jj-1) is (cc[1]-1)*fac
                //                  fine cell (jj-2) is  cc[2]   *fac
                entry_values[foff+irx      +irz*rr[0]*2] += fac* cc[2];
                entry_values[foff+irx+rr[0]+irz*rr[0]*2] += fac*(cc[1]-Real(1.0));

                // The stencil and non-stencil coarse cells need updates too.
                Real x = Real(-0.5) + (irx+Real(0.5))/Real(rr[0]);
                Real z = Real(-0.5) + (irz+Real(0.5))/Real(rr[2]);
                Real fac0 = fac*cc[0];
                Real s0 = Real(1.0);

                if (!fine_mask(i-1,j,k) && !fine_mask(i+1,j,k)) {
                    s0 -= x*x;
                    stencil(i,j,k)[1] += fac0*Real(0.5)*x*(x-Real(1.0));
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x*(x+Real(1.0));
                } else if (!fine_mask(i-1,j,k)) {
                    s0 += Real(0.5)*x;
                    stencil(i,j,k)[1] += fac0*Real(-0.5)*x;
                } else {
                    s0 += Real(-0.5)*x;
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x;
                }

                if (!fine_mask(i,j,k-1) && !fine_mask(i,j,k+1)) {
                    s0 -= z*z;
                    stencil(i,j,k)[5] += fac0*Real(0.5)*z*(z-Real(1.0));
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z*(z+Real(1.0));
                } else if (!fine_mask(i,j,k-1)) {
                    s0 += Real(0.5)*z;
                    stencil(i,j,k)[5] += fac0*Real(-0.5)*z;
                } else {
                    s0 += Real(-0.5)*z;
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z;
                }

                stencil(i,j,k)[0] += fac0*s0;

                if (corner[1]) {
                    int offset = tooff + (corner[0] ? 1 : 0);
                    entry_values[offset++] += fac0*Real( 0.25)*x*z;
                    entry_values[offset++] += fac0*Real(-0.25)*x*z;
                    if (corner[0]) { offset += 2; }
                    if (corner[2]) { offset += 4; }
                    entry_values[offset++] += fac0*Real(-0.25)*x*z;
                    entry_values[offset  ] += fac0*Real( 0.25)*x*z;
                }
            }}
            foff += 2*rr[0]*rr[2];
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
            for (int irz = 0; irz < rr[2]; ++irz) {
            for (int iry = 0; iry < rr[1]; ++iry) {
                Real bxm = bx ? bx[offset_bx(i,j,k)+iry+irz*rr[1]] : Real(1.0);
                Real fac = sb*dxcinv*dxfinv*bxm/Real(rr[1]*rr[2]);
                // int ii = i*rr[0]
                // int jj = j*rr[1] + iry
                // int kk = k*rr[2]
                // dphi/dx = (phi_interp - phi_fine(ii-1)) / dx_fine
                //         = (phi_coarse*cc[0] + phi_fine(ii-1)*(cc[1]-1)
                //                             + phi_fine(ii-2)* cc[2]) / dx_fine
                // So the entry for fine cell(ii-1) is (cc[1]-1)*fac
                //                  fine cell(ii-2) is  cc[2]   *fac
                entry_values[foff++] = fac* cc[2];
                entry_values[foff++] = fac*(cc[1] - Real(1.0));

                // The stencil and non-stencil coarse cells need updates too.
                Real y = Real(-0.5) + (iry+Real(0.5))/Real(rr[1]);
                Real z = Real(-0.5) + (irz+Real(0.5))/Real(rr[2]);
                Real fac0 = fac*cc[0];
                Real s0 = Real(1.0);

                if (!fine_mask(i,j-1,k) && !fine_mask(i,j+1,k)) {
                    s0 -= y*y;
                    stencil(i,j,k)[3] += fac0*Real(0.5)*y*(y-Real(1.0));
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y*(y+Real(1.0));
                } else if (!fine_mask(i,j-1,k)) {
                    s0 += Real(0.5)*y;
                    stencil(i,j,k)[3] += fac0*Real(-0.5)*y;
                } else {
                    s0 += Real(-0.5)*y;
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y;
                }

                if (!fine_mask(i,j,k-1) && !fine_mask(i,j,k+1)) {
                    s0 -= z*z;
                    stencil(i,j,k)[5] += fac0*Real(0.5)*z*(z-Real(1.0));
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z*(z+Real(1.0));
                } else if (!fine_mask(i,j,k-1)) {
                    s0 += Real(0.5)*z;
                    stencil(i,j,k)[5] += fac0*Real(-0.5)*z;
                } else {
                    s0 += Real(-0.5)*z;
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z;
                }

                stencil(i,j,k)[0] += fac0*s0;

                if (corner[0]) {
                    int offset = tooff;
                    entry_values[offset++] += fac0*Real( 0.25)*y*z;
                    if (corner[1]) { offset += 2; }
                    entry_values[offset++] += fac0*Real(-0.25)*y*z;
                    if (corner[2]) { offset += 4; }
                    entry_values[offset++] += fac0*Real(-0.25)*y*z;
                    if (corner[1]) { offset += 2; }
                    entry_values[offset  ] += fac0*Real( 0.25)*y*z;
                }
            }}
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
            for (int irz = 0; irz < rr[2]; ++irz) {
            for (int iry = 0; iry < rr[1]; ++iry) {
                Real bxp = bx ? bx[offset_bx(i+1,j,k)+iry+irz*rr[1]] : Real(1.0);
                Real fac = sb*dxcinv*dxfinv*bxp/Real(rr[1]*rr[2]);
                // int ii = i*rr[0] + (rr[0]-1)
                // int jj = j*rr[1] + iry
                // -dphi/dx = (phi_interp - phi_fine(ii+1)) / dx_fine
                //          = (phi_coarse*cc[0] + phi_fine(ii+1)*(cc[1]-1)
                //                              + phi_fine(ii+2)* cc[2]) / dx_fine
                // So the entry for fine cell(ii+1) is (cc[1]-1)*fac
                //                  fine cell(ii+2) is  cc[2]   *fac
                entry_values[foff++] = fac*(cc[1] - Real(1.0));
                entry_values[foff++] = fac* cc[2];

                // The stencil and non-stencil coarse cells need updates too.
                Real y = Real(-0.5) + (iry+Real(0.5))/Real(rr[1]);
                Real z = Real(-0.5) + (irz+Real(0.5))/Real(rr[2]);
                Real fac0 = fac*cc[0];
                Real s0 = Real(1.0);

                if (!fine_mask(i,j-1,k) && !fine_mask(i,j+1,k)) {
                    s0 -= y*y;
                    stencil(i,j,k)[3] += fac0*Real(0.5)*y*(y-Real(1.0));
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y*(y+Real(1.0));
                } else if (!fine_mask(i,j-1,k)) {
                    s0 += Real(0.5)*y;
                    stencil(i,j,k)[3] += fac0*Real(-0.5)*y;
                } else {
                    s0 += Real(-0.5)*y;
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y;
                }

                if (!fine_mask(i,j,k-1) && !fine_mask(i,j,k+1)) {
                    s0 -= z*z;
                    stencil(i,j,k)[5] += fac0*Real(0.5)*z*(z-Real(1.0));
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z*(z+Real(1.0));
                } else if (!fine_mask(i,j,k-1)) {
                    s0 += Real(0.5)*z;
                    stencil(i,j,k)[5] += fac0*Real(-0.5)*z;
                } else {
                    s0 += Real(-0.5)*z;
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z;
                }

                stencil(i,j,k)[0] += fac0*s0;

                if (corner[0]) {
                    int offset = tooff;
                    entry_values[offset++] += fac0*Real( 0.25)*y*z;
                    if (corner[1]) { offset += 2; }
                    entry_values[offset++] += fac0*Real(-0.25)*y*z;
                    if (corner[2]) { offset += 4; }
                    entry_values[offset++] += fac0*Real(-0.25)*y*z;
                    if (corner[1]) { offset += 2; }
                    entry_values[offset  ] += fac0*Real( 0.25)*y*z;
                }
            }}
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
            for (int irz = 0; irz < rr[2]; ++irz) {
            for (int irx = 0; irx < rr[0]; ++irx) {
                Real byp = by ? by[offset_by(i,j+1,k)+irx+irz*rr[0]] : Real(1.0);
                Real fac = sb*dycinv*dyfinv*byp/Real(rr[0]*rr[2]);
                // int ii = i*rr[0] + irx
                // int jj = j*rr[1] + (rr[1]-1)
                // int kk = k*rr[2] + irz;
                // -dphi/dy = (phi_interp - phi_fine(jj+1)) / dy_fine
                //          = (phi_coarse*cc[0] + phi_fine(jj+1)*(cc[1]-1)
                //                              + phi_fine(jj+2)* cc[2]) / dy_fine
                // So the entry for fine cell (jj+1) is (cc[1]-1)*fac
                //                  fine cell (jj+2) is  cc[2]   *fac
                entry_values[foff+irx      +irz*rr[0]*2] += fac*(cc[1]-Real(1.0));
                entry_values[foff+irx+rr[0]+irz*rr[0]*2] += fac* cc[2];

                // The stencil and non-stencil coarse cells need updates too.
                Real x = Real(-0.5) + (irx+Real(0.5))/Real(rr[0]);
                Real z = Real(-0.5) + (irz+Real(0.5))/Real(rr[2]);
                Real fac0 = fac*cc[0];
                Real s0 = Real(1.0);

                if (!fine_mask(i-1,j,k) && !fine_mask(i+1,j,k)) {
                    s0 -= x*x;
                    stencil(i,j,k)[1] += fac0*Real(0.5)*x*(x-Real(1.0));
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x*(x+Real(1.0));
                } else if (!fine_mask(i-1,j,k)) {
                    s0 += Real(0.5)*x;
                    stencil(i,j,k)[1] += fac0*Real(-0.5)*x;
                } else {
                    s0 += Real(-0.5)*x;
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x;
                }

                if (!fine_mask(i,j,k-1) && !fine_mask(i,j,k+1)) {
                    s0 -= z*z;
                    stencil(i,j,k)[5] += fac0*Real(0.5)*z*(z-Real(1.0));
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z*(z+Real(1.0));
                } else if (!fine_mask(i,j,k-1)) {
                    s0 += Real(0.5)*z;
                    stencil(i,j,k)[5] += fac0*Real(-0.5)*z;
                } else {
                    s0 += Real(-0.5)*z;
                    stencil(i,j,k)[6] += fac0*Real(0.5)*z;
                }

                stencil(i,j,k)[0] += fac0*s0;

                if (corner[1]) {
                    int offset = tooff + (corner[0] ? 1 : 0);
                    entry_values[offset++] += fac0*Real( 0.25)*x*z;
                    entry_values[offset++] += fac0*Real(-0.25)*x*z;
                    if (corner[0]) { offset += 2; }
                    if (corner[2]) { offset += 4; }
                    entry_values[offset++] += fac0*Real(-0.25)*x*z;
                    entry_values[offset  ] += fac0*Real( 0.25)*x*z;
                }
            }}
            foff += 2*rr[0]*rr[2];
        }

        if (fine_mask(i,j,k+1)) {
            stencil(i,j,k)[0] += stencil(i,j,k)[6];
            stencil(i,j,k)[6] = Real(0.0);
            // Reflux: sb/h^2*bz(i,j,k+1)*(phi(i,j,k)-phi(i,j,k+1)) is replaced by
            // sb/h*sum_{fine_faces}(-dphi/dz*bz)/n_fine_faces
            Real dzf = dx[2] / Real(rr[2]);
            Real dzcinv = Real(1.0) / dx[2];
            Real dzfinv = Real(1.0) / dzf;
            Real cc[3];
            Real zz[3] = {dx[2]*Real(-0.5), dzf*Real(0.5), dzf*Real(1.5)};
            poly_interp_coeff<3>(dzf*Real(-0.5), zz, cc);
            for (int iry = 0; iry < rr[1]; ++iry) {
            for (int irx = 0; irx < rr[0]; ++irx) {
                Real bzp = bz ? bz[offset_bz(i,j,k+1)+irx+iry*rr[0]] : Real(1.0);
                Real fac = sb*dzcinv*dzfinv*bzp/Real(rr[0]*rr[1]);
                // int ii = i*rr[0] + irx
                // int jj = j*rr[1] + iry
                // int kk = k*rr[2] + (rr[2]-1)
                // -dphi/dz = (phi_interp - phi_fine(kk+1))/dz_fine
                //          = (phi_coarse*cc[0] + phi_fine(kk+1)*(cc[1]-1)
                //                              + phi_fine(kk+2)* cc[2]) / dz_fine
                // So the entry for fine cell (kk+1) is (cc[1]-1)*fac
                //                  fine cell (kk+2) is  cc[2]   *fac
                entry_values[foff+irx+iry*rr[0]            ] += fac*(cc[1]-Real(1.0));
                entry_values[foff+irx+iry*rr[0]+rr[0]*rr[1]] += fac* cc[2];

                // The stencil and non-stencil coarse cells need updates too.
                Real x = Real(-0.5) + (irx+Real(0.5))/Real(rr[0]);
                Real y = Real(-0.5) + (iry+Real(0.5))/Real(rr[1]);
                Real fac0 = fac*cc[0];
                Real s0 = Real(1.0);

                if (!fine_mask(i-1,j,k) && !fine_mask(i+1,j,k)) {
                    s0 -= x*x;
                    stencil(i,j,k)[1] += fac0*Real(0.5)*x*(x-Real(1.0));
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x*(x+Real(1.0));
                } else if (!fine_mask(i-1,j,k)) {
                    s0 += Real(0.5)*x;
                    stencil(i,j,k)[1] += fac0*Real(-0.5)*x;
                } else {
                    s0 += Real(-0.5)*x;
                    stencil(i,j,k)[2] += fac0*Real(0.5)*x;
                }

                if (!fine_mask(i,j-1,k) && !fine_mask(i,j+1,k)) {
                    s0 -= y*y;
                    stencil(i,j,k)[3] += fac0*Real(0.5)*y*(y-Real(1.0));
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y*(y+Real(1.0));
                } else if (!fine_mask(i,j-1,k)) {
                    s0 += Real(0.5)*y;
                    stencil(i,j,k)[3] += fac0*Real(-0.5)*y;
                } else {
                    s0 += Real(-0.5)*y;
                    stencil(i,j,k)[4] += fac0*Real(0.5)*y;
                }

                stencil(i,j,k)[0] += fac0*s0;

                if (corner[2]) {
                    int offset = tooff + (corner[0] ? 2 : 0) + (corner[1] ? 2 : 0);
                    entry_values[offset++] += fac0*Real( 0.25)*x*y;
                    entry_values[offset++] += fac0*Real(-0.25)*x*y;
                    entry_values[offset++] += fac0*Real(-0.25)*x*y;
                    entry_values[offset  ] += fac0*Real( 0.25)*x*y;
                }
            }}
            // no need to foff += 2*rr[0]*rr[1];
        }
    }
}

}

#endif
