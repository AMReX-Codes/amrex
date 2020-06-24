#include <AMReX_Habec_3D_K.H>

namespace amrex {

void amrex_hpacoef (Box const& box, 
                    BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                    const FArrayBox& a,
                    Real& sa)
{
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real const> const& a_arr = a.const_array();
    if (sa == 0.0){
        AMREX_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = 0.0;
        });
    }else{
        AMREX_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = sa * a_arr(i,j,k);
        });
    }
}

void amrex_hpbcoef (Box const& box,
                    BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                    const FArrayBox& b,
                    Real& sb,
                    const Real* dx,
                    int& idim)
{
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real const> const& b_arr = b.const_array();
    Real fac = sb / (dx[idim]*dx[idim]);

    if (idim == 0) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + fac * (b_arr(i,j,k) + b_arr(i+1,j,k));
            mat_arr(i,j,k)[1] = - fac * b_arr(i,j,k);
            mat_arr(i,j,k)[2] = - fac * b_arr(i+1,j,k);
        });
    }else if (idim == 1) {
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + fac * (b_arr(i,j,k) + b_arr(i,j+1,k));
            mat_arr(i,j,k)[3] = - fac * b_arr(i,j,k);
            mat_arr(i,j,k)[4] = - fac * b_arr(i,j+1,k);
        });
    }else{
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + fac * (b_arr(i,j,k) + b_arr(i,j,k+1));
            mat_arr(i,j,k)[5] = - fac * b_arr(i,j,k);
            mat_arr(i,j,k)[6] = - fac * b_arr(i,j,k+1);
        });
    }
}

void amrex_hpmat (Box const& box,
                  BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                  const FArrayBox& b,
                  const Mask& msk,
                  Real& sb,
                  const Real* dx,
                  int& cdir,
                  const int& bct,
                  const Real& bcl,
                  const int& bho) 
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real const> const& b_arr = b.const_array();
    Array4<int const> const& msk_arr = msk.const_array();

    Real h = 0.0;
    if (cdir == 0 || cdir == 3) {
        h = dx[0];
    }else if (cdir == 1 || cdir == 4) {
        h = dx[1];
    }else{
        h = dx[2];
    }
    Real fac = sb/(h*h);

    Real h2 = 0.0;
    Real h3 = 0.0;
    Real bf1 = 0.0;
    Real bf2 = 0.0;
    if (bct == AMREX_LO_DIRICHLET) {
        h2 = 0.5*h;
        if (bho >= 1) {
            h3 = 3.0*h2;
            bf1 = fac*((h3 - bcl)/(bcl + h2) - 1.0);
            bf2 = fac*(bcl - h2)/(bcl + h3);
        }else{
            bf1 = fac*( h/(bcl + h2) - 1.0);
            bf2 = 0.0;
        }
    }else if (bct == AMREX_LO_NEUMANN) {
        bf1 = -fac;
        bf2 = 0.0;
    }else{
        Abort("hpmat: unsupported boundary type");
    }

    if (cdir == 0){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            if (msk_arr(lo.x-1,j,k) > 0 && i==lo.x){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j,k);
                mat_arr(i,j,k)[1] = 0.0;
                mat_arr(i,j,k)[2] = mat_arr(i,j,k)[2] + bf2*b_arr(i,j,k);
            }
        });
    }else if (cdir == 3){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            if (msk_arr(hi.x+1,j,k) > 0 && i==hi.x){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i+1,j,k);
                mat_arr(i,j,k)[2] = 0.0;
                mat_arr(i,j,k)[1] = mat_arr(i,j,k)[1] + bf2*b_arr(i+1,j,k);
            }
        });
   }else if (cdir == 1){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            if (msk_arr(i,lo.y-1,k) > 0 && j==lo.y){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j,k);
                mat_arr(i,j,k)[3] = 0.0;
                mat_arr(i,j,k)[4] = mat_arr(i,j,k)[4] + bf2*b_arr(i,j,k);
            }
        });
   }else if (cdir == 4){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            if (msk_arr(i,hi.y+1,k) > 0 && j==hi.y){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j+1,k);
                mat_arr(i,j,k)[4] = 0.0;
                mat_arr(i,j,k)[3] = mat_arr(i,j,k)[3] + bf2*b_arr(i,j+1,k);
            }
        });
   }else if (cdir == 2){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            if (msk_arr(i,j,lo.z-1) > 0 && k==lo.z){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j,k);
                mat_arr(i,j,k)[5] = 0.0;
                mat_arr(i,j,k)[6] = mat_arr(i,j,k)[6] + bf2*b_arr(i,j,k);
            }
        });   
   }else if (cdir == 5){
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
        {
            if (msk_arr(i,j,hi.z+1) > 0 && k==hi.z){
                mat_arr(i,j,k)[0] = mat_arr(i,j,k)[0] + bf1*b_arr(i,j,k+1);
                mat_arr(i,j,k)[6] = 0.0;
                mat_arr(i,j,k)[5] = mat_arr(i,j,k)[5] + bf2*b_arr(i,j,k+1);
            }
        });
   }else{
        Abort("hpmat: impossible face orientation");
   }
}

void amrex_hpdiag (Box const& box,
                   BaseFab<GpuArray<Real,2*AMREX_SPACEDIM + 1>>& mat,
                   FArrayBox& diag) 
{
    Array4<GpuArray<Real,2*AMREX_SPACEDIM + 1>> const& mat_arr = mat.array();
    Array4<Real> const& diag_arr = diag.array();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {  
        diag_arr(i,j,k) = 1.0/mat_arr(i,j,k)[0]; 
        for (int ii=0; ii<2*AMREX_SPACEDIM + 1; ii++){
            mat_arr(i,j,k)[ii] = mat_arr(i,j,k)[ii]*diag_arr(i,j,k);
        }
    });
}

void amrex_hpijmatrix (Box const& box,
                       const HYPRE_Int& nrows, 
                       HYPRE_Int* ncols,
                       HYPRE_Int* rows, 
                       HYPRE_Int* colsg,
                       Real* matg, 
                       BaseFab<int>& cell_id,
                       HYPRE_Int& offset, FArrayBox& diaginv,
                       FArrayBox& acoefs, FArrayBox& bcoefsx,
                       FArrayBox& bcoefsy, FArrayBox& bcoefsz,
                       Real& sa, Real& sb, const Real* dx, 
                       GpuArray<int,AMREX_SPACEDIM*2>& bct,
                       GpuArray<Real,AMREX_SPACEDIM*2> bcl,
                       const int& bho) 
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    GpuArray<Real,AMREX_SPACEDIM> fac;
    GpuArray<Real,AMREX_SPACEDIM*2> bf1;
    GpuArray<Real,AMREX_SPACEDIM*2> bf2;

    for (int i=0; i<AMREX_SPACEDIM; i++){
        fac[i] = sb/(dx[i]*dx[i]);
    }
   
    int idim = 0;
    for (int cdir=0; cdir<=AMREX_SPACEDIM*2-1; cdir++){ 
       if (cdir == 0 || cdir == 3){
          idim = 0;
       }else if (cdir == 1 || cdir == 4){
          idim = 1;
       }else{
          idim = 2;
       }

       Real h = dx[idim];

       if (bct[cdir] == AMREX_LO_DIRICHLET) {
          Real h2 = 0.5*h;
          if (bho>=1) {
             Real h3 = 3.0*h2;
             bf1[cdir] = fac[idim] * ((h3 - bcl[cdir]) / (bcl[cdir] + h2) - 1.0);
             bf2[cdir] = fac[idim] * (bcl[cdir] - h2) / (bcl[cdir] + h3);
          }else{
             bf1[cdir] = fac[idim] * ( h / (bcl[cdir] + h2) - 1.0);
             bf2[cdir] = 0.0;
          }
       }else if (bct[cdir] == AMREX_LO_NEUMANN) {
          bf1[cdir] = -fac[idim];
          bf2[cdir] = 0.0;
       }
    }

    Gpu::DeviceVector<int> iter(2,0);

    Array4<int> const& cell_id_arr = cell_id.array();
    Gpu::DeviceVector<int> cols_tmp(AMREX_SPACEDIM*2+1,0.0);
    Gpu::DeviceVector<Real> mat_tmp(AMREX_SPACEDIM*2+1,0.0);
    Array4<Real const> const& a_arr = acoefs.const_array();
    Array4<Real const> const& bx_arr = bcoefsx.const_array();
    Array4<Real const> const& by_arr = bcoefsy.const_array();
    Array4<Real const> const& bz_arr = bcoefsz.const_array();
    Array4<Real> const& diag_arr = diaginv.array();

    auto cols_tmpPtr = cols_tmp.dataPtr();
    auto mat_tmpPtr = mat_tmp.dataPtr();
    auto iterPtr = iter.dataPtr();

    AMREX_PARALLEL_FOR_3D (box, i, j, k,
    {
        int irow = iterPtr[0];
        rows[irow]  = cell_id_arr(i,j,k);
        ncols[irow] = 0;

        cols_tmpPtr[0] = cell_id_arr(i,j,k);
        mat_tmpPtr[0]  = sa*a_arr(i,j,k) + fac[0]*(bx_arr(i,j,k)+bx_arr(i+1,j,k)) 
                                      + fac[1]*(by_arr(i,j,k)+by_arr(i,j+1,k)) 
                                      + fac[2]*(bz_arr(i,j,k)+bz_arr(i,j,k+1));

        cols_tmpPtr[1] = cell_id_arr(i-1,j,k);
        mat_tmpPtr[1]  = -fac[0]*bx_arr(i,j,k);

        cols_tmpPtr[2] = cell_id_arr(i+1,j,k);
        mat_tmpPtr[2]  = -fac[0]*bx_arr(i+1,j,k);

        cols_tmpPtr[3] = cell_id_arr(i,j-1,k);
        mat_tmpPtr[3]  = -fac[1]*by_arr(i,j,k);

        cols_tmpPtr[4] = cell_id_arr(i,j+1,k);
        mat_tmpPtr[4]  = -fac[1]*by_arr(i,j+1,k);

        cols_tmpPtr[5] = cell_id_arr(i,j,k-1);
        mat_tmpPtr[5]  = -fac[2]*bz_arr(i,j,k);

        cols_tmpPtr[6] = cell_id_arr(i,j,k+1);
        mat_tmpPtr[6]  = -fac[2]*bz_arr(i,j,k+1);

        if (i == lo.x && cell_id_arr(i-1,j,k)<0) {
           mat_tmpPtr[0] += bf1[0]*bx_arr(i,j,k);
           mat_tmpPtr[2] += bf2[0]*bx_arr(i,j,k);
        }

        if (i == hi.x && cell_id_arr(i+1,j,k)<0) {
           mat_tmpPtr[0] += bf1[3]*bx_arr(i+1,j,k);
           mat_tmpPtr[1] += bf2[3]*bx_arr(i+1,j,k);
        }

        if (j == lo.y && cell_id_arr(i,j-1,k)<0) {
           mat_tmpPtr[0] += bf1[1]*by_arr(i,j,k);
           mat_tmpPtr[4] += bf2[1]*by_arr(i,j,k);
        }

        if (j == hi.y && cell_id_arr(i,j+1,k)<0) {
           mat_tmpPtr[0] += bf1[4]*by_arr(i,j+1,k);
           mat_tmpPtr[3] += bf2[4]*by_arr(i,j+1,k);
        }

        if (k == lo.z && cell_id_arr(i,j,k-1)<0) {
           mat_tmpPtr[0] += bf1[2]*bz_arr(i,j,k);
           mat_tmpPtr[6] += bf2[2]*bz_arr(i,j,k);
        }

        if (k == hi.z && cell_id_arr(i,j,k+1)<0) {
           mat_tmpPtr[0] += bf1[5]*bz_arr(i,j,k+1);
           mat_tmpPtr[5] += bf2[5]*bz_arr(i,j,k+1);
        }

        diag_arr(i,j,k) = 1.0/mat_tmpPtr[0];

        for (int ic=0; ic<=6; ic++){
            if (cols_tmpPtr[ic] >= 0){
                int imat = iterPtr[1];
                ncols[irow]  = ncols[irow] +1;
                colsg[imat]  = cols_tmpPtr[ic];
                matg[imat]   = mat_tmpPtr[ic]*diag_arr(i,j,k);
                Gpu::Atomic::Add(&iterPtr[1], 1);
            }
        }
        Gpu::Atomic::Add(&iterPtr[0], 1);
    });
}

#ifdef AMREX_USE_EB

void amrex_hpeb_fill_cellid (Box const& box,
                             int& nrows,
                             BaseFab<HYPRE_Int>& cell_id, 
                             const EBCellFlagFab& flag) 
{
    Array4<int> const& cell_id_arr = cell_id.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();

    nrows = 0;
    Gpu::DeviceScalar<int> nrows_gpu(nrows);
    int* nrowsg = nrows_gpu.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        if (!flag_arr(i,j,k).isCovered()){
            cell_id_arr(i,j,k) = *nrowsg;
            Gpu::Atomic::Add(nrowsg, 1);
        }
    });
    nrows = nrows_gpu.dataValue();
}

void amrex_hpeb_copy_from_vec (Box const& box,
                               FArrayBox& a,
                               Real* v,
                               const EBCellFlagFab& flag) 
{
    Array4<Real> const& a_arr = a.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();

    int nrows = 0;
    Gpu::DeviceScalar<int> nrows_gpu(nrows);
    int* nrowsg = nrows_gpu.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        if (!flag_arr(i,j,k).isCovered()){
            a_arr(i,j,k) = v[*nrowsg];
            Gpu::Atomic::Add(nrowsg, 1);
        }
    });
}

void amrex_hpeb_copy_to_vec (Box const& box,
                             FArrayBox& a,
                             Real* v,
                             const EBCellFlagFab& flag) 
{
    Array4<Real> const& a_arr = a.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();

    int nrows = 0;
    Gpu::DeviceScalar<int> nrows_gpu(nrows);
    int* nrowsg = nrows_gpu.dataPtr();

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (box, i, j, k,
    {
        if (!flag_arr(i,j,k).isCovered()){
            v[*nrowsg] = a_arr(i,j,k);
            Gpu::Atomic::Add(nrowsg, 1);
        }
    });
}

void amrex_hpeb_ijmatrix (Box const& box,
                          const HYPRE_Int& nrows,
                          HYPRE_Int* ncols,
                          HYPRE_Int* rows,
                          HYPRE_Int* colsg,
                          Real* matg,
                          BaseFab<int>& cell_id,
                          HYPRE_Int& offset, FArrayBox& diaginv,
                          FArrayBox& acoefs, FArrayBox& bcoefsx,
                          FArrayBox& bcoefsy, FArrayBox& bcoefsz,
                          const EBCellFlagFab& flag,
                          const FArrayBox& vfrc,
                          const FArrayBox& apx, const FArrayBox& apy,
                          const FArrayBox& apz, const FArrayBox& fcx,
                          const FArrayBox& fcy, const FArrayBox& fcz,
                          const FArrayBox& ba, const FArrayBox& bcen,
                          const FArrayBox& beb,
                          const int& is_eb_dirichlet,
                          Real& sa, Real& sb, const Real* dx,
                          GpuArray<int,AMREX_SPACEDIM*2>& bct,
                          GpuArray<Real,AMREX_SPACEDIM*2> bcl,
                          const int& bho) 
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    GpuArray<Real,AMREX_SPACEDIM> fac;
    GpuArray<Real,AMREX_SPACEDIM*2> bf1;
    GpuArray<Real,AMREX_SPACEDIM*2> bf2;
    GpuArray<Real,AMREX_SPACEDIM*2> bflo;
 
    int is_dirichlet = (is_eb_dirichlet != 0) ? 1 : 0;

    for (int i=0; i<AMREX_SPACEDIM; i++){
        fac[i] = sb/(dx[i]*dx[i]);
    }
   
    int idim = 0;
    for (int cdir=0; cdir<2*AMREX_SPACEDIM; cdir++){ 
       if (cdir == 0 || cdir == 3){
          idim = 0;
       }else if (cdir == 1 || cdir == 4){
          idim = 1;
       }else{
          idim = 2;
       }

       Real h = dx[idim];

       if (bct[cdir] == AMREX_LO_DIRICHLET) {
          Real h2 = 0.5*h;
          bflo[cdir] = fac[idim] * ( h / (bcl[cdir] + h2) - 1.0); 
          if (bho>=1) {
             Real h3 = 3.0*h2;
             bf1[cdir] = fac[idim] * ((h3 - bcl[cdir]) / (bcl[cdir] + h2) - 1.0);
             bf2[cdir] = fac[idim] * (bcl[cdir] - h2) / (bcl[cdir] + h3);
          }else{
             bf1[cdir] = bflo[cdir];
             bf2[cdir] = 0.0;
          }
       }else if (bct[cdir] == AMREX_LO_NEUMANN) {
          bflo[cdir] = -fac[idim];
          bf1[cdir] = -fac[idim];
          bf2[cdir] = 0.0;
       }
    }

    Gpu::DeviceVector<int> iter(2,0);
    Array4<int> const& cell_id_arr = cell_id.array();
    Array4<Real const> const& a_arr = acoefs.const_array();
    Array4<Real const> const& bx_arr = bcoefsx.const_array();
    Array4<Real const> const& by_arr = bcoefsy.const_array();
    Array4<Real const> const& bz_arr = bcoefsz.const_array();
    Array4<Real> const& diag_arr = diaginv.array();
    Array4<const EBCellFlag> const& flag_arr = flag.array();
    Array4<Real const> const& apx_arr = apx.const_array();
    Array4<Real const> const& fcx_arr = fcx.const_array();
    Array4<Real const> const& apy_arr = apy.const_array();
    Array4<Real const> const& fcy_arr = fcy.const_array();
    Array4<Real const> const& apz_arr = apz.const_array();
    Array4<Real const> const& fcz_arr = fcz.const_array();
    Array4<Real const> const& bcen_arr = bcen.const_array();
    Array4<Real const> const& vfrc_arr = vfrc.const_array();
    Array4<Real const> const& ba_arr = ba.const_array();
    Array4<Real const> const& beb_arr = beb.const_array();

    auto iterPtr = iter.dataPtr();

    AMREX_HOST_DEVICE_FOR_3D (box, i, j, k,
    {
        hpeb_ij(i,j,k,nrows,ncols,rows,colsg,matg,fac,bf1,bf2,bflo,
                  lo,hi,sa,bho,is_dirichlet,
                  iterPtr,cell_id_arr,a_arr,bx_arr,
                  by_arr,bz_arr,diag_arr,flag_arr,apx_arr,fcx_arr,
                  apy_arr,fcy_arr,apz_arr,fcz_arr,bcen_arr,vfrc_arr,
                  ba_arr,beb_arr);
    });
}

#endif

}

