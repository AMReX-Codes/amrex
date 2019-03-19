#include <WarpX.H>
#include <BilinearFilter.H>
#include <Filter_Kernels.H>
#include <WarpX_f.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

void compute_stencil(Vector<Real> &stencil, int npass){
    Vector<Real> old_s(1+npass,0.);
    Vector<Real> new_s(1+npass,0.);

    old_s[0] = 1.;
    int jmax = 1;
    amrex::Real loc;
    // Convolve the filter with itself npass times
    for(int ipass=1; ipass<npass+1; ipass++){
        // element 0 has to be treated in its own way
        new_s[0] = 0.5 * old_s[0];
        if (1<jmax) new_s[0] += 0.5 * old_s[1];
        loc = 0.;
        // For each element j, apply the filter to 
        // old_s to get new_s[j]. loc stores the tmp 
        // filtered value.
        for(int j=1; j<jmax+1; j++){
            loc = 0.5 * old_s[j];
            loc += 0.25 * old_s[j-1];
            if (j<jmax) loc += 0.25 * old_s[j+1];
            new_s[j] = loc;
        }
        // copy new_s into old_s
        old_s = new_s;
        // extend the stencil length for next iteration
        jmax += 1;
    }
    // we use old_s here to make sure the stencil
    // is corrent even when npass = 0
    stencil = old_s;
}

void BilinearFilter::ComputeStencils(){
    Print()<<"npass_each_dir "<<npass_each_dir<<'\n';
    stencil_length_each_dir = npass_each_dir;
    stencil_length_each_dir += 1.;
#if (AMREX_SPACEDIM == 3)
    // npass_each_dir = npass_x npass_y npass_z
    stencil_x.resize( 1 + npass_each_dir[0] );
    stencil_y.resize( 1 + npass_each_dir[1] );
    stencil_z.resize( 1 + npass_each_dir[2] );
    compute_stencil(stencil_x, npass_each_dir[0]);
    compute_stencil(stencil_y, npass_each_dir[1]);
    compute_stencil(stencil_z, npass_each_dir[2]);
#elif (AMREX_SPACEDIM == 2)
    // npass_each_dir = npass_x npass_z
    stencil_x.resize( 1 + npass_each_dir[0] );
    stencil_z.resize( 1 + npass_each_dir[1] );
    compute_stencil(stencil_x, npass_each_dir[0]);
    compute_stencil(stencil_z, npass_each_dir[1]);
#endif
}


void
BilinearFilter::ApplyStencil (MultiFab& dstmf, const MultiFab& srcmf, int scomp, int dcomp, int ncomp)
{
    ncomp = std::min(ncomp, srcmf.nComp());
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Print()<<"stencil_x: "; for(int i=0; i<stencil_length_each_dir[0];i++){Print()<<stencil_x[i]<<" ";} Print()<<'\n';
        Print()<<"stencil_z: "; for(int i=0; i<stencil_length_each_dir[1];i++){Print()<<stencil_z[i]<<" ";} Print()<<'\n';
        FArrayBox tmpfab;
        for (MFIter mfi(dstmf,true); mfi.isValid(); ++mfi){
            const auto& srcfab = srcmf[mfi];
            auto& dstfab = dstmf[mfi];
            const Box& tbx = mfi.growntilebox();
            const Box& gbx = amrex::grow(tbx,1);
            tmpfab.resize(gbx,ncomp);
            tmpfab.setVal(0.0, gbx, 0, ncomp);
            const Box& ibx = gbx & srcfab.box();
            tmpfab.copy(srcfab, ibx, scomp, ibx, 0, ncomp);
            

            //#define BL_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), AMREX_ARLIM_ANYD((x).loVect()), AMREX_ARLIM_ANYD((x).hiVect())

            //WRPX_FILTER(BL_TO_FORTRAN_BOX(tbx),
            //            BL_TO_FORTRAN_ANYD(tmpfab),
            //            BL_TO_FORTRAN_N_ANYD(dstfab,dcomp),
            //            ncomp);
            filter_2d(tbx, gbx, tmpfab, dstfab, ncomp);
        }
    }
}

void BilinearFilter::filter_2d(const Box& tbx, const Box& gbx, FArrayBox &tmpfab, FArrayBox &dstfab, int ncomp)
{
    //    ncomp = std::min(ncomp, srcmf.nComp());
    Print()<<"ncomp "<<ncomp<<'\n';
    //const int* loVect = gbx.loVect();
    //const int* hiVect = gbx.hiVect();
    const int* loVector = gbx.loVect();
    const int* hiVector = gbx.hiVect();
    //const int* loVect = amrex::lbound(tbx);
    //const int* hiVect = amrex::hbound(tbx);
    Print()<<"loVector: "; for(int i=0; i<2;i++){Print()<<loVector[i]<<" ";} Print()<<'\n';
    Print()<<"hiVector: "; for(int i=0; i<2;i++){Print()<<hiVector[i]<<" ";} Print()<<'\n';
    Print()<<"stencil_length_each_dir: "; for(int i=0; i<2;i++){Print()<<stencil_length_each_dir[i]<<" ";} Print()<<'\n';
    Print()<<"stencil_x: "; for(int i=0; i<stencil_length_each_dir[0];i++){Print()<<stencil_x[i]<<" ";} Print()<<'\n';
    Print()<<"stencil_z: "; for(int i=0; i<stencil_length_each_dir[1];i++){Print()<<stencil_z[i]<<" ";} Print()<<'\n';
    Print()<<"above toto\n";
    Array4<Real> const& tmparr = tmpfab.array();
    Array4<Real> const& dstarr = dstfab.array();
    //    for(int c=0; c<ncomp; c++){
    for(int i=loVector[0]; i<hiVector[0]; i++){
    for(int j=loVector[1]; j<hiVector[1]; j++){
        for (int ix=-stencil_length_each_dir[0]+1; ix<stencil_length_each_dir[0]; ix++){
        for (int iz=-stencil_length_each_dir[1]+1; iz<stencil_length_each_dir[1]; iz++){
            //Print()<<"ix "<<ix<<" "<<stencil_x[abs(ix)]<<" "<<"iz "<<iz<<" "<<stencil_z[abs(iz)]<<"\n";
            // dstarr(i,j,0) = tmparr(i,j,0) + stencil_x[abs(ix)]*stencil_z[abs(iz)]*tmparr(i+ix,j+iz,0);
            dstarr(i,j,0) = stencil_x[abs(ix)]*stencil_z[abs(iz)]*tmparr(i+ix,j+iz,0);
            //Print()<<tmparr(i+ix,j+iz,c)<<' '<<dstarr(i,j,c)<<'\n';
        }
        }
    }
    }
    //}
    Print()<<"here --- 6\n";
}

/*
void
WarpX::applyFilter (MultiFab& dstmf, const MultiFab& srcmf, int scomp, int dcomp, int ncomp)
{
    ncomp = std::min(ncomp, srcmf.nComp());
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox tmpfab;
        for (MFIter mfi(dstmf,true); mfi.isValid(); ++mfi)
            {
                const auto& srcfab = srcmf[mfi];
                auto& dstfab = dstmf[mfi];
                const Box& tbx = mfi.growntilebox();
                const Box& gbx = amrex::grow(tbx,1);
                tmpfab.resize(gbx,ncomp);
                tmpfab.setVal(0.0, gbx, 0, ncomp);
                const Box& ibx = gbx & srcfab.box();
                tmpfab.copy(srcfab, ibx, scomp, ibx, 0, ncomp);
                WRPX_FILTER(BL_TO_FORTRAN_BOX(tbx),
                            BL_TO_FORTRAN_ANYD(tmpfab),
                            BL_TO_FORTRAN_N_ANYD(dstfab,dcomp),
                            ncomp);
            }
    }
}

*/
 /*
void BilinearFilter::ApplyStencils(amrex::MultiFab* dstmf, amrex::MultiFab* srcmf){
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        //FArrayBox tmpfab;
        for ( MFIter mfi(*dstmf, TilingIfNotGPU()); mfi.isValid(); ++mfi ){
            //const auto& srcfab = srcmf[mfi];
            //auto& dstfab = dstmf[mfi];
            const Box& tbx = mfi.growntilebox();
            //const Box& gbx = amrex::grow(tbx,1);
            //tmpfab.resize(gbx,ncomp);
            //tmpfab.setVal(0.0, gbx, 0, ncomp);
            //const Box& ibx = gbx & srcfab.box();
            //tmpfab.copy(srcfab, ibx, scomp, ibx, 0, ncomp);
            auto const& dstmf_fab = dstmf->array(mfi);
            auto const& srcmf_fab = srcmf->array(mfi);

            amrex::ParallelFor(tbx,
                [=] AMREX_GPU_DEVICE (int j, int k, int l){
                warpx_apply_filter_2d(j,k,l,dstmf_fab, srcmf_fab,
                                      stencil_x,stencil_z,
                                      stencil_length_each_dir[0],
                                      stencil_length_each_dir[1]);
                               });
        }
    }
}
 */
