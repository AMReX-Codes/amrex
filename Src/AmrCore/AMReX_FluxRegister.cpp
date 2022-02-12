
#include <AMReX_BArena.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

FluxRegister::FluxRegister ()
{
    fine_level = ncomp = -1;
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
}

FluxRegister::FluxRegister (const BoxArray&            fine_boxes,
                            const DistributionMapping& dm,
                            const IntVect&             ref_ratio,
                            int                        fine_lev,
                            int                        nvar)
{
    define(fine_boxes,dm,ref_ratio,fine_lev,nvar);
}

const IntVect&
FluxRegister::refRatio () const noexcept
{
    return ratio;
}

int
FluxRegister::fineLevel () const noexcept
{
    return fine_level;
}

int
FluxRegister::crseLevel () const noexcept
{
    return fine_level-1;
}

int
FluxRegister::nComp () const noexcept
{
    return ncomp;
}

const BoxArray&
FluxRegister::coarsenedBoxes () const noexcept
{
    return grids;
}

void
FluxRegister::define (const BoxArray&            fine_boxes,
                      const DistributionMapping& dm,
                      const IntVect&             ref_ratio,
                      int                        fine_lev,
                      int                        nvar)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids = fine_boxes;
    grids.coarsen(ratio);

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

        IndexType typ(IndexType::TheCellType());

        typ.setType(dir,IndexType::NODE);

        BndryRegister::define(lo_face,typ,0,1,0,nvar,dm);
        BndryRegister::define(hi_face,typ,0,1,0,nvar,dm);
    }
}

void
FluxRegister::clear ()
{
    BndryRegister::clear();
}

FluxRegister::~FluxRegister () {}

Real
FluxRegister::SumReg (int comp) const
{
    Real sum = 0.0;

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        const FabSet& lofabs = bndry[Orientation(dir,Orientation::low) ];
        const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            sum += amrex::ReduceSum(lofabs.m_mf, 0,
                   [=] AMREX_GPU_HOST_DEVICE (Box const& tbx,
                                              Array4<Real const> const& lofab) -> Real
                   {
                       Real r = 0.0;
                       AMREX_LOOP_3D(tbx, i, j, k,
                       {
                           r += lofab(i,j,k,comp);
                       });
                       return r;
                   });
            sum += amrex::ReduceSum(hifabs.m_mf, 0,
                   [=] AMREX_GPU_HOST_DEVICE (Box const& tbx,
                                              Array4<Real const> const& hifab) -> Real
                   {
                       Real r = 0.0;
                       AMREX_LOOP_3D(tbx, i, j, k,
                       {
                           r -= hifab(i,j,k,comp);
                       });
                       return r;
                   });
        }
        else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(+:sum)
#endif
            for (FabSetIter fsi(lofabs); fsi.isValid(); ++fsi)
            {
                Array4<Real const> const& lofab = lofabs.const_array(fsi);
                Box lobx(lofab);
                AMREX_LOOP_3D(lobx, i, j, k,
                {
                    sum += lofab(i,j,k,comp);
                });
                Array4<Real const> const& hifab = hifabs.const_array(fsi);
                Box hibx(hifab);
                AMREX_LOOP_3D(hibx, i, j, k,
                {
                    sum -= hifab(i,j,k,comp);
                });
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        const MultiFab& area,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);

    MultiFab mf(mflx.boxArray(),mflx.DistributionMap(),numcomp,0,
                MFInfo(), mflx.Factory());

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mflx.isFusingCandidate()) {
        auto const& dma = mf.arrays();
        auto const& sma = mflx.const_arrays();
        auto const& ama = area.const_arrays();
        ParallelFor(mf, IntVect(0), numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dma[box_no](i,j,k,n) = sma[box_no](i,j,k,n+srccomp)*mult*ama[box_no](i,j,k);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mflx,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto       dfab =   mf.array(mfi);
            auto const sfab = mflx.const_array(mfi);
            auto const afab = area.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dfab(i,j,k,n) = sfab(i,j,k,n+srccomp)*mult*afab(i,j,k);
            });
        }
    }

    for (int pass = 0; pass < 2; pass++)
    {
        const Orientation face = ((pass == 0) ? face_lo : face_hi);

        if (op == FluxRegister::COPY)
        {
            bndry[face].copyFrom(mf,0,0,destcomp,numcomp);
        }
        else
        {
            FabSet fs(bndry[face].boxArray(),bndry[face].DistributionMap(),numcomp);

            fs.setVal(0);

            fs.copyFrom(mf,0,0,0,numcomp);

#ifdef AMREX_USE_GPU
            using Tag = Array4PairTag<Real>;
            Vector<Tag> tags;
            tags.reserve(mflx.local_size()*AMREX_SPACEDIM*2);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                auto const sfab = fs.const_array(mfi);
                auto       dfab = bndry[face].array(mfi);
#ifdef AMREX_USE_GPU
                if (Gpu::inLaunchRegion()) {
                    tags.push_back({dfab, sfab, bx});
                } else
#endif
                {
                    AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+destcomp) += sfab(i,j,k,n);
                    });
                }
            }

#ifdef AMREX_USE_GPU
            ParallelFor(tags, numcomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Tag const& tag) noexcept
            {
                tag.dfab(i,j,k,n+destcomp) += tag.sfab(i,j,k,n);
            });
#endif
        }
    }
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    MultiFab area(mflx.boxArray(), mflx.DistributionMap(), 1, 0,
                  MFInfo(), mflx.Factory());

    area.setVal(1, 0, 1, 0);

    CrseInit(mflx,area,dir,srccomp,destcomp,numcomp,mult,op);
}

//void
//FluxRegister::CrseInit_DG ( const MultiFab& SurfaceFlux,
//                            int             iDimX,
//                            int             SrcComp,
//                            int             DestComp,
//                            int             nComp,
//                            Real            mult,
//                            MultiFab&       FluxIncrement,
//                            int             nDOFX_K,
//                            int             nNodesX1,
//                            Real            LX_Dn_1D[],
//                            Real            LX_Up_1D[],
//                            Real            WeightsX[],
//                            Real            dX[],
//                            FrOp            op)
//{
//
//    /* Consistency checks */
////    BL_ASSERT( SrcComp  >= 0 && SrcComp  + nComp <= SurfaceFlux.nComp() );
////    BL_ASSERT( DestComp >= 0 && DestComp + nComp <= ncomp );
//
//    /* face_lo = (0), face_hi = (1) */
//    const Orientation face_lo( iDimX, Orientation::low  );
//    const Orientation face_hi( iDimX, Orientation::high );
//
//    /* Define destination MultiFab */
//    MultiFab mf( FluxIncrement.boxArray(), FluxIncrement.DistributionMap(),
//                 nComp, 0, MFInfo(), FluxIncrement.Factory() );
////    mf.setVal( (Real)0.0 );
//
//const Real One = 1.0;
//
///* Send these in */
//const int nCF = 6;
//const int swX[3] = { 1, 0, 0 };
//int nNodesX[3] = { 2, 1, 1 };
//const int nDOFX_X1 = 1;
//const int nDOFX_X2 = 1;
//const int nDOFX_X3 = 1;
//const int nDOFX_X = nDOFX_X1;
//const int nDOFX = 2;
///* */
//
//int iNX;
//int iNX1;
//int iNX2;
//int iNX3;
//int iNX_X;
//int iNX_X1;
//int iNX_X2;
//int iNX_X3;
//
//int iX_B0[3];
//int iX_E0[3];
//int iX1, iX2, iX3;
//
//int iCF;
//
//int NodeNumberTableX[3][nDOFX];
//int NodeNumberTableX3D[2][2][2];
//int NodeNumberTableX_X[nDOFX];
//
//iNX = -1;
//for( iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
//for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
//for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
//  iNX += 1;
//  NodeNumberTableX[0][iNX] = iNX1;
//  NodeNumberTableX[1][iNX] = iNX2;
//  NodeNumberTableX[2][iNX] = iNX3;
//  NodeNumberTableX3D[iNX1][iNX2][iNX3] = iNX;
//}}}
//
///* Break these up into dimension-dependent individual functions */
//if( iDimX == 0 )
//{
//  for( iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
//  for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
//  for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
//    iNX = NodeNumberTableX3D[iNX1][iNX2][iNX3];
//    NodeNumberTableX_X[iNX] = iNX / nNodesX[0];
//  }}}
//}
//if( iDimX == 1 )
//{
//  for( iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
//  for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
//  for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
//    iNX = NodeNumberTableX3D[iNX1][iNX2][iNX3];
//    NodeNumberTableX_X[iNX] = iNX % nNodesX[0] + nNodesX[0] * iNX3;
//  }}}
//}
//if( iDimX == 2 )
//{
//  for( iNX3 = 0; iNX3 < nNodesX[2]; iNX3++ ) {
//  for( iNX2 = 0; iNX2 < nNodesX[1]; iNX2++ ) {
//  for( iNX1 = 0; iNX1 < nNodesX[0]; iNX1++ ) {
//    iNX = NodeNumberTableX3D[iNX1][iNX2][iNX3];
//    NodeNumberTableX_X[iNX] = iNX % nDOFX_X1;
//  }}}
//}
//
//    /* Populate destination MultiFab */
//#ifdef AMREX_USE_OMP
//#pragma omp parallel if (Gpu::notInLaunchRegion())
//#endif
//    for ( MFIter mfi( SurfaceFlux, TilingIfNotGPU() ); mfi.isValid(); ++mfi )
//    {
//        const Box& bx   = mfi.tilebox();
//        auto       dfab = mf.array( mfi );
//        auto const sfab = SurfaceFlux.const_array( mfi );
//
//        AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, nComp, i, j, k, n,
//        {
////          dfab(i,j,k,n) = sfab(i,j,k,n+SrcComp) * mult;
//          for( iCF = 0; iCF < nCF; iCF++ )
//          {
//            for( iNX = 0; iNX < nDOFX; iNX++ )
//              {
//
//                iNX1 = NodeNumberTableX[0][iNX];
//
//                iNX_X1 = NodeNumberTableX_X[iNX];
//
//                dfab(i,j,k,iNX+iCF*nDOFX)
//                  = sfab(iX1,iX2,iX3,iNX_X1+iCF*nDOFX_X1) * mult;
//
//              } /* iNX */
//          } /* iCF */
//
//        });
//
////        /* bx.loVect() has only amrex_spacedim components */
////
////        iX_B0[0] = bx.loVect()[0];
////        iX_E0[0] = bx.hiVect()[0];
////
////        iX_B0[1] = 0;
////        iX_E0[1] = 1;
////        #if( BL_SPACEDIM > 1 )
////        iX_B0[1] = bx.loVect()[1];
////        iX_E0[1] = bx.hiVect()[1];
////        #endif
////
////        iX_B0[2] = 0;
////        iX_E0[2] = 1;
////        #if( BL_SPACEDIM > 2 )
////        iX_B0[2] = bx.loVect()[2];
////        iX_E0[2] = bx.hiVect()[2];
////        #endif
////
////        for( iX1 = iX_B0[0]; iX1 < iX_E0[0]; iX1++ ) {
////        for( iX2 = iX_B0[1]; iX2 < iX_E0[1]; iX2++ ) {
////        for( iX3 = iX_B0[2]; iX3 < iX_E0[2]; iX3++ ) {
////
////          for( iCF = 0; iCF < nCF; iCF++ )
////          {
////            for( iNX = 0; iNX < nDOFX; iNX++ )
////              {
////
////                iNX1 = NodeNumberTableX[0][iNX];
////
////                iNX_X1 = NodeNumberTableX_X[iNX];
////
//////                dfab(iX1,iX2,iX3,iNX+iCF*nDOFX)
//////                  = -One / ( WeightsX[iNX1] * dX[iX1] )
//////                      * (   LX_Up_1D[iNX1]
//////                              * sfab(iX1+1,iX2,iX3,iNX_X1+iCF*nDOFX_X1)
//////                          - LX_Dn_1D[iNX1]
//////                              * sfab(iX1  ,iX2,iX3,iNX_X1+iCF*nDOFX_X1) );
////
////                dfab(iX1,iX2,iX3,iNX+iCF*nDOFX)
////                  = +One / ( WeightsX[iNX1] * dX[iX1] )
////                      * LX_Dn_1D[iNX1] * sfab(iX1,iX2,iX3,iNX_X1+iCF*nDOFX_X1);
////
////              } /* iNX */
////          } /* iCF */
////
////        }}}
//
//    } /* END for MFIter  */
//
//    /* pass <==> which side of face */
//    for ( int pass = 0; pass < 2; pass++ )
//    {
//        /* if pass == 0 then face = face_lo; else face = face_hi */
//        const Orientation face = ( ( pass == 0 ) ? face_lo : face_hi );
//
//        if ( op == FluxRegister::COPY )
//        {
//            bndry[face].copyFrom( mf, 0, 0, DestComp, nComp );
//        }
//// This `else` never happens because `op` is always set to the default
//// (i.e., FluxRegister::COPY) in the Fortran interface, so it has been
//// commented out
//        else
//        {
//std::cout<<"THIS SHOULD NEVER PRINT!\n";
///*
//            FabSet fs( bndry[face].boxArray(), bndry[face].DistributionMap(),
//                       nComp );
//
//            fs.setVal( 0 );
//
//            fs.copyFrom( mf, 0, 0, 0, nComp );
//
//#ifdef AMREX_USE_OMP
//#pragma omp parallel if (Gpu::notInLaunchRegion())
//#endif
//            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
//            {
//                const Box& bx = mfi.validbox();
//                auto const sfab = fs.const_array(mfi);
//                auto       dfab = bndry[face].array(mfi);
//                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (bx, nComp, i, j, k, n,
//                {
//                    dfab(i,j,k,n+DestComp) += sfab(i,j,k,n);
//                });
//            }
//*/
//        }
//    } /* END for ( int pass = 0; pass < 2; pass++ ) */
//
//} /* END void FluxRegister::CrseInit_DG */
//
void
FluxRegister::CrseInit_DG ( const MultiFab& SurfaceFlux,
                            int             iDimX,
                            int             nFields,
                            int             nDOFX_X,
                            Real            WeightsX_X[],
                            int             SrcComp,
                            int             DestComp,
                            Real            mult,
                            FrOp            op)
{

    int swX[3];
    if( iDimX == 0 ) { swX[0] = 1; swX[1] = 0; swX[2] = 0; }
    if( iDimX == 1 ) { swX[0] = 0; swX[1] = 1; swX[2] = 0; }
    if( iDimX == 2 ) { swX[0] = 0; swX[1] = 0; swX[2] = 1; }

    int nComp = nDOFX_X * nFields;

    /* Consistency checks */
    BL_ASSERT( SrcComp  >= 0 && SrcComp  + nComp <= SurfaceFlux.nComp() );
    BL_ASSERT( DestComp >= 0 && DestComp + nComp <= nComp );

    /* face_lo = (0), face_hi = (1) */
    const Orientation face_lo( iDimX, Orientation::low  );
    const Orientation face_hi( iDimX, Orientation::high );

    /* Define destination MultiFab */
    MultiFab mf( SurfaceFlux.boxArray(), SurfaceFlux.DistributionMap(),
                 nComp, 0, MFInfo(), SurfaceFlux.Factory() );
    mf.setVal( 0.0 );

    int iX_B0[3];
    int iX_E0[3];

    /* Populate destination MultiFab */
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi( SurfaceFlux, TilingIfNotGPU() ); mfi.isValid(); ++mfi )
    {
        const Box& bx   = mfi.tilebox();
        auto       dfab = mf.array( mfi );
        auto const sfab = SurfaceFlux.const_array( mfi );

//        AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, nComp, i, j, k, n,
//        {
//          dfab(i,j,k,n) = sfab(i,j,k,n+SrcComp) * mult;
//        });

        /* bx.loVect() has only amrex_spacedim components */

        iX_B0[0] = bx.loVect()[0];
        iX_E0[0] = bx.hiVect()[0];

        iX_B0[1] = 0;
        iX_E0[1] = 1;
        if( AMREX_SPACEDIM > 1 ) {
        iX_B0[1] = bx.loVect()[1];
        iX_E0[1] = bx.hiVect()[1];
        }

        iX_B0[2] = 0;
        iX_E0[2] = 1;
        if( AMREX_SPACEDIM > 2 ) {
        iX_B0[2] = bx.loVect()[2];
        iX_E0[2] = bx.hiVect()[2];
        }

        for( int i = iX_B0[0]; i < iX_E0[0]+swX[0]; i++ ) {
        for( int j = iX_B0[1]; j < iX_E0[1]+swX[1]; j++ ) {
        for( int k = iX_B0[2]; k < iX_E0[2]+swX[2]; k++ ) {

            for( int iField = 0; iField < nFields; iField++ )
            {
                for( int iNX_X = 0; iNX_X < nDOFX_X; iNX_X++ )
                {
                  dfab(i,j,k,iNX_X+iField*nDOFX_X)
                    = mult * WeightsX_X[iNX_X]
                        * sfab(i,j,k,SrcComp+iNX_X+iField*nDOFX_X);
                } /* iNX_X */
            } /* iField */

        }}} /* i, j, k */

    } /* END for MFIter  */

    /* pass <==> which side of face */
    for ( int pass = 0; pass < 2; pass++ )
    {
        /* if pass == 0 then face = face_lo; else face = face_hi */
        const Orientation face = ( ( pass == 0 ) ? face_lo : face_hi );

        if ( op == FluxRegister::COPY )
        {
            bndry[face].copyFrom( mf, 0, 0, DestComp, nComp );
        }
// This `else` never happens because `op` is always set to the default
// (i.e., FluxRegister::COPY) in the Fortran interface, so it has been
// commented out
        else
        {
std::cout<<"THIS SHOULD NEVER PRINT!\n";
/*
            FabSet fs( bndry[face].boxArray(), bndry[face].DistributionMap(),
                       nComp );

            fs.setVal( 0 );

            fs.copyFrom( mf, 0, 0, 0, nComp );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                auto const sfab = fs.const_array(mfi);
                auto       dfab = bndry[face].array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (bx, nComp, i, j, k, n,
                {
                    dfab(i,j,k,n+DestComp) += sfab(i,j,k,n);
                });
            }
*/
        }
    } /* END for ( int pass = 0; pass < 2; pass++ ) */

} /* END void FluxRegister::CrseInit_DG */

void
FluxRegister::CrseAdd (const MultiFab& mflx,
                       const MultiFab& area,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult,
                       const Geometry& geom)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);

    MultiFab mf(mflx.boxArray(),mflx.DistributionMap(),numcomp,0,
                MFInfo(), mflx.Factory());

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mflx.isFusingCandidate()) {
        auto const& dma = mf.arrays();
        auto const& sma = mflx.const_arrays();
        auto const& ama = area.const_arrays();
        ParallelFor(mf, IntVect(0), numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dma[box_no](i,j,k,n) = sma[box_no](i,j,k,n+srccomp)*mult*ama[box_no](i,j,k);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mflx,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto       dfab =   mf.array(mfi);
            auto const sfab = mflx.const_array(mfi);
            auto const afab = area.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dfab(i,j,k,n) = sfab(i,j,k,n+srccomp)*mult*afab(i,j,k);
            });
        }
    }

    for (int pass = 0; pass < 2; pass++)
    {
        const Orientation face = ((pass == 0) ? face_lo : face_hi);
        bndry[face].plusFrom(mf,0,0,destcomp,numcomp,geom.periodicity());
    }
}

void
FluxRegister::CrseAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult,
                       const Geometry& geom)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    MultiFab area(mflx.boxArray(), mflx.DistributionMap(), 1, 0,
                  MFInfo(), mflx.Factory());

    area.setVal(1, 0, 1, 0);

    CrseAdd(mflx,area,dir,srccomp,destcomp,numcomp,mult,geom);
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],dir,k,srccomp,destcomp,numcomp,mult,RunOn::Gpu);
    }
}

void
FluxRegister::FineAdd_DG (const MultiFab& SurfaceFluxes,
                          int             iDimX,
                          int             nF,
                          int             nDOFX_X,
                          Real            WeightsX_X[],
                          Real            LX_X[],
                          int             SrcComp,
                          int             DestComp,
                          Real            mult)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(SurfaceFluxes); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd_DG( SurfaceFluxes[mfi], iDimX, k, nF, nDOFX_X,
                    WeightsX_X, LX_X, SrcComp, DestComp, mult, RunOn::Gpu);
    }
} /* END void FluxRegister::FineAdd_DG */

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       const MultiFab& area,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],area[mfi],dir,k,srccomp,destcomp,numcomp,mult,RunOn::Gpu);
    }
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult,
                       RunOn            runon) noexcept
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];
    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];
    const Box& lobox = loreg.box();
    const Box& hibox = hireg.box();

    Array4<Real> loarr = loreg.array();
    Array4<Real> hiarr = hireg.array();
    Array4<Real const> farr = flux.const_array();
    const Dim3 local_ratio = ratio.dim3();

    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion())
    {
        AMREX_LAUNCH_DEVICE_LAMBDA
        ( lobox, tlobx,
          {
              fluxreg_fineadd(tlobx, loarr, destcomp, farr, srccomp,
                              numcomp, dir, local_ratio, mult);
          },
          hibox, thibx,
          {
              fluxreg_fineadd(thibx, hiarr, destcomp, farr, srccomp,
                              numcomp, dir, local_ratio, mult);
          }
        );
    }
    else
    {
        fluxreg_fineadd(lobox, loarr, destcomp, farr, srccomp,
                        numcomp, dir, local_ratio, mult);
        fluxreg_fineadd(hibox, hiarr, destcomp, farr, srccomp,
                        numcomp, dir, local_ratio, mult);
    }
}

void
FluxRegister::FineAdd_DG (const FArrayBox& SurfaceFluxes,
                          int              iDimX,
                          int              boxno,
                          int              nF,
                          int              nDOFX_X,
                          Real             WeightsX_X[],
                          Real             LX_X[],
                          int              SrcComp,
                          int              DestComp,
                          Real             mult,
                          RunOn            runon) noexcept
{
//    /* Convert LX_X into 2D array */
//
//    Real ** LX_XX = new Real * [nDOFX_X];
//    for( int iNX = 0; iNX < nDOFX_X; iNX++ ) { LX_XX[iNX] = new Real[nDOFX_X]; }
//
//    int k = -1;
//    for( int iNX_X = 0; iNX_X < nDOFX_X                      ; iNX_X++ ) {
//    for( int iNX   = 0; iNX   < pow(2,AMREX_SPACEDIM)*nDOFX_X; iNX++   ) {
//        k += 1;
//        LX_XX[iNX_X][iNX] = LX_X[k];
//    }}

    int nComp = nDOFX_X * nF;

    BL_ASSERT( SrcComp  >= 0 && SrcComp  + nComp <= SurfaceFluxes.nComp() );
    BL_ASSERT( DestComp >= 0 && DestComp + nComp <= ncomp );

    FArrayBox& loreg = bndry[Orientation(iDimX,Orientation::low)][boxno];
    FArrayBox& hireg = bndry[Orientation(iDimX,Orientation::high)][boxno];
    const Box& lobox = loreg.box();
    const Box& hibox = hireg.box();

    Array4<Real> loarr = loreg.array();
    Array4<Real> hiarr = hireg.array();
    Array4<Real const> farr = SurfaceFluxes.const_array();
    const Dim3 local_ratio = ratio.dim3();

    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion())
    {
        AMREX_LAUNCH_DEVICE_LAMBDA
        ( lobox, tlobx,
          {
              fluxreg_fineadd_dg(tlobx, loarr, farr,
                                 iDimX, nF, nDOFX_X, WeightsX_X, LX_X,
                                 local_ratio, mult);
          },
          hibox, thibx,
          {
              fluxreg_fineadd_dg(thibx, hiarr, farr,
                                 iDimX, nF, nDOFX_X, WeightsX_X, LX_X,
                                 local_ratio, mult);
          }
        );
    }
    else
    {
        fluxreg_fineadd_dg(lobox, loarr, farr,
                           iDimX, nF, nDOFX_X, WeightsX_X, LX_X,
                           local_ratio, mult);
        fluxreg_fineadd_dg(hibox, hiarr, farr,
                           iDimX, nF, nDOFX_X, WeightsX_X, LX_X,
                           local_ratio, mult);
    }

//    for( int iNX = 0; iNX < nDOFX_X; iNX++ ) { delete LX_XX[iNX]; }
//    delete LX_XX;

} /* END void FluxRegister::FineAdd_DG */

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       const FArrayBox& area,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult,
                       RunOn            runon) noexcept
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];
    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];
    const Box& lobox = loreg.box();
    const Box& hibox = hireg.box();

    Array4<Real> loarr = loreg.array();
    Array4<Real> hiarr = hireg.array();
    Array4<Real const> farr = flux.const_array();
    Array4<Real const> aarr = area.const_array();
    const Dim3 local_ratio = ratio.dim3();

    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion())
    {
        AMREX_LAUNCH_DEVICE_LAMBDA
        ( lobox, tlobx,
          {
              fluxreg_fineareaadd(tlobx, loarr, destcomp,
                                  aarr, farr, srccomp,
                                  numcomp, dir, local_ratio, mult);
          },
          hibox, thibx,
          {
              fluxreg_fineareaadd(thibx, hiarr, destcomp,
                                  aarr, farr, srccomp,
                                  numcomp, dir, local_ratio, mult);
          }
        );
    }
    else
    {
        fluxreg_fineareaadd(lobox, loarr, destcomp,
                            aarr, farr, srccomp,
                            numcomp, dir, local_ratio, mult);
        fluxreg_fineareaadd(hibox, hiarr, destcomp,
                            aarr, farr, srccomp,
                            numcomp, dir, local_ratio, mult);
    }
}

void
FluxRegister::FineSetVal (int              dir,
                          int              boxno,
                          int              destcomp,
                          int              numcomp,
                          Real             val,
                          RunOn            runon) noexcept
{
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];
    BL_ASSERT(numcomp <= loreg.nComp());
    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion()) {
        loreg.setVal<RunOn::Device>(val, loreg.box(), destcomp, numcomp);
    } else {
        loreg.setVal<RunOn::Host>(val, loreg.box(), destcomp, numcomp);
    }

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];
    BL_ASSERT(numcomp <= hireg.nComp());
    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion()) {
        hireg.setVal<RunOn::Device>(val, hireg.box(), destcomp, numcomp);
    } else {
        hireg.setVal<RunOn::Host>(val, hireg.box(), destcomp, numcomp);
    }
}

void
FluxRegister::Reflux (MultiFab&       mf,
                      const MultiFab& volume,
                      Real            scale,
                      int             scomp,
                      int             dcomp,
                      int             nc,
                      const Geometry& geom)
{
    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation& face = fi();
        Reflux(mf, volume, face, scale, scomp, dcomp, nc, geom);
    }
}

void
FluxRegister::Reflux (MultiFab&       mf,
                      const MultiFab& volume,
                      int             dir,
                      Real            scale,
                      int             scomp,
                      int             dcomp,
                      int             nc,
                      const Geometry& geom)
{
    for (int s = 0; s < 2; ++s)
    {
        Orientation::Side side = (s==0) ? Orientation::low : Orientation::high;
        Orientation face(dir, side);
        Reflux(mf, volume, face, scale, scomp, dcomp, nc, geom);
    }
}

void
FluxRegister::Reflux_DG ( MultiFab&       mf,
                          int             nComp,
                          const Geometry& geom,
                          int             nFields,
                          int             nDOFX,
                          int             nNodesX[],
                          Real            WeightsX_q[],
                          Real            dX1[],
                          Real            dX2[],
                          Real            dX3[],
                          Real            LX_X1_Dn[],
                          Real            LX_X1_Up[],
                          Real            LX_X2_Dn[],
                          Real            LX_X2_Up[],
                          Real            LX_X3_Dn[],
                          Real            LX_X3_Up[] )
{
    for( OrientationIter fi; fi; ++fi )
    {
        const Orientation& face = fi();
        Reflux_DG( mf, face, nComp, geom,
                   nFields, nDOFX, nNodesX, WeightsX_q, dX1, dX2, dX3,
                   LX_X1_Dn, LX_X1_Up, LX_X2_Dn, LX_X2_Up, LX_X3_Dn, LX_X3_Up );
    }
} /* END void FluxRegister::Reflux_DG */

void
FluxRegister::Reflux (MultiFab&       mf,
                      Real            scale,
                      int             scomp,
                      int             dcomp,
                      int             nc,
                      const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    MultiFab volume(mf.boxArray(), mf.DistributionMap(), 1, 0,
                    MFInfo(), mf.Factory());

    volume.setVal(AMREX_D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, 0);

    Reflux(mf,volume,scale,scomp,dcomp,nc,geom);
}

void
FluxRegister::Reflux (MultiFab&       mf,
                      int             dir,
                      Real            scale,
                      int             scomp,
                      int             dcomp,
                      int             nc,
                      const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    MultiFab volume(mf.boxArray(), mf.DistributionMap(), 1, 0,
                    MFInfo(), mf.Factory());

    volume.setVal(AMREX_D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, 0);

    Reflux(mf,volume,dir,scale,scomp,dcomp,nc,geom);
}

void
FluxRegister::Reflux (MultiFab& mf, const MultiFab& volume, Orientation face,
                      Real scale, int scomp, int dcomp, int nc, const Geometry& geom)
{
    BL_PROFILE("FluxRegister::Reflux()");

    int idir = face.coordDir();

    MultiFab flux(amrex::convert(mf.boxArray(), IntVect::TheDimensionVector(idir)),
                  mf.DistributionMap(), nc, 0, MFInfo(), mf.Factory());
    flux.setVal(0.0);

    bndry[face].copyTo(flux, 0, scomp, 0, nc, geom.periodicity());

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mf.isFusingCandidate()) {
        auto const& sma = mf.arrays();
        auto const& fma = flux.const_arrays();
        auto const& vma = volume.const_arrays();
        ParallelFor(mf, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            fluxreg_reflux(Box(IntVect(AMREX_D_DECL(i,j,k)),IntVect(AMREX_D_DECL(i,j,k))),
                           sma[box_no], dcomp, fma[box_no], vma[box_no], nc, scale, face);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& sfab = mf.array(mfi);
            Array4<Real const> const& ffab = flux.const_array(mfi);
            Array4<Real const> const& vfab = volume.const_array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
            {
                fluxreg_reflux(tbx, sfab, dcomp, ffab, vfab, nc, scale, face);
            });
        }
    }
}

void
FluxRegister::Reflux_DG
  ( MultiFab& mf, Orientation face,
    int nComp, const Geometry& geom,
    int nFields, int nDOFX, int nNodesX[],
    Real WeightsX_q[], Real dX1[], Real dX2[], Real dX3[],
    Real LX_X1_Dn[], Real LX_X1_Up[],
    Real LX_X2_Dn[], Real LX_X2_Up[],
    Real LX_X3_Dn[], Real LX_X3_Up[] )
{
    BL_PROFILE("FluxRegister::Reflux_DG()");

    int iDimX = face.coordDir();

    /* What is this doing? */
    MultiFab flux( amrex::convert( mf.boxArray(),
                                   IntVect::TheDimensionVector(iDimX) ),
                   mf.DistributionMap(), nComp, 0, MFInfo(), mf.Factory() );

    bndry[face].copyTo( flux, 0, 0, 0, nComp, geom.periodicity() );

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for( MFIter mfi( mf, TilingIfNotGPU() ); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& u = mf.array(mfi);
        Array4<Real const> const& ffab = flux.const_array(mfi);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
        {
            fluxreg_reflux_dg( tbx, u, ffab, face,
                               nFields, nDOFX, nNodesX,
                               WeightsX_q, dX1, dX2, dX3,
                               LX_X1_Dn, LX_X1_Up,
                               LX_X2_Dn, LX_X2_Up,
                               LX_X3_Dn, LX_X3_Up );
        });
    }
} /* END void FluxRegister::Reflux_DG */

void
FluxRegister::ClearInternalBorders (const Geometry& geom)
{
    int nc = this->nComp();
    const Box& domain = geom.Domain();

    Vector<Array4BoxTag<Real> > tags;

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        Orientation lo(dir, Orientation::low);
        Orientation hi(dir, Orientation::high);

        FabSet& frlo = (*this)[lo];
        FabSet& frhi = (*this)[hi];

        const BoxArray& balo = frlo.boxArray();
        const BoxArray& bahi = frhi.boxArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (FabSetIter fsi(frlo); fsi.isValid(); ++fsi) {
                const Box& bx = fsi.validbox();
                const std::vector< std::pair<int,Box> >& isects = bahi.intersections(bx);
                auto const& frarr = frlo[fsi].array();
                for (int ii = 0; ii < static_cast<int>(isects.size()); ++ii) {
                    if (Gpu::inLaunchRegion()) {
                        tags.emplace_back(Array4BoxTag<Real>{frarr, isects[ii].second});
                    } else {
                        frlo[fsi].setVal<RunOn::Host>(0.0, isects[ii].second, 0, nc);
                    }
                }
                if (geom.isPeriodic(dir)) {
                    if (bx.smallEnd(dir) == domain.smallEnd(dir)) {
                        const Box& sbx = amrex::shift(bx, dir, domain.length(dir));
                        const std::vector<std::pair<int,Box> >& isects2 = bahi.intersections(sbx);
                        for (int ii = 0; ii < static_cast<int>(isects2.size()); ++ii) {
                            const Box& bx2 = amrex::shift(isects2[ii].second, dir, -domain.length(dir));
                            if (Gpu::inLaunchRegion()) {
                                tags.emplace_back(Array4BoxTag<Real>{frarr, bx2});
                            } else {
                                frlo[fsi].setVal<RunOn::Host>(0.0, bx2, 0, nc);
                            }
                        }
                    }
                }
            }

            for (FabSetIter fsi(frhi); fsi.isValid(); ++fsi) {
                const Box& bx = fsi.validbox();
                const std::vector< std::pair<int,Box> >& isects = balo.intersections(bx);
                auto const& frarr = frhi[fsi].array();
                for (int ii = 0; ii < static_cast<int>(isects.size()); ++ii) {
                    if (Gpu::inLaunchRegion()) {
                        tags.emplace_back(Array4BoxTag<Real>{frarr, isects[ii].second});
                    } else {
                        frhi[fsi].setVal<RunOn::Host>(0.0, isects[ii].second, 0, nc);
                    }
                }
                if (geom.isPeriodic(dir)) {
                    if (bx.bigEnd(dir) == domain.bigEnd(dir)) {
                        const Box& sbx = amrex::shift(bx, dir, -domain.length(dir));
                        const std::vector<std::pair<int,Box> >& isects2 = balo.intersections(sbx);
                        for (int ii = 0; ii < static_cast<int>(isects2.size()); ++ii) {
                            const Box& bx2 = amrex::shift(isects2[ii].second, dir, domain.length(dir));
                            if (Gpu::inLaunchRegion()) {
                                tags.emplace_back(Array4BoxTag<Real>{frarr, bx2});
                            } else {
                                frhi[fsi].setVal<RunOn::Host>(0.0, bx2, 0, nc);
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    // There is Gpu::synchronize in Parallelfor below internally.
    ParallelFor(tags, nc,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Array4BoxTag<Real> const& tag)
    {
        tag.dfab(i,j,k,n) = 0.0;
    });
#endif
}

void
FluxRegister::OverwriteFlux (Array<MultiFab*,AMREX_SPACEDIM> const& crse_fluxes,
                             Real scale, int srccomp, int destcomp, int numcomp,
                             const Geometry& crse_geom)
{
    BL_PROFILE("FluxRegister::OverwriteFlux()");

    const auto& cperiod = crse_geom.periodicity();

    Box cdomain = crse_geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (crse_geom.isPeriodic(idim)) cdomain.grow(idim, 1);
    }

    bool run_on_gpu = Gpu::inLaunchRegion();

    // cell-centered mask:
    constexpr int crse_cell = 0;
    constexpr int fine_cell = 1; // covered by fine
    constexpr int phbc_cell = 2;
    const BoxArray& cba = amrex::convert(crse_fluxes[0]->boxArray(), IntVect::TheCellVector());
    iMultiFab cc_mask(cba, crse_fluxes[0]->DistributionMap(), 1, 1);

    bool inited = false;
#ifdef AMREX_USE_GPU
    if (run_on_gpu && cc_mask.isFusingCandidate()) {
        auto const& ma = cc_mask.arrays();
        ParallelFor(cc_mask, IntVect(1),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            if (cdomain.contains(i,j,k)) {
                ma[box_no](i,j,k) = crse_cell;
            } else {
                ma[box_no](i,j,k) = phbc_cell;
            }
        });
        // No need to call streamSynchronize here
        inited = true;
    }
#endif

    {
        const std::vector<IntVect>& pshifts = cperiod.shiftIntVect();
        Vector<Array4BoxTag<int> > tags;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (!run_on_gpu)
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
            {
                auto const& fab = cc_mask.array(mfi);
                const Box& bx = mfi.fabbox();

                if (!inited) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        if (cdomain.contains(i,j,k)) {
                            fab(i,j,k) = crse_cell;
                        } else {
                            fab(i,j,k) = phbc_cell;
                        }
                    });
                }

                for (const auto& iv : pshifts)
                {
                    grids.intersections(bx+iv, isects);
                    for (const auto& is : isects)
                    {
                        Box const& b = is.second-iv;
                        if (run_on_gpu) {
                            tags.push_back({fab,b});
                        } else {
                            cc_mask[mfi].setVal<RunOn::Host>(fine_cell, b, 0, 1);
                        }
                    }
                }
            }
        }

#ifdef AMREX_USE_GPU
        amrex::ParallelFor(tags, 1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                              Array4BoxTag<int> const& tag) noexcept
        {
            tag.dfab(i,j,k,n) = fine_cell;
        });
#endif
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        MultiFab& crse_flux = *crse_fluxes[idim];
        MultiFab fine_flux(crse_flux.boxArray(), crse_flux.DistributionMap(),
                           numcomp, 0);
        fine_flux.setVal(0.0);

        Orientation lo_face(idim, Orientation::low);
        Orientation hi_face(idim, Orientation::high);

        bndry[lo_face].copyTo(fine_flux, 0, srccomp, 0, numcomp, cperiod);
        bndry[hi_face].plusTo(fine_flux, 0, srccomp, 0, numcomp, cperiod);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (!run_on_gpu)
#endif
        for (MFIter mfi(crse_flux,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& cflx = crse_flux[mfi].array(destcomp);
            auto const& fflx = fine_flux[mfi].const_array();
            auto const& msk = cc_mask[mfi].const_array();
            if (idim == 0) {
                AMREX_HOST_DEVICE_FOR_4D(bx, numcomp, i, j, k, n,
                {
                    if ( (msk(i-1,j,k) == crse_cell && msk(i,j,k) == fine_cell) ||
                         (msk(i-1,j,k) == fine_cell && msk(i,j,k) == crse_cell) ) {
                        cflx(i,j,k,n) = scale*fflx(i,j,k,n);
                    }
                });
            } else if (idim == 1) {
                AMREX_HOST_DEVICE_FOR_4D(bx, numcomp, i, j, k, n,
                {
                    if ( (msk(i,j-1,k) == crse_cell && msk(i,j,k) == fine_cell) ||
                         (msk(i,j-1,k) == fine_cell && msk(i,j,k) == crse_cell) ) {
                        cflx(i,j,k,n) = scale*fflx(i,j,k,n);
                    }
                });
            } else if(idim == 2) {
                AMREX_HOST_DEVICE_FOR_4D(bx, numcomp, i, j, k, n,
                {
                    if ( (msk(i,j,k-1) == crse_cell && msk(i,j,k) == fine_cell) ||
                         (msk(i,j,k-1) == fine_cell && msk(i,j,k) == crse_cell) ) {
                        cflx(i,j,k,n) = scale*fflx(i,j,k,n);
                    }
                });
            }
        }
    }
}

void
FluxRegister::write (const std::string& name, std::ostream& os) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        os << ratio      << '\n';
        os << fine_level << '\n';
        os << ncomp      << '\n';
    }

    const BndryRegister* br = this;

    br->write(name,os);
}


void
FluxRegister::read (const std::string& name, std::istream& is)
{
    if (ncomp < 0) {
        amrex::Abort("FluxRegister::read: FluxRegister not defined");
    }

    IntVect ratio_in;
    int fine_level_in;
    int ncomp_in;

    is >> ratio_in;
    is >> fine_level_in;
    is >> ncomp_in;

    if (ratio_in != ratio || fine_level_in != fine_level || ncomp_in != ncomp) {
        amrex::Abort("FluxRegister::read: predefined FluxRegister does not match the one in istream");
    }

    BndryRegister* br = this;

    br->read(name,is);
}

}
