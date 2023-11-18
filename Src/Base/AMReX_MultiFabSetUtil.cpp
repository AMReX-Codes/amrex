// #include <AMReX_MultiFabSetUtil.H>

// namespace amrex
// {

// // void average_down (const MultiFabSet& fine, MultiFabSet& crse, const IntVect& ratio,
// //                    const IntVect& ngcrse)
// // {
// //     for (int iSet = 0; iSet < fine.nSet(); iSet++) {
// //         amrex::average_down_same_type(fine[iSet], crse[iSet], ratio, ngcrse);
// //     }
// // }

// void average_down_same_type (const MultiFab& fine, MultiFab& crse, const IntVect& ratio,
//                              const IntVect& ngcrse)
// {
//     BL_PROFILE("average_down_same_type");

//     AMREX_ASSERT(crse.nComp() == fine.nComp());
//     AMREX_ASSERT(fine.ixType() == crse.ixType());

//     int ncomp = fine.nComp();
//     const auto type = fine.ixType();

// // #ifdef AMREX_USE_GPU
//     if (Gpu::inLaunchRegion() && crse.isFusingCandidate()) {
//         const auto& crsema = crse.arrays();
//         const auto& finema = fine.const_arrays();
//         ParallelFor(crse, ngcrse, ncomp,
//             [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
//             {
//                 amrex_avgdown_general(i,j,k,n,crsema[box_no],finema[box_no],ratio,type);
//             });
//         if (!Gpu::inNoSyncRegion()) {
//             Gpu::streamSynchronize();
//         }
//     } else
// // #endif
//     {
// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//         for (MFIter mfi(crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             //  NOTE: The tilebox is defined at the coarse level.
//             const Box& bx = mfi.growntilebox(ngcrse);
//             const auto& crsearr = crse.array(mfi);
//             const auto& finearr = fine.const_array(mfi);
//             AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
//             {
//                 amrex_avgdown_general(i,j,k,n,crsearr,finearr,ratio,type);
//             });
//         }
//     }
// }

// // void ComputeIndexRange

// void amrex_avgdown_general(int i, int j, int k, int n,
//                            const Array4<Real>& crse, const Array4<const Real>& fine,
//                            const IntVect& ratio, const IndexType ixType)
// {
//     IntVect ijk(AMREX_D_DECL(i,j,k));
//     Dim3 indexFineMin{i,j,k};
//     Dim3 indexChangeMax{1,1,1};
//     int denom = 1;

//     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
//         indexFineMin[dir] = ratio[dir]*ijk[dir];
//         if (ixType.cellCentered(dir)) {
//             indexChangeMax[dir] = ratio[dir];
//             denom *= ratio[dir];
//         }
//     }

//     Real c = Real(0.0);
//     for (int kk = indexFineMin[2]; kk < indexFineMin[2] + indexChangeMax[2]; ++kk) {
//     for (int jj = indexFineMin[1]; jj < indexFineMin[1] + indexChangeMax[1]; ++jj) {
//     for (int ii = indexFineMin[0]; ii < indexFineMin[1] + indexChangeMax[0]; ++ii) {
//         c += fine(ii,jj,kk,n);
//     }}}

//     crse(i,j,k,n) = c / Real(denom);

// }

// } // namespace amrex