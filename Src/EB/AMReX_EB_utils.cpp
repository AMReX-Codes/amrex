#include <AMReX_EB_F.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB_utils.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>


namespace amrex {

    void FillEBNormals(MultiFab & normals, const EBFArrayBoxFactory & eb_factory,
                       const Geometry & geom) {

        BL_PROFILE("amrex::FillEBNormals()")

        BoxArray ba = normals.boxArray();
        DistributionMapping dm = normals.DistributionMap();
        int n_grow = normals.nGrow();

        // Dummy array for MFIter
        MultiFab dummy(ba, dm, 1, n_grow, MFInfo(), eb_factory);
        // Area fraction data
        std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory.getAreaFrac();

        const auto & flags = eb_factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(dummy, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.growntilebox();
            const int * lo = tile_box.loVect();
            const int * hi = tile_box.hiVect();

            const auto & flag = flags[mfi];

            if (flag.getType(tile_box) == FabType::singlevalued) {
                // Target for compute_normals(...)
                auto & norm_tile = normals[mfi];
                // Area fractions in x, y, and z directions
                const auto & af_x_tile = (* areafrac[0])[mfi];
                const auto & af_y_tile = (* areafrac[1])[mfi];
                const auto & af_z_tile = (* areafrac[2])[mfi];

                amrex_eb_compute_normals(lo, hi,
                                         BL_TO_FORTRAN_3D(flag),
                                         BL_TO_FORTRAN_3D(norm_tile),
                                         BL_TO_FORTRAN_3D(af_x_tile),
                                         BL_TO_FORTRAN_3D(af_y_tile),
                                         BL_TO_FORTRAN_3D(af_z_tile)  );
            }
        }

        normals.FillBoundary(geom.periodicity());
    }

}
