/* Copyright 2016-2020 Andrew Myers, Ann Almgren, Axel Huebl
 *                     David Grote, Jean-Luc Vay, Remi Lehe
 *                     Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <iostream>


int main(int argc, char* argv[])
{
    using namespace amrex;

#if defined(AMREX_USE_MPI)
#   if defined(_OPENMP) && defined(WARPX_USE_PSATD)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    AMREX_ALWAYS_ASSERT(provided >= MPI_THREAD_FUNNELED);
#   else
    MPI_Init(&argc, &argv);
#   endif
#endif

    amrex::Initialize(argc,argv);

    ConvertLabParamsToBoost();

    WARPX_PROFILE_VAR("main()", pmain);

    const Real strt_total = amrex::second();

    {
        WarpX warpx;

        warpx.InitData();

        warpx.Evolve();

        Real end_total = amrex::second() - strt_total;

        ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
        if (warpx.Verbose()) {
            amrex::Print() << "Total Time                     : " << end_total << '\n';
            amrex::Print() << "WarpX Version: " << WarpX::Version() << '\n';
            amrex::Print() << "PICSAR Version: " << WarpX::PicsarVersion() << '\n';
        }
    }

    WARPX_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
#if defined(AMREX_USE_MPI)
    MPI_Finalize();
#endif
}
