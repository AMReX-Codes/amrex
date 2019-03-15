
#include <WarpX.H>
#include <FieldIO.H>

using namespace amrex;

void
PackPlotDataPtrs (Vector<const MultiFab*>& pmf,
                  const std::array<std::unique_ptr<MultiFab>,3>& data)
{
    BL_ASSERT(pmf.size() == AMREX_SPACEDIM);
#if (AMREX_SPACEDIM == 3)
    pmf[0] = data[0].get();
    pmf[1] = data[1].get();
    pmf[2] = data[2].get();
#elif (AMREX_SPACEDIM == 2)
    pmf[0] = data[0].get();
    pmf[1] = data[2].get();
#endif
}

/** \brief Takes an array of 3 MultiFab `vector_field`
 * (representing the x, y, z components of a vector),
 * averages it to the cell center, and stores the
 * resulting MultiFab in mf_avg (in the components dcomp to dcomp+2)
 */
void
AverageAndPackVectorField( MultiFab& mf_avg,
                           const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
                           const int dcomp, const int ngrow )
{
    // The object below is temporary, and is needed because
    // `average_edge_to_cellcenter` requires fields to be passed as Vector
    Vector<const MultiFab*> srcmf(AMREX_SPACEDIM);

    // Check the type of staggering of the 3-component `vector_field`
    // and average accordingly:
    // - Fully cell-centered field (no average needed; simply copy)
    if ( vector_field[0]->is_cell_centered() ){

        MultiFab::Copy( mf_avg, *vector_field[0], 0, dcomp  , 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[1], 0, dcomp+1, 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[2], 0, dcomp+2, 1, ngrow);

        // - Fully nodal
    } else if ( vector_field[0]->is_nodal() ){

        amrex::average_node_to_cellcenter( mf_avg, dcomp  ,
                                          *vector_field[0], 0, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                          *vector_field[1], 0, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+2,
                                          *vector_field[2], 0, 1, ngrow);

        // - Face centered, in the same way as B on a Yee grid
    } else if ( vector_field[0]->is_nodal(0) ){

        PackPlotDataPtrs(srcmf, vector_field);
        amrex::average_face_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
#if (AMREX_SPACEDIM == 2)
        MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[1], 0, dcomp+1, 1, ngrow);
#endif

        // - Edge centered, in the same way as E on a Yee grid
    } else if ( !vector_field[0]->is_nodal(0) ){

        PackPlotDataPtrs(srcmf, vector_field);
        amrex::average_edge_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
#if (AMREX_SPACEDIM == 2)
        MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                          *vector_field[1], 0, 1, ngrow);
#endif

    } else {
        amrex::Abort("Unknown staggering.");
    }
}

/** \brief Take a MultiFab `scalar_field`
 * averages it to the cell center, and stores the
 * resulting MultiFab in mf_avg (in the components dcomp)
 */
void
AverageAndPackScalarField( MultiFab& mf_avg,
                           const MultiFab & scalar_field,
                           const int dcomp, const int ngrow )
{
    // Check the type of staggering of the 3-component `vector_field`
    // and average accordingly:
    // - Fully cell-centered field (no average needed; simply copy)
    if ( scalar_field.is_cell_centered() ){

        MultiFab::Copy( mf_avg, scalar_field, 0, dcomp, 1, ngrow);

        // - Fully nodal
    } else if ( scalar_field.is_nodal() ){

        amrex::average_node_to_cellcenter( mf_avg, dcomp, scalar_field, 0, 1, ngrow);

    } else {
        amrex::Abort("Unknown staggering.");
    }
}

/** \brief Write the different fields that are meant for output,
 * into the vector of MultiFab `mf_avg` (one MultiFab per level)
 * after averaging them to the cell centers.
 */
void
WarpX::AverageAndPackFields ( Vector<std::string>& varnames,
    amrex::Vector<MultiFab>& mf_avg, const int ngrow) const
{
    // Count how many different fields should be written (ncomp)
    const int ncomp = 3*3
        + static_cast<int>(plot_part_per_cell)
        + static_cast<int>(plot_part_per_grid)
        + static_cast<int>(plot_part_per_proc)
        + static_cast<int>(plot_proc_number)
        + static_cast<int>(plot_divb)
        + static_cast<int>(plot_dive)
        + static_cast<int>(plot_rho)
        + static_cast<int>(plot_F)
        + static_cast<int>(plot_finepatch)*6
        + static_cast<int>(plot_crsepatch)*6
        + static_cast<int>(costs[0] != nullptr);

    // Loop over levels of refinement
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Allocate pointers to the `ncomp` fields that will be added
        mf_avg.push_back( MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

        // Go through the different fields, pack them into mf_avg[lev],
        // add the corresponding names to `varnames` and increment dcomp
        int dcomp = 0;
        AverageAndPackVectorField(mf_avg[lev], current_fp[lev], dcomp, ngrow);
        if(lev==0) for(auto name:{"jx","jy","jz"}) varnames.push_back(name);
        dcomp += 3;
        AverageAndPackVectorField(mf_avg[lev], Efield_aux[lev], dcomp, ngrow);
        if(lev==0) for(auto name:{"Ex","Ey","Ez"}) varnames.push_back(name);
        dcomp += 3;
        AverageAndPackVectorField(mf_avg[lev], Bfield_aux[lev], dcomp, ngrow);
        if(lev==0) for(auto name:{"Bx","By","Bz"}) varnames.push_back(name);
        dcomp += 3;

        if (plot_part_per_cell)
        {
            MultiFab temp_dat(grids[lev],mf_avg[lev].DistributionMap(),1,0);
            temp_dat.setVal(0);

            // MultiFab containing number of particles in each cell
            mypc->Increment(temp_dat, lev);
            AverageAndPackScalarField( mf_avg[lev], temp_dat, dcomp, ngrow );
            if(lev==0) varnames.push_back("part_per_cell");
            dcomp += 1;
        }

        if (plot_part_per_grid)
        {
            const Vector<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);
            // MultiFab containing number of particles per grid
            // (stored as constant for all cells in each grid)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(mf_avg[lev]); mfi.isValid(); ++mfi) {
                (mf_avg[lev])[mfi].setVal(static_cast<Real>(npart_in_grid[mfi.index()]),
                                           dcomp);
            }
            if(lev==0) varnames.push_back("part_per_grid");
            dcomp += 1;
        }

        if (plot_part_per_proc)
        {
            const Vector<long>& npart_in_grid = mypc->NumberOfParticlesInGrid(lev);
            // MultiFab containing number of particles per process
            // (stored as constant for all cells in each grid)
            long n_per_proc = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:n_per_proc)
#endif
            for (MFIter mfi(mf_avg[lev]); mfi.isValid(); ++mfi) {
                n_per_proc += npart_in_grid[mfi.index()];
            }
            mf_avg[lev].setVal(static_cast<Real>(n_per_proc), dcomp,1);
            if(lev==0) varnames.push_back("part_per_proc");
            dcomp += 1;
        }

        if (plot_proc_number)
        {
            // MultiFab containing the Processor ID
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(mf_avg[lev]); mfi.isValid(); ++mfi) {
                (mf_avg[lev])[mfi].setVal(static_cast<Real>(ParallelDescriptor::MyProc()),
                                           dcomp);
            }
            if(lev==0) varnames.push_back("proc_number");
            dcomp += 1;
        }

        if (plot_divb)
        {
            if (do_nodal) amrex::Abort("TODO: do_nodal && plot_divb");
            ComputeDivB(mf_avg[lev], dcomp,
                        {Bfield_aux[lev][0].get(),
                                Bfield_aux[lev][1].get(),
                                Bfield_aux[lev][2].get()},
                        WarpX::CellSize(lev)
                        );
            if(lev == 0) varnames.push_back("divB");
            dcomp += 1;
        }

        if (plot_dive)
        {
            if (do_nodal) amrex::Abort("TODO: do_nodal && plot_dive");
            const BoxArray& ba = amrex::convert(boxArray(lev),IntVect::TheUnitVector());
            MultiFab dive(ba,DistributionMap(lev),1,0);
            ComputeDivE( dive, 0,
                         {Efield_aux[lev][0].get(),
                                 Efield_aux[lev][1].get(),
                                 Efield_aux[lev][2].get()},
                         WarpX::CellSize(lev)
                         );
            AverageAndPackScalarField( mf_avg[lev], dive, dcomp, ngrow );
            if(lev == 0) varnames.push_back("divE");
            dcomp += 1;
        }

        if (plot_rho)
        {
            AverageAndPackScalarField( mf_avg[lev], *rho_fp[lev], dcomp, ngrow );
            if(lev == 0) varnames.push_back("rho");
            dcomp += 1;
        }

        if (plot_F)
        {
            AverageAndPackScalarField( mf_avg[lev], *F_fp[lev], dcomp, ngrow);
            if(lev == 0) varnames.push_back("F");
            dcomp += 1;
        }

        if (plot_finepatch)
        {
            AverageAndPackVectorField( mf_avg[lev], Efield_fp[lev], dcomp, ngrow );
            if(lev==0) for(auto name:{"Ex_fp","Ey_fp","Ez_fp"}) varnames.push_back(name);
            dcomp += 3;
            AverageAndPackVectorField( mf_avg[lev], Bfield_fp[lev], dcomp, ngrow );
            if(lev==0) for(auto name:{"Bx_fp","By_fp","Bz_fp"}) varnames.push_back(name);
            dcomp += 3;
        }

        if (plot_crsepatch)
        {
            if (lev == 0)
            {
                mf_avg[lev].setVal(0.0, dcomp, 3, ngrow);
            }
            else
            {
                if (do_nodal) amrex::Abort("TODO: do_nodal && plot_crsepatch");
                std::array<std::unique_ptr<MultiFab>, 3> E = getInterpolatedE(lev);
                AverageAndPackVectorField( mf_avg[lev], E, dcomp, ngrow );

            }
            if(lev==0) for(auto name:{"Ex_cp","Ey_cp","Ez_cp"}) varnames.push_back(name);
            dcomp += 3;

            // now the magnetic field
            if (lev == 0)
            {
                mf_avg[lev].setVal(0.0, dcomp, 3, ngrow);
            }
            else
            {
                if (do_nodal) amrex::Abort("TODO: do_nodal && plot_crsepatch");
                std::array<std::unique_ptr<MultiFab>, 3> B = getInterpolatedB(lev);
                AverageAndPackVectorField( mf_avg[lev], B, dcomp, ngrow );
            }
            if(lev==0) for(auto name:{"Bx_cp","By_cp","Bz_cp"}) varnames.push_back(name);
            dcomp += 3;
        }

        if (costs[0] != nullptr)
        {
            AverageAndPackScalarField( mf_avg[lev], *costs[lev], dcomp, ngrow );
            if(lev==0) varnames.push_back("costs");
            dcomp += 1;
        }

    BL_ASSERT(dcomp == ncomp);
    } // end loop over levels of refinement
};
