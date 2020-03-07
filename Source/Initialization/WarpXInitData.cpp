/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 * Axel Huebl, Burlen Loring, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Filter/BilinearFilter.H"
#include "Filter/NCIGodfreyFilter.H"
#include "Parser/GpuParser.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif


using namespace amrex;

void
WarpX::InitData ()
{
    WARPX_PROFILE("WarpX::InitData()");

    if (restart_chkfile.empty())
    {
        ComputeDt();
        InitFromScratch();
    }
    else
    {
        InitFromCheckpoint();
        if (is_synchronized) {
            ComputeDt();
        }
        PostRestart();
    }

    ComputePMLFactors();

    if (WarpX::use_fdtd_nci_corr) {
        WarpX::InitNCICorrector();
    }

    if (WarpX::use_filter) {
        WarpX::InitFilter();
    }

    BuildBufferMasks();

    InitDiagnostics();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nGrids Summary:\n";
        printGridSummary(std::cout, 0, finestLevel());
    }

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge = new amrex::AmrMeshInSituBridge;
    insitu_bridge->setEnabled(insitu_int > 0 ? 1 : 0);
    insitu_bridge->setConfig(insitu_config);
    insitu_bridge->setPinMesh(insitu_pin_mesh);
    if (insitu_bridge->initialize())
    {
        amrex::ErrorStream()
            << "WarpX::InitData : Failed to initialize the in situ bridge."
            << std::endl;

        amrex::Abort();
    }
    insitu_bridge->setFrequency(1);
#endif

    if (restart_chkfile.empty())
    {
        if (plot_int > 0)
            WritePlotFile();

        if (openpmd_int > 0)
            WriteOpenPMDFile();

        if (check_int > 0)
            WriteCheckPointFile();

        if ((insitu_int > 0) && (insitu_start == 0))
            UpdateInSitu();

        // Write reduced diagnostics before the first iteration.
        if (reduced_diags->m_plot_rd != 0)
        {
            reduced_diags->ComputeDiags(-1);
            reduced_diags->WriteToFile(-1);
        }
    }
}

void
WarpX::InitDiagnostics () {
    if (do_back_transformed_diagnostics) {
        const Real* current_lo = geom[0].ProbLo();
        const Real* current_hi = geom[0].ProbHi();
        Real dt_boost = dt[0];
        // Find the positions of the lab-frame box that corresponds to the boosted-frame box at t=0
        Real zmin_lab = current_lo[moving_window_dir]/( (1.+beta_boost)*gamma_boost );
        Real zmax_lab = current_hi[moving_window_dir]/( (1.+beta_boost)*gamma_boost );
        myBFD.reset(new BackTransformedDiagnostic(zmin_lab,
                                               zmax_lab,
                                               moving_window_v, dt_snapshots_lab,
                                               num_snapshots_lab,
                                               dt_slice_snapshots_lab,
                                               num_slice_snapshots_lab,
                                               gamma_boost, t_new[0], dt_boost,
                                               moving_window_dir, geom[0],
                                               slice_realbox,
                                               particle_slice_width_lab));
    }
}

void
WarpX::InitFromScratch ()
{
    const Real time = 0.0;

    AmrCore::InitFromScratch(time);  // This will call MakeNewLevelFromScratch

    mypc->AllocData();
    mypc->InitData();

    // Loop through species and calculate their space-charge field
    bool const reset_fields = false; // Do not erase previous user-specified values on the grid
    ComputeSpaceChargeField(reset_fields);

    InitPML();
}

void
WarpX::InitPML ()
{
    if (do_pml)
    {
        amrex::IntVect do_pml_Lo_corrected = do_pml_Lo;

#ifdef WARPX_DIM_RZ
        do_pml_Lo_corrected[0] = 0; // no PML at r=0, in cylindrical geometry
#endif
        pml[0].reset(new PML(boxArray(0), DistributionMap(0), &Geom(0), nullptr,
                             pml_ncell, pml_delta, 0,
#ifdef WARPX_USE_PSATD
                             dt[0], nox_fft, noy_fft, noz_fft, do_nodal,
#endif
                             do_dive_cleaning, do_moving_window,
                             pml_has_particles, do_pml_in_domain,
                             do_pml_Lo_corrected, do_pml_Hi));
        for (int lev = 1; lev <= finest_level; ++lev)
        {
            amrex::IntVect do_pml_Lo_MR = amrex::IntVect::TheUnitVector();
#ifdef WARPX_DIM_RZ
            //In cylindrical geometry, if the edge of the patch is at r=0, do not add PML
            if ((max_level > 0) && (fine_tag_lo[0]==0.)) {
                do_pml_Lo_MR[0] = 0;
            }
#endif
            pml[lev].reset(new PML(boxArray(lev), DistributionMap(lev),
                                   &Geom(lev), &Geom(lev-1),
                                   pml_ncell, pml_delta, refRatio(lev-1)[0],
#ifdef WARPX_USE_PSATD
                                   dt[lev], nox_fft, noy_fft, noz_fft, do_nodal,
#endif
                                   do_dive_cleaning, do_moving_window,
                                   pml_has_particles, do_pml_in_domain,
                                   do_pml_Lo_MR, amrex::IntVect::TheUnitVector()));
        }
    }
}

void
WarpX::ComputePMLFactors ()
{
    if (do_pml)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            pml[lev]->ComputePMLFactors(dt[lev]);
        }
    }
}

void
WarpX::InitNCICorrector ()
{
    if (WarpX::use_fdtd_nci_corr)
    {
        for (int lev = 0; lev <= max_level; ++lev)
        {
            const Geometry& gm = Geom(lev);
            const Real* dx = gm.CellSize();
            amrex::Real dz, cdtodz;
            if (AMREX_SPACEDIM == 3){
                dz = dx[2];
            }else{
                dz = dx[1];
            }
            cdtodz = PhysConst::c * dt[lev] / dz;

            // Initialize Godfrey filters
            // Same filter for fields Ex, Ey and Bz
            const bool nodal_gather = (l_lower_order_in_v == 0);
            nci_godfrey_filter_exeybz[lev].reset( new NCIGodfreyFilter(godfrey_coeff_set::Ex_Ey_Bz, cdtodz, nodal_gather) );
            // Same filter for fields Bx, By and Ez
            nci_godfrey_filter_bxbyez[lev].reset( new NCIGodfreyFilter(godfrey_coeff_set::Bx_By_Ez, cdtodz, nodal_gather) );
            // Compute Godfrey filters stencils
            nci_godfrey_filter_exeybz[lev]->ComputeStencils();
            nci_godfrey_filter_bxbyez[lev]->ComputeStencils();
        }
    }
}

void
WarpX::InitFilter (){
    if (WarpX::use_filter){
        WarpX::bilinear_filter.npass_each_dir = WarpX::filter_npass_each_dir;
        WarpX::bilinear_filter.ComputeStencils();
    }
}

void
WarpX::PostRestart ()
{
#ifdef WARPX_USE_PSATD
    amrex::Abort("WarpX::PostRestart: TODO for PSATD");
#endif
    mypc->PostRestart();
}


void
WarpX::InitLevelData (int lev, Real /*time*/)
{

    ParmParse pp("warpx");

    // default values of E_external_grid and B_external_grid
    // are used to set the E and B field when "constant" or
    // "parser" is not explicitly used in the input.
    pp.query("B_ext_grid_init_style", B_ext_grid_s);
    std::transform(B_ext_grid_s.begin(),
                   B_ext_grid_s.end(),
                   B_ext_grid_s.begin(),
                   ::tolower);

    pp.query("E_ext_grid_init_style", E_ext_grid_s);
    std::transform(E_ext_grid_s.begin(),
                   E_ext_grid_s.end(),
                   E_ext_grid_s.begin(),
                   ::tolower);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    if (B_ext_grid_s == "constant")
        pp.getarr("B_external_grid", B_external_grid);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    if (E_ext_grid_s == "constant")
        pp.getarr("E_external_grid", E_external_grid);

    for (int i = 0; i < 3; ++i) {
        current_fp[lev][i]->setVal(0.0);
        if (lev > 0)
           current_cp[lev][i]->setVal(0.0);

        if (B_ext_grid_s == "constant" || B_ext_grid_s == "default") {
           Bfield_fp[lev][i]->setVal(B_external_grid[i]);
           if (lev > 0) {
              Bfield_aux[lev][i]->setVal(B_external_grid[i]);
              Bfield_cp[lev][i]->setVal(B_external_grid[i]);
           }
        }
        if (E_ext_grid_s == "constant" || E_ext_grid_s == "default") {
           Efield_fp[lev][i]->setVal(E_external_grid[i]);
           if (lev > 0) {
              Efield_aux[lev][i]->setVal(E_external_grid[i]);
              Efield_cp[lev][i]->setVal(E_external_grid[i]);
           }
        }
    }

    // if the input string for the B-field is "parse_b_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (B_ext_grid_s == "parse_b_ext_grid_function") {

#ifdef WARPX_DIM_RZ
       amrex::Abort("E and B parser for external fields does not work with RZ -- TO DO");
#endif
       Store_parserString(pp, "Bx_external_grid_function(x,y,z)",
                                                    str_Bx_ext_grid_function);
       Store_parserString(pp, "By_external_grid_function(x,y,z)",
                                                    str_By_ext_grid_function);
       Store_parserString(pp, "Bz_external_grid_function(x,y,z)",
                                                    str_Bz_ext_grid_function);

       Bxfield_parser.reset(new ParserWrapper<3>(
                                makeParser(str_Bx_ext_grid_function,{"x","y","z"})));
       Byfield_parser.reset(new ParserWrapper<3>(
                                makeParser(str_By_ext_grid_function,{"x","y","z"})));
       Bzfield_parser.reset(new ParserWrapper<3>(
                                makeParser(str_Bz_ext_grid_function,{"x","y","z"})));

       // Initialize Bfield_fp with external function
       InitializeExternalFieldsOnGridUsingParser(Bfield_fp[lev][0].get(),
                                                 Bfield_fp[lev][1].get(),
                                                 Bfield_fp[lev][2].get(),
                                                 Bxfield_parser.get(),
                                                 Byfield_parser.get(),
                                                 Bzfield_parser.get(),
                                                 Bx_nodal_flag, By_nodal_flag,
                                                 Bz_nodal_flag, lev);
       if (lev > 0) {
          InitializeExternalFieldsOnGridUsingParser(Bfield_aux[lev][0].get(),
                                                    Bfield_aux[lev][1].get(),
                                                    Bfield_aux[lev][2].get(),
                                                    Bxfield_parser.get(),
                                                    Byfield_parser.get(),
                                                    Bzfield_parser.get(),
                                                    Bx_nodal_flag, By_nodal_flag,
                                                    Bz_nodal_flag, lev);

          InitializeExternalFieldsOnGridUsingParser(Bfield_cp[lev][0].get(),
                                                    Bfield_cp[lev][1].get(),
                                                    Bfield_cp[lev][2].get(),
                                                    Bxfield_parser.get(),
                                                    Byfield_parser.get(),
                                                    Bzfield_parser.get(),
                                                    Bx_nodal_flag, By_nodal_flag,
                                                    Bz_nodal_flag, lev);
       }
    }

    // if the input string for the E-field is "parse_e_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (E_ext_grid_s == "parse_e_ext_grid_function") {

#ifdef WARPX_DIM_RZ
       amrex::Abort("E and B parser for external fields does not work with RZ -- TO DO");
#endif
       Store_parserString(pp, "Ex_external_grid_function(x,y,z)",
                                                    str_Ex_ext_grid_function);
       Store_parserString(pp, "Ey_external_grid_function(x,y,z)",
                                                    str_Ey_ext_grid_function);
       Store_parserString(pp, "Ez_external_grid_function(x,y,z)",
                                                    str_Ez_ext_grid_function);

       Exfield_parser.reset(new ParserWrapper<3>(
                                makeParser(str_Ex_ext_grid_function,{"x","y","z"})));
       Eyfield_parser.reset(new ParserWrapper<3>(
                                makeParser(str_Ey_ext_grid_function,{"x","y","z"})));
       Ezfield_parser.reset(new ParserWrapper<3>(
                                makeParser(str_Ez_ext_grid_function,{"x","y","z"})));

       // Initialize Efield_fp with external function
       InitializeExternalFieldsOnGridUsingParser(Efield_fp[lev][0].get(),
                                                 Efield_fp[lev][1].get(),
                                                 Efield_fp[lev][2].get(),
                                                 Exfield_parser.get(),
                                                 Eyfield_parser.get(),
                                                 Ezfield_parser.get(),
                                                 Ex_nodal_flag, Ey_nodal_flag,
                                                 Ez_nodal_flag, lev);
       if (lev > 0) {
          InitializeExternalFieldsOnGridUsingParser(Efield_aux[lev][0].get(),
                                                    Efield_aux[lev][1].get(),
                                                    Efield_aux[lev][2].get(),
                                                    Exfield_parser.get(),
                                                    Eyfield_parser.get(),
                                                    Ezfield_parser.get(),
                                                    Ex_nodal_flag, Ey_nodal_flag,
                                                    Ez_nodal_flag, lev);

          InitializeExternalFieldsOnGridUsingParser(Efield_cp[lev][0].get(),
                                                    Efield_cp[lev][1].get(),
                                                    Efield_cp[lev][2].get(),
                                                    Exfield_parser.get(),
                                                    Eyfield_parser.get(),
                                                    Ezfield_parser.get(),
                                                    Ex_nodal_flag, Ey_nodal_flag,
                                                    Ez_nodal_flag, lev);
       }
    }

    if (F_fp[lev]) {
        F_fp[lev]->setVal(0.0);
    }

    if (rho_fp[lev]) {
        rho_fp[lev]->setVal(0.0);
    }

    if (F_cp[lev]) {
        F_cp[lev]->setVal(0.0);
    }

    if (rho_cp[lev]) {
        rho_cp[lev]->setVal(0.0);
    }

    if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
        if (costs[lev]) {
            costs[lev]->setVal(0.0);
        }
    } else if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic) {
        if (costs_heuristic[lev]) {
            std::fill((*costs_heuristic[lev]).begin(),
                      (*costs_heuristic[lev]).end(),
                      0.0);
        }
    }
}

#ifdef WARPX_USE_PSATD_HYBRID

void
WarpX::InitLevelDataFFT (int lev, Real time)
{

    Efield_fp_fft[lev][0]->setVal(0.0);
    Efield_fp_fft[lev][1]->setVal(0.0);
    Efield_fp_fft[lev][2]->setVal(0.0);
    Bfield_fp_fft[lev][0]->setVal(0.0);
    Bfield_fp_fft[lev][1]->setVal(0.0);
    Bfield_fp_fft[lev][2]->setVal(0.0);
    current_fp_fft[lev][0]->setVal(0.0);
    current_fp_fft[lev][1]->setVal(0.0);
    current_fp_fft[lev][2]->setVal(0.0);
    rho_fp_fft[lev]->setVal(0.0);

    if (lev > 0)
    {
        Efield_cp_fft[lev][0]->setVal(0.0);
        Efield_cp_fft[lev][1]->setVal(0.0);
        Efield_cp_fft[lev][2]->setVal(0.0);
        Bfield_cp_fft[lev][0]->setVal(0.0);
        Bfield_cp_fft[lev][1]->setVal(0.0);
        Bfield_cp_fft[lev][2]->setVal(0.0);
        current_cp_fft[lev][0]->setVal(0.0);
        current_cp_fft[lev][1]->setVal(0.0);
        current_cp_fft[lev][2]->setVal(0.0);
        rho_cp_fft[lev]->setVal(0.0);
    }

}

#endif

void
WarpX::InitializeExternalFieldsOnGridUsingParser (
       MultiFab *mfx, MultiFab *mfy, MultiFab *mfz,
       ParserWrapper<3> *xfield_parser, ParserWrapper<3> *yfield_parser,
       ParserWrapper<3> *zfield_parser, IntVect x_nodal_flag,
       IntVect y_nodal_flag, IntVect z_nodal_flag,
       const int lev)
{

    const auto dx_lev = geom[lev].CellSizeArray();
    const RealBox& real_box = geom[lev].ProbDomain();
    for ( MFIter mfi(*mfx, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& tbx = mfi.tilebox(x_nodal_flag);
       const Box& tby = mfi.tilebox(y_nodal_flag);
       const Box& tbz = mfi.tilebox(z_nodal_flag);

       auto const& mfxfab = mfx->array(mfi);
       auto const& mfyfab = mfy->array(mfi);
       auto const& mfzfab = mfz->array(mfi);

       auto const& mfx_IndexType = (*mfx).ixType();
       auto const& mfy_IndexType = (*mfy).ixType();
       auto const& mfz_IndexType = (*mfz).ixType();

       // Initialize IntVect based on the index type of multiFab
       // 0 if cell-centered, 1 if node-centered.
       IntVect mfx_type(AMREX_D_DECL(0,0,0));
       IntVect mfy_type(AMREX_D_DECL(0,0,0));
       IntVect mfz_type(AMREX_D_DECL(0,0,0));

       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
           mfx_type[idim] = mfx_IndexType.nodeCentered(idim);
           mfy_type[idim] = mfy_IndexType.nodeCentered(idim);
           mfz_type[idim] = mfz_IndexType.nodeCentered(idim);
       }

       amrex::ParallelFor (tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Shift required in the x-, y-, or z- position
                // depending on the index type of the multifab
                Real fac_x = (1.0 - mfx_type[0]) * dx_lev[0]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
#if (AMREX_SPACEDIM==2)
                Real y = 0.0;
                Real fac_z = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#else
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                Real fac_z = (1.0 - mfx_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the x-component of the field.
                mfxfab(i,j,k) = (*xfield_parser)(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fac_x = (1.0 - mfy_type[0]) * dx_lev[0]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
#if (AMREX_SPACEDIM==2)
                Real y = 0.0;
                Real fac_z = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif (AMREX_SPACEDIM==3)
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                Real fac_z = (1.0 - mfx_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the y-component of the field.
                mfyfab(i,j,k)  = (*yfield_parser)(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fac_x = (1.0 - mfz_type[0]) * dx_lev[0]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
#if (AMREX_SPACEDIM==2)
                Real y = 0.0;
                Real fac_z = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif (AMREX_SPACEDIM==3)
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                Real fac_z = (1.0 - mfz_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the z-component of the field.
                mfzfab(i,j,k) = (*zfield_parser)(x,y,z);
            }
        );
    }

}
