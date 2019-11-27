
#include <WarpX.H>
#include <WarpX_f.H>
#include <BilinearFilter.H>
#include <NCIGodfreyFilter.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif
#include <GpuParser.H>


using namespace amrex;

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

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

        if (check_int > 0)
            WriteCheckPointFile();

        if ((insitu_int > 0) && (insitu_start == 0))
            UpdateInSitu();
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
    for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
        WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
        if (species.initialize_self_fields) {
            InitSpaceChargeField(species);
        }
    }

    InitPML();

#ifdef WARPX_DO_ELECTROSTATIC
    if (do_electrostatic) {
        getLevelMasks(masks);

        // the plus one is to convert from num_cells to num_nodes
        getLevelMasks(gather_masks, n_buffer + 1);
    }
#endif // WARPX_DO_ELECTROSTATIC
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


namespace {
WarpXParser makeParser (std::string const& parse_function)
{
    std::cout << " in make parser " << parse_function << std::endl;
    WarpXParser parser(parse_function);
    parser.registerVariables({"x","y","z"});
    ParmParse pp("my_constants");
    std::set<std::string> symbols = parser.symbols();
    symbols.erase("x");
    symbols.erase("y");
    symbols.erase("z");
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;
        if (pp.query(it->c_str(), v)) {
           parser.setConstant(*it, v);
           it = symbols.erase(it);
        } else {
           ++it;
        }
    }
    for (auto const& s : symbols) {
        amrex::Abort(" ExternalEBFieldOnGrid::makeParser::Unknown symbol "+s);
    }
    return parser;
}
}




/* \brief
 *  This function initializes E, B, rho, and F, at all the levels
 *  of the multifab. rho and F are initialized with 0.
 *  The E and B fields are initialized using user-defined inputs.
 *  The initialization type is set using "B_ext_grid_init_style"
 *  and "E_ext_grid_init_style". The initialization style is set to "default"
 *  if not explicitly defined by the user, and the E and B fields are
 *  initialized with E_external_grid and B_external_grid, respectively, each with
 *  a default value of 0.
 *  If the initialization type for the E and B field is "constant",
 *  then, the E and B fields at all the levels are initialized with
 *  user-defined values for E_external_grid and B_external_grid.
 *  If the initialization type for B-field is set to
 *  "parse_B_ext_grid_function", then, the parser is used to read
 *  Bx_external_grid_function(x,y,z), By_external_grid_function(x,y,z),
 *  and Bz_external_grid_function(x,y,z).
 *  Similarly, if the E-field initialization type is set to
 *  "parse_E_ext_grid_function", then, the parser is used to read
 *  Ex_external_grid_function(x,y,z), Ey_external_grid_function(x,y,z),
 *  and Ex_external_grid_function(x,y,z). The parser for the E and B
 *  initialization assumes that the function has three independent
 *  variables, at max, namely, x, y, z. However, any number of constants
 *  can be used in the function used to define the E and B fields on the grid.
 */

void
WarpX::InitLevelData (int lev, Real time)
{

    ParmParse pp("warpx");

    std::string B_ext_grid_s;
    std::string E_ext_grid_s;
    // default values of E_external_grid and B_external_grid
    // are used to set the E and B field when "constant" or
    // "parser" is not explicitly used in the input.
    B_ext_grid_s = "default";
    E_ext_grid_s = "default";

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
        if (B_ext_grid_s == "constant" || B_ext_grid_s == "default")
           Bfield_fp[lev][i]->setVal(B_external_grid[i]);
        if (E_ext_grid_s == "constant" || E_ext_grid_s == "default")
           Efield_fp[lev][i]->setVal(E_external_grid[i]);
    }

    if (B_ext_grid_s == "parse_b_ext_grid_function") {

       std::vector<std::string> f;
       // Parse Bx_external_grid_function
       pp.getarr("Bx_external_grid_function(x,y,z)", f);
       str_Bx_ext_grid_function.clear();
       for (auto const& s : f) {
           str_Bx_ext_grid_function += s;
       }
       f.clear();

       // Parse By_external_grid_function
       pp.getarr("By_external_grid_function(x,y,z)", f);
       str_By_ext_grid_function.clear();
       for (auto const& s : f) {
            str_By_ext_grid_function += s;
       }
       f.clear();

       // Parse Bz_external_grid_function
       pp.getarr("Bz_external_grid_function(x,y,z)", f);
       str_Bz_ext_grid_function.clear();
       for (auto const& s : f) {
           str_Bz_ext_grid_function += s;
       }
       f.clear();

       // Initialize Bfield_fp with external function
       MultiFab *Bx, *By, *Bz;
       Bx = Bfield_fp[lev][0].get();
       By = Bfield_fp[lev][1].get();
       Bz = Bfield_fp[lev][2].get();

       bool B_flag = 1;
       InitializeExternalFieldsOnGridUsingParser(Bx, By, Bz, lev, B_flag);
    }

    if (E_ext_grid_s == "parse_e_ext_grid_function") {

       std::vector<std::string> f;
       // Parse Ex_external_grid_function
       pp.getarr("Ex_external_grid_function(x,y,z)", f);
       str_Ex_ext_grid_function.clear();
       for (auto const& s : f) {
           str_Ex_ext_grid_function += s;
       }
       f.clear();

       // Parse Ey_external_grid_function
       pp.getarr("Ey_external_grid_function(x,y,z)", f);
       str_Ey_ext_grid_function.clear();
       for (auto const& s : f) {
           str_Ey_ext_grid_function += s;
       }
       f.clear();

       // Parse Ez_external_grid_function
       pp.getarr("Ez_external_grid_function(x,y,z)", f);
       str_Ez_ext_grid_function.clear();
       for (auto const& s : f) {
           str_Ez_ext_grid_function += s;
       }
       f.clear();

       // Initialize Efield_fp with external function
       MultiFab *Ex, *Ey, *Ez;
       Ex = Efield_fp[lev][0].get();
       Ey = Efield_fp[lev][1].get();
       Ez = Efield_fp[lev][2].get();

       bool B_flag = 0;
       InitializeExternalFieldsOnGridUsingParser(Ex, Ey, Ez, lev, B_flag);

    }

    if (lev > 0) {
        for (int i = 0; i < 3; ++i) {
            current_cp[lev][i]->setVal(0.0);
            if (B_ext_grid_s == "constant" || B_ext_grid_s == "default") {
               Bfield_aux[lev][i]->setVal(B_external_grid[i]);
               Bfield_cp[lev][i]->setVal(B_external_grid[i]);
            }
            if (E_ext_grid_s == "constant" || E_ext_grid_s == " default") {
               Efield_aux[lev][i]->setVal(E_external_grid[i]);
               Efield_cp[lev][i]->setVal(E_external_grid[i]);
            }
        }
        if (B_ext_grid_s == "parse_b_ext_grid_function") {

            MultiFab *Bx_aux, *By_aux, *Bz_aux;
            Bx_aux = Bfield_aux[lev][0].get();
            By_aux = Bfield_aux[lev][1].get();
            Bz_aux = Bfield_aux[lev][2].get();

            // Setting b_flag to 1 since we are initializing
            // B_external on the grid.
            bool B_flag = 1;
            InitializeExternalFieldsOnGridUsingParser(Bx_aux, By_aux,
                                                      Bz_aux, lev, B_flag);

            MultiFab *Bx_cp, *By_cp, *Bz_cp;
            Bx_cp = Bfield_cp[lev][0].get();
            By_cp = Bfield_cp[lev][1].get();
            Bz_cp = Bfield_cp[lev][2].get();

            InitializeExternalFieldsOnGridUsingParser(Bx_cp, By_cp,
                                                      Bz_cp, lev, B_flag);

        }
        if (E_ext_grid_s == "parse_e_ext_grid_function") {

            MultiFab *Ex_aux, *Ey_aux, *Ez_aux;
            Ex_aux = Efield_aux[lev][0].get();
            Ey_aux = Efield_aux[lev][1].get();
            Ez_aux = Efield_aux[lev][2].get();
            // Setting b_flag to zero since we are initializing
            // E_external on the grid here.
            bool B_flag = 0;
            InitializeExternalFieldsOnGridUsingParser(Ex_aux, Ey_aux,
                                                      Ez_aux, lev, B_flag);

            MultiFab *Ex_cp, *Ey_cp, *Ez_cp;
            Ex_cp = Efield_cp[lev][0].get();
            Ey_cp = Efield_cp[lev][1].get();
            Ez_cp = Efield_cp[lev][2].get();

            InitializeExternalFieldsOnGridUsingParser(Ex_cp, Ey_cp,
                                                      Ez_cp, lev, B_flag);

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

    if (costs[lev]) {
        costs[lev]->setVal(0.0);
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

/* \brief
 * This function initializes the E and B fields on each level
 * using the parser and the user-defined function for the external fields.
 * Depending on the bool value of the B_flag, the subroutine will parse
 * the Bx_/By_/Bz_external_grid_function (if B_flag==1)
 * or parse the Ex_/Ey_/Ez_external_grid_function (if B_flag==0).
 * Then, the B or E multifab is initialized based on the (x,y,z) position
 * on the staggered yee-grid or cell-centered grid.
 */

void
WarpX::InitializeExternalFieldsOnGridUsingParser (
       MultiFab *mfx, MultiFab *mfy, MultiFab *mfz,
       const int lev, const bool B_flag)
{
    std::unique_ptr<ParserWrapper> xfield_parsewrap;
    std::unique_ptr<ParserWrapper> yfield_parsewrap;
    std::unique_ptr<ParserWrapper> zfield_parsewrap;

    if (B_flag == 1) {
        xfield_parsewrap.reset(new ParserWrapper
                          (makeParser(str_Bx_ext_grid_function)));
        yfield_parsewrap.reset(new ParserWrapper
                          (makeParser(str_By_ext_grid_function)));
        zfield_parsewrap.reset(new ParserWrapper
                          (makeParser(str_Bz_ext_grid_function)));
    } else {
        xfield_parsewrap.reset(new ParserWrapper
                          (makeParser(str_Ex_ext_grid_function)));
        yfield_parsewrap.reset(new ParserWrapper
                          (makeParser(str_Ey_ext_grid_function)));
        zfield_parsewrap.reset(new ParserWrapper
                           (makeParser(str_Ez_ext_grid_function)));
    }

    ParserWrapper *xfield_wrap = xfield_parsewrap.get();
    ParserWrapper *yfield_wrap = yfield_parsewrap.get();
    ParserWrapper *zfield_wrap = zfield_parsewrap.get();

    const auto dx_lev = geom[lev].CellSizeArray();
    const RealBox& real_box = geom[lev].ProbDomain();
    for ( MFIter mfi(*mfx, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       IntVect x_nodal_flag, y_nodal_flag, z_nodal_flag;
       if (B_flag == 1) {
          x_nodal_flag = Bx_nodal_flag;
          y_nodal_flag = By_nodal_flag;
          z_nodal_flag = Bz_nodal_flag;
       } else {
          x_nodal_flag = Ex_nodal_flag;
          y_nodal_flag = Ey_nodal_flag;
          z_nodal_flag = Ez_nodal_flag;
       }
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
                Real fac_y = (1.0 - mfx_type[1]) * dx_lev[1]*0.5;

                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
#if (AMREX_SPACEDIM==2)
                Real z = 0.0;
#else
                Real fac_z = (1.0 - mfx_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the x-component of the field.
                mfxfab(i,j,k) = xfield_wrap->getField(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fac_x = (1.0 - mfy_type[0]) * dx_lev[0]*0.5;
                Real fac_y = (1.0 - mfy_type[1]) * dx_lev[1]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
#if (AMREX_SPACEDIM==2)
                Real z = 0.0;
#elif (AMREX_SPACEDIM==3)
                Real fac_z = (1.0 - mfy_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the y-component of the field.
                mfyfab(i,j,k)  = yfield_wrap->getField(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fac_x = (1.0 - mfz_type[0]) * dx_lev[0]*0.5;
                Real fac_y = (1.0 - mfz_type[1]) * dx_lev[1]*0.5;
                Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
#if (AMREX_SPACEDIM==2)
                Real z = 0.0;
#elif (AMREX_SPACEDIM==3)
                Real fac_z = (1.0 - mfz_type[2]) * dx_lev[2]*0.5;
                Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the z-component of the field.
                mfzfab(i,j,k) = zfield_wrap->getField(x,y,z);
            },
            amrex::Gpu::numThreadsPerBlockParallelFor() * sizeof(double) * 3
        );

    }

}
