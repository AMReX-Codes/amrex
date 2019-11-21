
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

#include <FieldInit.H>

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






void
WarpX::InitLevelData (int lev, Real time)
{

    amrex::Print() << " at initlevel data " << lev << " time " << time << "\n";
    ParmParse pp("warpx");

    std::string B_ext_grid_s;
    // default is constant so that the default values of 
    // E_external_grid[i] and B_external_grid[i] can 
    // be used for setting the value of the E & B fields.
    B_ext_grid_s = "constant";
    
    pp.query("B_ext_grid_init_style", B_ext_grid_s);
    std::transform(B_ext_grid_s.begin(),
                   B_ext_grid_s.end(),
                   B_ext_grid_s.begin(),
                   ::tolower);
    amrex::Print() << " init style " << B_ext_grid_s << "\n";

    std::string E_ext_grid_s;
    pp.query("E_ext_grid_init_style", E_ext_grid_s);
    std::transform(E_ext_grid_s.begin(),
                   E_ext_grid_s.end(),
                   E_ext_grid_s.begin(),
                   ::tolower);

    if (B_ext_grid_s == "constant") {
        pp.getarr("B_external_grid", B_external_grid);     
    }
    if (B_ext_grid_s == "parse_b_ext_grid_function") 
    {
       amrex::Print() << " parse function " << "\n";
//       std::vector<std::string> f;
//       pp.getarr("Bx_external_grid_function(x,y,z)", f);
//       for (auto const& s : f) {
//           str_Bx_ext_grid_function += s;
//       }
//       f.clear();
//
    }

    if (E_ext_grid_s == "constant") {
        pp.getarr("E_external_grid", E_external_grid);     
    }

    for (int i = 0; i < 3; ++i) {
        current_fp[lev][i]->setVal(0.0);
        Bfield_fp[lev][i]->setVal(0.0);
        Efield_fp[lev][i]->setVal(0.0);
        if (B_ext_grid_s == "constant") {
           Bfield_fp[lev][i]->setVal(B_external_grid[i]);
        }
        if (E_ext_grid_s == "constant") {
           Efield_fp[lev][i]->setVal(E_external_grid[i]);
        }
    }
    if (B_ext_grid_s == "parse_b_ext_grid_function") {
       
       MultiFab *Bx, *By, *Bz;
       Bx = Bfield_fp[lev][0].get();
       By = Bfield_fp[lev][1].get();
       Bz = Bfield_fp[lev][2].get();
       const auto dx_lev = geom[lev].CellSizeArray();
       const RealBox& real_box = geom[lev].ProbDomain();
       amrex::Print() << " cell size " << dx_lev[0] << " " << dx_lev[1] << " \n";
          std::string str_Bx_function;
          std::vector<std::string> f;
          pp.getarr("Bx_external_grid_function(x,y,z)", f);
          str_Bx_function.clear();
          for (auto const& s : f) {
              str_Bx_function += s;
          }
          f.clear();

          std::unique_ptr<ParserWrapper> Bx_parsewrap;
          Bx_parsewrap.reset(new ParserWrapper(makeParser(str_Bx_function)));
          ParserWrapper* Bx_wrap = (Bx_parsewrap.get());

       for ( MFIter mfi(*Bx, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
       {
          const Box& tbx = mfi.tilebox(Bx_nodal_flag);
          const Box& tby = mfi.tilebox(By_nodal_flag);
          const Box& tbz = mfi.tilebox(Bz_nodal_flag);
          // Bx 
          auto const& Bxfab = Bx->array(mfi);
          auto const& lo = lbound(tbx);
          auto const& hi = ubound(tbx);
          amrex::Print() << " lo " << lo << " hi " << hi << "\n";
          auto const& Bx_IndexType = (*Bx).ixType();
          IntVect Bx_type(AMREX_D_DECL(0,0,0));
          for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
              Bx_type[idim] = Bx_IndexType.nodeCentered(idim);
          }
          amrex::Print() << " ix type " << Bx_IndexType << "\n";
          amrex::Print() << " type : " << Bx_type << "\n";
          //ParserWrapper *Bx_parsewrap = new ParserWrapper(Bx_parser);
//          amrex::Print() << " parsewrap " << Bx_parsewrap->m_parser << "\n";
//          amrex::Print() <<  " parse wrap " << Bx_wrap->m_parser << "\n";
          Real xx = 1.0; Real yy = 0.0; Real zz = 0.0;
          Real g = Bx_wrap->m_parser(xx,yy,zz);
          amrex::Print() << "g " << g << "\n";
          Real g2 = Bx_wrap->getField(xx,yy,zz);
          amrex::Print() << "g2 " << g2 << "\n";
          amrex::ParallelFor (tbx, 
              [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                  int xdim = 0; int ydim = 1; int zdim = 2;
                  
                  //Real fac_x = (1.0 - Bx_type[0]) * dx_lev[0]*0.5;
                  //Real x = i*dx_lev[xdim] + real_box.lo(xdim) + fac_x;
                  //Real fac_y = (1.0 - Bx_type[1]) * dx_lev[1]*0.5;
                  //Real y = j*dx_lev[ydim] + real_box.lo(ydim) + fac_y;
                  Real x = 0.0, y = 0.0,  z = 0.0;
                  Bxfab(i,j,k) = Bx_wrap->getField(x,y,z);
          //        Bxfab(i,j,k) = 0.0;
                  
              }
          );
         
       //   for (int k = lo.z; k <= hi.z; ++k) {
       //   for (int j = lo.y; j <= hi.y; ++j) {
       //   for (int i = lo.x; i <= hi.x; ++i) {
       //       //std::vector<std::string> f;
       //       //pp.getarr("Bx_external_grid_function(x,y,z)", f);
       //       //str_Bx_ext_grid_function.clear();
       //       //for (auto const& s : f) {
       //       //    str_Bx_ext_grid_function += s;
       //       //}
       //       //f.clear();
       //       //WarpXParser const& Bx_parser = makeParser(str_Bx_ext_grid_function);
       //       // 
       //       //GpuParser m_Bx_parser = Bx_parser;

       //       //Real x,y,z;
       //       int xdim = 0; int ydim = 1; int zdim = 2;
       //       //
       //       Real fac_x = (1.0 - Bx_type[0]) * dx_lev[0]*0.5;
       //       Real x = i*dx_lev[xdim] + real_box.lo(xdim) + fac_x;
       //       Real fac_y = (1.0 - Bx_type[1]) * dx_lev[1]*0.5;
       //       Real y = j*dx_lev[ydim] + real_box.lo(ydim) + fac_y;
       //       Real z = 0.0;
       //       //Bxfab(i,j,k) = m_Bx_parser(x,y,z);
       //       //amrex::Print() << " facx " << fac_x << " facy " << fac_y << "\n";
       //       //amrex::Print() << " Bx at " << i << " " << j << " " << k << " is " << Bxfab(i,j,k) << " x is " << x << " z is "  << y <<  "\n";
       //       amrex::Print() << " x " << x << " " << y  << "\n";
       //       //Real compare = 200*std::cosh(x/2.0) + 10*y;
       //       Real compare = 200*std::cos(x/2.0) + 10*y;
       //       amrex::Print() << " Bx at " << i << " " << j << " " << k << " is " << Bxfab(i,j,k) << " compare " << compare << "\n";
       //   }}}          
       //   // By
       //   auto const& By_IndexType = (*By).ixType();
       //   auto const& Byfab = By->array(mfi);
       //   auto const& lo_y = lbound(tby);
       //   auto const& hi_y = ubound(tby);
       //   for (int k = lo_y.z; k <= hi_y.z; ++k) {
       //   for (int j = lo_y.y; j <= hi_y.y; ++j) {
       //   for (int i = lo_y.x; i <= hi_y.x; ++i) {
//     //         Byfab(i,j,k) = B_external_grid[1];
       //       Byfab(i,j,k) = 0;
       //       amrex::Print() << " By at " << i << " " << j << " " << k << " is " << Byfab(i,j,k) << "\n";
       //   }}}          
       //   
       //   // Bz
       //   auto const& Bz_IndexType = (*Bz).ixType();
       //   auto const& Bzfab = Bz->array(mfi);
       //   auto const& lo_z = lbound(tbz);
       //   auto const& hi_z = ubound(tbz);
       //   for (int k = lo_z.z; k <= hi_z.z; ++k) {
       //   for (int j = lo_z.y; j <= hi_z.y; ++j) {
       //   for (int i = lo_z.x; i <= hi_z.x; ++i) {
       //       //Bzfab(i,j,k) = B_external_grid[2];
       //       Bzfab(i,j,k) = 0;
       //       amrex::Print() << " Bz at " << i << " " << j << " " << k << " is " << Bzfab(i,j,k) << "\n";
       //   }}}          
       //   
       }
    }

    if (lev > 0) {
        for (int i = 0; i < 3; ++i) {
            current_cp[lev][i]->setVal(0.0);
            if (B_ext_grid_s == "constant") {
               Bfield_aux[lev][i]->setVal(B_external_grid[i]);
               Bfield_cp[lev][i]->setVal(B_external_grid[i]);
            }
            else {
               Bfield_aux[lev][i]->setVal(0.0);
               Bfield_cp[lev][i]->setVal(0.0);
            }
            if (E_ext_grid_s == "constant") {
               Efield_aux[lev][i]->setVal(E_external_grid[i]);
               Efield_cp[lev][i]->setVal(E_external_grid[i]);
            } else {
               Efield_aux[lev][i]->setVal(0.0);
               Efield_cp[lev][i]->setVal(0.0);
            }
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
