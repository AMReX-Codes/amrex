#include "Diagnostics.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/PartPerCellFunctor.H"
#include "ComputeDiagFunctors/PartPerGridFunctor.H"
#include "ComputeDiagFunctors/DivBFunctor.H"
#include "ComputeDiagFunctors/DivEFunctor.H"
#include "WarpX.H"
#include "Utils/Average.H"
#include "Utils/WarpXUtil.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : m_diag_name(name), m_diag_index(i)
{
    ReadParameters();
}

Diagnostics::~Diagnostics ()
{
    delete m_flush_format;
}

void
Diagnostics::ReadParameters ()
{
    auto & warpx = WarpX::GetInstance();
    // Read list of fields requested by the user.
    ParmParse pp(m_diag_name);
    m_file_prefix = "diags/" + m_diag_name;
    pp.query("file_prefix", m_file_prefix);
    pp.query("period", m_period);
    pp.query("format", m_format);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile" || m_format == "openpmd",
        "<diag>.format must be plotfile or openpmd");
    pp.query("plot_raw_fields", m_plot_raw_fields);
    pp.query("plot_raw_fields_guards", m_plot_raw_fields_guards);
    if (!pp.queryarr("fields_to_plot", m_varnames)){
        m_varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    }
    // set plot_rho to true of the users requests it, so that
    // rho is computed at each iteration.
    if (WarpXUtilStr::is_in(m_varnames, "rho")) warpx.setplot_rho(true);
    // Sanity check if user requests to plot F
    if (WarpXUtilStr::is_in(m_varnames, "F")){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            warpx.do_dive_cleaning,
            "plot F only works if warpx.do_dive_cleaning = 1");
    }
    // If user requests to plot proc_number for a serial run,
    // delete proc_number from fields_to_plot
    if (ParallelDescriptor::NProcs() == 1){
        m_varnames.erase(
            std::remove(m_varnames.begin(), m_varnames.end(), "proc_number"),
            m_varnames.end());
    }


    // Read user-defined physical extents for the output and store in m_lo and m_hi.
    m_lo.resize(AMREX_SPACEDIM);
    m_hi.resize(AMREX_SPACEDIM);

    if (!pp.queryarr("diag_lo", m_lo) ) {
       for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
            m_lo[idim] = warpx.Geom(0).ProbLo(idim);
       }
    }
    if (! pp.queryarr("diag_hi", m_hi) ) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
            m_hi[idim] = warpx.Geom(0).ProbHi(idim);
       }
    }

    // Initialize cr_ratio with default value of 1 for each dimension.
    Vector<int> cr_ratio(AMREX_SPACEDIM, 1);
    // Read user-defined coarsening ratio for the output MultiFab.
    if (pp.queryarr("coarsening_ratio", cr_ratio) ) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
           m_crse_ratio[idim] = cr_ratio[idim];
       }
    }

}

void
Diagnostics::InitData ()
{
    Print()<<"Diagnostics::InitData\n";
    auto & warpx = WarpX::GetInstance();
    nlev = warpx.finestLevel() + 1;
    // Initialize vector of pointers to the fields requested by the user.
    m_all_field_functors.resize( nlev );
    m_mf_output.resize( nlev );

    for ( int lev=0; lev<nlev; lev++ ){
        m_all_field_functors[lev].resize( m_varnames.size() );
        // Fill vector of functors for all components except individual cylindrical modes.
        for (int comp=0, n=m_all_field_functors[lev].size(); comp<n; comp++){
            if        ( m_varnames[comp] == "Ex" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 0), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "Ey" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 1), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "Ez" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 2), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "Bx" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 0), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "By" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 1), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "Bz" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 2), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "jx" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 0), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "jy" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 1), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "jz" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 2), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "rho" ){
                // rho_new is stored in component 1 of rho_fp when using PSATD
#ifdef WARPX_USE_PSATD
                MultiFab* rho_new = new MultiFab(*warpx.get_pointer_rho_fp(lev), amrex::make_alias, 1, 1);
                m_all_field_functors[lev][comp] = new CellCenterFunctor(rho_new, lev, m_crse_ratio);
#else
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_rho_fp(lev), lev, m_crse_ratio);
#endif
            } else if ( m_varnames[comp] == "F" ){
                m_all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_F_fp(lev), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "part_per_cell" ){
                m_all_field_functors[lev][comp] = new PartPerCellFunctor(nullptr, lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "part_per_grid" ){
                m_all_field_functors[lev][comp] = new PartPerGridFunctor(nullptr, lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "divB" ){
                m_all_field_functors[lev][comp] = new DivBFunctor(warpx.get_array_Bfield_aux(lev), lev, m_crse_ratio);
            } else if ( m_varnames[comp] == "divE" ){
                m_all_field_functors[lev][comp] = new DivEFunctor(warpx.get_array_Efield_aux(lev), lev, m_crse_ratio);
            }
        }

        AddRZModesToDiags( lev );

        // At this point, m_varnames.size() >= m_all_field_functors[0].size()

        // Initialize member variable m_mf_output depending on m_crse_ratio, m_lo and m_hi
        DefineDiagMultiFab( lev );

    }
    // Construct Flush class.
    if        (m_format == "plotfile"){
        m_flush_format = new FlushFormatPlotfile;
    } else if (m_format == "openpmd"){
#ifdef WARPX_USE_OPENPMD
        m_flush_format = new FlushFormatOpenPMD(m_diag_name);
#else
        amrex::Abort("To use openpmd output format, need to compile with USE_OPENPMD=TRUE");
#endif
    } else {
        amrex::Abort("unknown output format");
    }
}

void
Diagnostics::ComputeAndPack ()
{
    // First, make sure all guard cells are properly filled
    // Probably overkill/unnecessary, but safe and shouldn't happen often !!
    auto & warpx = WarpX::GetInstance();
    warpx.FillBoundaryE(warpx.getngE(), warpx.getngExtra());
    warpx.FillBoundaryB(warpx.getngE(), warpx.getngExtra());
#ifndef WARPX_USE_PSATD
    warpx.FillBoundaryAux(warpx.getngUpdateAux());
#endif
    warpx.UpdateAuxilaryData();

    warpx.FieldGather();

    // cell-center fields and store result in m_mf_output.
    int icomp_dst = 0;
    for(int lev=0; lev<nlev; lev++){
        for (int icomp=0, n=m_all_field_functors[0].size(); icomp<n; icomp++){
            // Call all functors in m_all_field_functors[lev]. Each of them computes
            // a diagnostics and writes in one or more components of the output
            // multifab m_mf_output[lev].
            m_all_field_functors[lev][icomp]->operator()(m_mf_output[lev], icomp_dst);
            // update the index of the next component to fill
            icomp_dst += m_all_field_functors[lev][icomp]->nComp();
        }
    }
    // Check that the proper number of components of m_mf_output were updated.
    AMREX_ALWAYS_ASSERT( icomp_dst == m_varnames.size() );
}

void
Diagnostics::Flush ()
{
    auto & warpx = WarpX::GetInstance();
    m_flush_format->WriteToFile(
        m_varnames, GetVecOfConstPtrs(m_mf_output), warpx.Geom(), warpx.getistep(),
        warpx.gett_new(0), warpx.GetPartContainer(), nlev, m_file_prefix,
        m_plot_raw_fields, m_plot_raw_fields_guards, m_plot_raw_rho, m_plot_raw_F);
}

void
Diagnostics::FlushRaw () {}

bool
Diagnostics::DoDump (int step, bool force_flush)
{
    if (force_flush) return true;
    if ( m_period>0 && (step+1)%m_period==0 ) return true;
    return false;
}

void
Diagnostics::AddRZModesToDiags (int lev)
{
#ifdef WARPX_DIM_RZ
    auto & warpx = WarpX::GetInstance();
    int ncomp_multimodefab = warpx.get_pointer_Efield_aux(0, 0)->nComp();
    // Make sure all multifabs have the same number of components
    for (int dim=0; dim<3; dim++){
        AMREX_ALWAYS_ASSERT(
            warpx.get_pointer_Efield_aux(lev, dim)->nComp() == ncomp_multimodefab );
        AMREX_ALWAYS_ASSERT(
            warpx.get_pointer_Bfield_aux(lev, dim)->nComp() == ncomp_multimodefab );
        AMREX_ALWAYS_ASSERT(
            warpx.get_pointer_current_fp(lev, dim)->nComp() == ncomp_multimodefab );
    }

    // First index of m_all_field_functors[lev] where RZ modes are stored
    int icomp = m_all_field_functors[0].size();
    const std::array<std::string, 3> coord {"r", "theta", "z"};

    // Er, Etheta, Ez, Br, Btheta, Bz, jr, jtheta, jz
    // Each of them being a multi-component multifab
    m_all_field_functors[lev].resize( m_all_field_functors[0].size() + 9 );
    // E
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        m_all_field_functors[lev][icomp] = new
            CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        AddRZModesToOutputNames(std::string("E") + coord[dim],
                                warpx.get_pointer_Efield_aux(0, 0)->nComp());
        icomp += 1;
    }
    // B
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        m_all_field_functors[lev][icomp] = new
            CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        AddRZModesToOutputNames(std::string("B") + coord[dim],
                                warpx.get_pointer_Bfield_aux(0, 0)->nComp());
        icomp += 1;
    }
    // j
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        m_all_field_functors[lev][icomp] = new
            CellCenterFunctor(warpx.get_pointer_current_fp(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        icomp += 1;
        AddRZModesToOutputNames(std::string("J") + coord[dim],
                                warpx.get_pointer_current_fp(0, 0)->nComp());
    }
    // Sum the number of components in input vector m_all_field_functors
    // and check that it corresponds to the number of components in m_varnames
    // and m_mf_output
    int ncomp_from_src = 0;
    for (int i=0; i<m_all_field_functors[0].size(); i++){
        ncomp_from_src += m_all_field_functors[lev][i]->nComp();
    }
    AMREX_ALWAYS_ASSERT( ncomp_from_src == m_varnames.size() );
#endif
}

void
Diagnostics::AddRZModesToOutputNames (const std::string& field, int ncomp){
#ifdef WARPX_DIM_RZ
    // In cylindrical geometry, real and imag part of each mode are also
    // dumped to file separately, so they need to be added to m_varnames
    m_varnames.push_back( field + "_0_real" );
    for (int ic=1 ; ic < ncomp ; ic += 2) {
        m_varnames.push_back( field + "_" + std::to_string(ic) + "_real" );
        m_varnames.push_back( field + "_" + std::to_string(ic) + "_imag" );
    }
#endif
}

void
Diagnostics::DefineDiagMultiFab ( int lev ) {
    auto & warpx = WarpX::GetInstance();
    amrex::RealBox diag_dom;
    bool use_warpxba = true;
    const IntVect blockingFactor = warpx.blockingFactor( lev );

    // Default BoxArray and DistributionMap for initializing the output MultiFab, m_mf_output.
    BoxArray ba = warpx.boxArray(lev);
    DistributionMapping dmap = warpx.DistributionMap(lev);

    // Check if warpx BoxArray is coarsenable.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE (
        ba.coarsenable(m_crse_ratio),
        "Invalid coarsening ratio for warpx boxArray. Must be an integer divisor of the blocking factor."
    );

    // Find if user-defined physical dimensions are different from the simulation domain.
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
         // To ensure that the diagnostic lo and hi are within the domain defined at level, lev.
        diag_dom.setLo(idim, max(m_lo[idim],warpx.Geom(lev).ProbLo(idim)) );
        diag_dom.setHi(idim, min(m_hi[idim],warpx.Geom(lev).ProbHi(idim)) );
        if ( fabs(warpx.Geom(lev).ProbLo(idim) - diag_dom.lo(idim))
                               >  warpx.Geom(lev).CellSize(idim) )
             use_warpxba = false;
        if ( fabs(warpx.Geom(lev).ProbHi(idim) - diag_dom.hi(idim))
                               > warpx.Geom(lev).CellSize(idim) )
             use_warpxba = false;

        // User-defined value for coarsening should be an integer divisor of
        // blocking factor at level, lev.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( blockingFactor[idim] % m_crse_ratio[idim]==0,
                       " coarsening ratio must be integer divisor of blocking factor");
    }


    if (use_warpxba == false) {
        // Following are the steps to compute the lo and hi index corresponding to user-defined
        // m_lo and m_hi using the same resolution as the simulation at level, lev.
        IntVect lo(0);
        IntVect hi(1);
        for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
            // lo index with same cell-size as simulation at level, lev.
            lo[idim] = max( static_cast<int>( floor (
                          ( diag_dom.lo(idim) - warpx.Geom(lev).ProbLo(idim)) /
                            warpx.Geom(lev).CellSize(idim)) ), 0 );
            // hi index with same cell-size as simulation at level, lev.
            hi[idim] = max( static_cast<int> ( ceil (
                          ( diag_dom.hi(idim) - warpx.Geom(lev).ProbLo(idim)) /
                            warpx.Geom(lev).CellSize(idim) ) ), 0) - 1 ;
            // if hi<=lo, then hi = lo + 1, to ensure one cell in that dimension
            if ( hi[idim] <= lo[idim] ) {
                 hi[idim]  = lo[idim] + 1;
                 AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    m_crse_ratio[idim]==1, "coarsening ratio in reduced dimension must be 1."
                 );
            }
        }

        // Box for the output MultiFab corresponding to the user-defined physical co-ordinates at lev.
        Box diag_box( lo, hi );
        // Define box array
        BoxArray diag_ba;
        diag_ba.define(diag_box);
        ba = diag_ba.maxSize( warpx.maxGridSize( lev ) );
        // At this point in the code, the BoxArray, ba, is defined with the same index space and
        // resolution as the simulation, at level, lev.
        // Coarsen and refine so that the new BoxArray is coarsenable.
        ba.coarsen(m_crse_ratio).refine(m_crse_ratio);

        // Update the physical co-ordinates m_lo and m_hi using the final index values
        // from the coarsenable, cell-centered BoxArray, ba.
        for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_lo[idim] = warpx.Geom(lev).ProbLo(idim) + warpx.Geom(lev).CellSize(idim)/2.0_rt +
                ba.getCellCenteredBox(0).smallEnd(idim) * warpx.Geom(lev).CellSize(idim);
            m_hi[idim] = warpx.Geom(lev).ProbLo(idim) + warpx.Geom(lev).CellSize(idim)/2.0_rt +
                ba.getCellCenteredBox( ba.size()-1 ).bigEnd(idim) * warpx.Geom(lev).CellSize(idim);
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_crse_ratio.min() > 0, "Coarsening ratio must be non-zero.");
    // The BoxArray is coarsened based on the user-defined coarsening ratio.
    ba.coarsen(m_crse_ratio);
    // Generate a new distribution map if the physical m_lo and m_hi for the output
    // is different from the lo and hi physical co-ordinates of the simulation domain.
    if (use_warpxba == false) dmap = DistributionMapping{ba};
    // Allocate output MultiFab for diagnostics. The data will be stored at cell-centers.
    m_mf_output[lev] = MultiFab(ba, dmap, m_varnames.size(), 0);
}
