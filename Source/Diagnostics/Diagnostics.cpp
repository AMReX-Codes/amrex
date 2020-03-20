
#include "Diagnostics.H"
#include "WarpX.H"
#include "Average.H"
#include "Utils/WarpXUtil.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_index(i), diag_name(name)
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
    ParmParse pp(diag_name);
    if (!pp.queryarr("fields_to_plot", varnames)){
        varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    }
    ncomp = varnames.size();
    // set plot_rho to true of the users requests it, so that
    // rho is computed at each iteration.
    if (WarpXUtilStr::is_in(varnames, "rho")) warpx.setplot_rho(true);
    // Sanity check if user requests to plot F
    if (WarpXUtilStr::is_in(varnames, "F")){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            warpx.do_dive_cleaning,
            "plot F only works if warpx.do_dive_cleaning = 1");
    }
    // If user requests to plot proc_number for a serial run,
    // delete proc_number from fields_to_plot
    if (ParallelDescriptor::NProcs() == 1){
        varnames.erase(
            std::remove(varnames.begin(), varnames.end(), "proc_number"),
            varnames.end());
    }
}

void
Diagnostics::InitData ()
{
    Print()<<"Diagnostics::InitData\n";
    auto & warpx = WarpX::GetInstance();
    nlev = warpx.finestLevel() + 1;
    // Initialize vector of pointers to the fields requested by the user.
    allfields.resize( nlev );
    mf_avg.resize( nlev );
    for ( int lev=0; lev<nlev; lev++ ){
        allfields[lev].resize( ncomp );
        for (int comp=0; comp<ncomp; comp++){
            if        ( varnames[comp] == "Ex" ){
                allfields[lev][comp] = warpx.get_pointer_Efield_aux(lev, 0);
            } else if ( varnames[comp] == "Ey" ){
                allfields[lev][comp] = warpx.get_pointer_Efield_aux(lev, 1);
            } else if ( varnames[comp] == "Ez" ){
                allfields[lev][comp] = warpx.get_pointer_Efield_aux(lev, 2);
            } else if ( varnames[comp] == "Bx" ){
                allfields[lev][comp] = warpx.get_pointer_Bfield_aux(lev, 0);
            } else if ( varnames[comp] == "By" ){
                allfields[lev][comp] = warpx.get_pointer_Bfield_aux(lev, 1);
            } else if ( varnames[comp] == "Bz" ){
                allfields[lev][comp] = warpx.get_pointer_Bfield_aux(lev, 2);
            } else if ( varnames[comp] == "jx" ){
                allfields[lev][comp] = warpx.get_pointer_current_fp(lev, 0);
            } else if ( varnames[comp] == "jy" ){
                allfields[lev][comp] = warpx.get_pointer_current_fp(lev, 1);
            } else if ( varnames[comp] == "jz" ){
                allfields[lev][comp] = warpx.get_pointer_current_fp(lev, 2);
            } else if ( varnames[comp] == "rho" ){
                allfields[lev][comp] = warpx.get_pointer_rho_fp(lev);
            } else if ( varnames[comp] == "part_per_cell" ){
                amrex::Abort("plotting part_per_cell is not supported");
            } else {
                amrex::Abort("Unknown field in fields_to_plot");
            }
        }
        // Allocate output multifab
        // Note: default MultiFab constructor is cell-centered
        mf_avg[lev] = MultiFab(warpx.boxArray(lev),
                               warpx.DistributionMap(lev),
                               ncomp, 0);
    }
    // Construct Flush class. So far, only Plotfile is implemented.
    m_flush_format = new FlushFormatPlotfile;
}

void
Diagnostics::FilterAndPack ()
{
    // cell-center fields and store result in mf_avg.
    for(int lev=0; lev<nlev; lev++){
        for (int icomp=0; icomp<ncomp; icomp++){
            Average::ToCellCenter ( mf_avg[lev],
                                    *allfields[lev][icomp],
                                    icomp, 0 );
        }
    }
}

void
Diagnostics::Flush ()
{
    auto & warpx = WarpX::GetInstance();
    m_flush_format->WriteToFile(varnames, GetVecOfConstPtrs(mf_avg),
                                warpx.Geom(), warpx.getistep(), 0.,
                                warpx.GetPartContainer(), nlev);
}

void
Diagnostics::FlushRaw () {}
