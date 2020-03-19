
#include "Diagnostics.H"
#include "WarpX.H"
#include "Average.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_index(i), diag_name(name)
{}

Diagnostics::~Diagnostics()
{
    delete m_flush_format;
    /*
    for ( int lev=0; lev<nlev; lev++ ){
        for ( int i=0; i<ncomp; i++ ){
            allfields[lev][i].reset();
        }
    }
    */
}

void
Diagnostics::InitData ()
{
    Print()<<"Diagnostics::InitData\n";

    auto & warpx = WarpX::GetInstance();
    nlev = warpx.finestLevel() + 1;
    // Initialize vector of pointer to the fields requested by the user.
    allfields.resize( nlev );
    mf_avg.resize( nlev );
    for ( int lev=0; lev<nlev; lev++ ){
        allfields[lev].resize( ncomp );
        for ( int dim=0; dim<3; dim++ ){
            allfields[lev][dim  ] = warpx.get_pointer_Efield_aux(lev, dim);
            allfields[lev][dim+3] = warpx.get_pointer_Bfield_aux(lev, dim);
            allfields[lev][dim+6] = warpx.get_pointer_current_fp(lev, dim);
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
