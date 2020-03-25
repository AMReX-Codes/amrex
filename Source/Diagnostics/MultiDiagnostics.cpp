#include "MultiDiagnostics.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

MultiDiagnostics::MultiDiagnostics ()
{
    ReadParameters();
    /** Resize alldiags and initialize each element to a pointer to a
     * diagnostics. Calls the corresponding diagnostics constructor.
     */
    alldiags.resize( ndiags );
    for (int i=0; i<ndiags; i++){
        if ( diags_types[i] == DiagTypes::Full ){
            alldiags[i].reset( new FullDiagnostics(i, diags_names[i]) );
        } else {
            amrex::Abort("Unknown diagnostic type");
        }
    }
}

void
MultiDiagnostics::InitData ()
{
    for( auto& diag : alldiags ){
        diag->InitData();
    }
}

void
MultiDiagnostics::ReadParameters ()
{
    ParmParse pp("diagnostics");
    pp.query("ndiags", ndiags);
    diags_types.resize( ndiags );
    if (ndiags > 0) pp.getarr("diags_names", diags_names);
    for (int i=0; i<ndiags; i++){
        ParmParse ppd(diags_names[i]);
        std::string diag_type_str;
        ppd.get("diag_type", diag_type_str);
        if (diag_type_str == "Full") diags_types[i] = DiagTypes::Full;
    }
}

void
MultiDiagnostics::FilterComputePackFlush (int step, bool force_flush)
{
    for (auto& diag : alldiags){
        if ( !diag->DoDump( step, force_flush ) ) return;
        diag->ComputeAndPack();
        diag->Flush();
        diag->FlushRaw();
    }
}
