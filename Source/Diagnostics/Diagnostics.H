#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include "FlushFormat.H"
#include "FlushFormatPlotfile.H"
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
/**
 * \brief base class for diagnostics.
 * Contains main routines to filter, compute and flush diagnostics.
 *
 * Each specific diagnostics derives from this class.
 */
class Diagnostics
{
public:
    Diagnostics (int i, std::string name);
    ~Diagnostics ();
    void ReadParameters ();
    /** Pack (stack) all fields in the cell-centered output MultiFab mf_avg */
    void FilterAndPack ();
    /** Flush mf_avg and particles to file. */
    void Flush ();
    /** Flush raw data */
    void FlushRaw ();
    /** Initialize pointers to main fields and allocate output multifab
     * mf_avg*/
    void InitData ();
    /** whether to flush at this time step
     * \param[in] step current time step
     * \param[in] force_flush if true, return true for any step
     * \return bool, whether to flush
     */
    bool DoDump (int step, bool force_flush=false);
private:
    /** Name of diagnostics: runtime parameter given in the input file. */
    std::string diag_name;
    /** Prefix for output directories */
    std::string file_prefix = "diags/diag";
    /** Index of diagnostics in MultiDiagnostics::alldiags */
    int diag_index;
    /** Period at which diagnostics is writte to file */
    int m_period = -1;
    /** Number of scalar components requested by the user */
    int ncomp {9};
    /** Names of  each component requested by the user */
    amrex::Vector< std::string > varnames;
    /** Pointer to all fields, per level
     * levels    < Ex Ey Ez Bx etc.                 > */
    amrex::Vector< amrex::Vector <amrex::MultiFab const *> > allfields;
    /** output multifab, where all fields are cell-centered and stacked */
    amrex::Vector< amrex::MultiFab > mf_avg;
    /** number of levels to output */
    int nlev;
    /** This class is responsible for flushing the data to file */
    FlushFormat* m_flush_format;
};

#endif // WARPX_DIAGNOSTICS_H_
