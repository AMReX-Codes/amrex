#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include "FlushFormats/FlushFormat.H"
#include "FlushFormats/FlushFormatPlotfile.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
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
    /** Pack (stack) all fields in the cell-centered output MultiFab mf_avg.
     * Fields are computed (e.g., cell-centered) on-the-fly using a functor.*/
    void ComputeAndPack ();
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
    void ReadParameters ();
    /** Append varnames with names for all modes of a field
     * \param[in] field field name (includes component, e.g., Er)
     * \param[in] ncomp number of components (modes, real and imag)
     */
    void AddRZModesToOutputNames (const std::string& field, int ncomp);
    /** Append all_field_functors[lev] with all modes of all components of E B and j
     * and add name to varnames.
     * \param[in] lev level on which source multifabs are defined.
     */
    void AddRZModesToDiags (int lev);
    /** Name of diagnostics: runtime parameter given in the input file. */
    std::string diag_name;
    /** Prefix for output directories */
    std::string file_prefix;
    /** Index of diagnostics in MultiDiagnostics::alldiags */
    int diag_index;
    /** Period at which diagnostics is writte to file */
    int m_period = -1;
    /** Names of  each component requested by the user.
     * in cylindrical geometry, this list is appended with
     * automatically-constructed names for all modes of all fields */
    amrex::Vector< std::string > varnames;
    /** Vector of (pointers to) functors to compute output fields, per level,
     * per component. This allows for simple operations (averaging to
     * cell-center for standard EB fields) as well as more involved operations
     * (back-transformed diagnostics, filtering, reconstructing cartesian
     * fields in cylindrical). */
    amrex::Vector< amrex::Vector <ComputeDiagFunctor const *> > all_field_functors;
    /** output multifab, where all fields are cell-centered and stacked */
    amrex::Vector< amrex::MultiFab > mf_avg;
    int nlev; /**< number of levels to output */
    /** This class is responsible for flushing the data to file */
    FlushFormat* m_flush_format;
    /** Whether to plot raw (i.e., NOT cell-centered) fields */
    bool m_plot_raw_fields = false;
    /** Whether to plot guard cells of raw fields */
    bool m_plot_raw_fields_guards = false;
    /** Whether to plot charge density rho */
    bool m_plot_rho = false;
    /** Whether to plot F (charge conservation error) */
    bool m_plot_F = false;
};

#endif // WARPX_DIAGNOSTICS_H_
