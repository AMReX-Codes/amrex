#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include "FlushFormats/FlushFormat.H"
#include "FlushFormats/FlushFormatPlotfile.H"
#ifdef WARPX_USE_OPENPMD
#   include "FlushFormats/FlushFormatOpenPMD.H"
#endif
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
    /** Pack (stack) all fields in the cell-centered output MultiFab m_mf_output.
     * Fields are computed (e.g., cell-centered) on-the-fly using a functor.*/
    void ComputeAndPack ();
    /** Flush m_mf_output and particles to file. */
    void Flush ();
    /** Flush raw data */
    void FlushRaw ();
    /** Initialize pointers to main fields and allocate output multifab m_mf_output. */
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
    /** Append m_all_field_functors[lev] with all modes of all components of E B and j
     * and add name to varnames.
     * \param[in] lev level on which source multifabs are defined.
     */
    void AddRZModesToDiags (int lev);
    /** Name of diagnostics: runtime parameter given in the input file. */
    std::string m_diag_name;
    /** Prefix for output directories */
    std::string m_file_prefix;
    /** Index of diagnostics in MultiDiagnostics::alldiags */
    int m_diag_index;
    /** Period at which diagnostics is written to file */
    int m_period = -1;
    /** Names of  each component requested by the user.
     * in cylindrical geometry, this list is appended with
     * automatically-constructed names for all modes of all fields */
    amrex::Vector< std::string > m_varnames;
    /** Vector of (pointers to) functors to compute output fields, per level,
     * per component. This allows for simple operations (averaging to
     * cell-center for standard EB fields) as well as more involved operations
     * (back-transformed diagnostics, filtering, reconstructing cartesian
     * fields in cylindrical). */
    amrex::Vector< amrex::Vector <ComputeDiagFunctor const *> > m_all_field_functors;
    /** output multifab, where all fields are cell-centered and stacked */
    amrex::Vector< amrex::MultiFab > m_mf_output;
    int nlev; /**< number of levels to output */
    /** This class is responsible for flushing the data to file */
    FlushFormat* m_flush_format;
    /** Whether to plot raw (i.e., NOT cell-centered) fields */
    bool m_plot_raw_fields = false;
    /** Whether to plot guard cells of raw fields */
    bool m_plot_raw_fields_guards = false;
    /** Whether to plot charge density rho in raw fields */
    bool m_plot_raw_rho = false;
    /** Whether to plot F (charge conservation error) in raw fields */
    bool m_plot_raw_F = false;

    /** Coarsening ratio such that, fields are averaged to the coarsened grid.
      * The ratio should render the grid to be coarsenable (as defined by AMReX). */
    amrex::IntVect m_crse_ratio = amrex::IntVect(1);
    /** Lower corner of the diagnostics output, in physical coordinates */
    amrex::Vector< amrex::Real> m_lo;
    /** Higher corner of the diagnostics output, in physical coordinates */
    amrex::Vector< amrex::Real> m_hi;
    /** Define the cell-centered multifab m_mf_output depending on user-defined
      * lo and hi and coarsening ratio. This MultiFab may have a different BoxArray and
      * DistributionMap than field MultiFabs in the simulation.
      *
      * \param[in] lev level on which source multifabs are defined
      */
    void DefineDiagMultiFab ( int lev );

    /** format for output files, "plotfile" or "openpmd" */
    std::string m_format = "plotfile";
};

#endif // WARPX_DIAGNOSTICS_H_
