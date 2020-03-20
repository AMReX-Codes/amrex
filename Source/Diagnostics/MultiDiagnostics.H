#ifndef WARPX_MULTIDIAGNOSTICS_H_
#define WARPX_MULTIDIAGNOSTICS_H_

#include "FullDiagnostics.H"

/** All types of diagnostics. */
enum struct DiagTypes {Full};

/**
 * \brief This class contains a vector of all diagnostics in the simulation.
 */
class MultiDiagnostics
{
public:
    MultiDiagnostics ();
    /** \brief Read Input parameters. Called in constructor. */
    void ReadParameters ();
    /** \brief Loop over diags in alldiags and call their InitDiags */
    void InitData ();
    /** \brief Called at each iteration. Compute diags and flush. */
    void FilterComputePackFlush (int step, bool force_flush=false);
private:
    /** Vector of pointers to all diagnostics */
    amrex::Vector<std::unique_ptr<Diagnostics> > alldiags;
    int ndiags = 0; /**< number of different diagnostics */
    std::vector<std::string> diags_names;
    /**Type of each diagnostics*/
    std::vector<DiagTypes> diags_types;
};

#endif // WARPX_MULTIDIAGNOSTICS_H_
