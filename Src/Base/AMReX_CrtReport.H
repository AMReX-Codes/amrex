#ifndef AMREX_CRT_REPORT_H_
#define AMREX_CRT_REPORT_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_CrtReport.H
 * \brief Handling of reports of the debug Microsoft C runtime library.
 *
 * See the runtime parameter amrex.handle_crt_reports.
 */

namespace amrex::detail
{
    /** \brief Write CRT reports that would open a dialog box to stderr, then break.
     *
     * Errors and assertions of the debug Microsoft C runtime library that
     * would open a modal dialog box are written with the rank to stderr,
     * followed by a break into the debugger. Without a debugger, by default,
     * the process exits with STATUS_BREAKPOINT and Windows Error Reporting can
     * write a crash dump. Warnings and report modes without a dialog box are
     * kept. Does nothing in other builds, or if already set up. Call after
     * ParallelDescriptor::StartParallel.
     */
    void CrtReportInitialize ();

    /** \brief Skip the CRT report of abort(), or restore it.
     *
     * While skipped, abort() raises SIGABRT without a CRT report, for the
     * SIGABRT handler of AMReX. Does nothing in other builds.
     */
    void CrtReportSkipAbortReport (bool skip);

    /** \brief Undo CrtReportInitialize and CrtReportSkipAbortReport. */
    void CrtReportFinalize ();
}

#endif
