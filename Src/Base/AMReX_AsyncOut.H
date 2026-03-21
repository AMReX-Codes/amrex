#ifndef AMREX_ASYNCOUT_H_
#define AMREX_ASYNCOUT_H_
#include <AMReX_Config.H>

#include <functional>

/**
 * \file AMReX_AsyncOut.H
 * \brief Asynchronous output interface backed by a background writer thread.
 */

namespace amrex::AsyncOut {

/**
 * \struct amrex::AsyncOut::WriteInfo
 * \brief Slot assignment metadata for ranks sharing an output file.
 * \ingroup amrex_io
 */
struct WriteInfo {
    int ifile;
    int ispot;
    int nspots;
};

/**
 * \brief Initialize the asynchronous output infrastructure.
 * \ingroup amrex_io
 */
void Initialize ();
/**
 * \brief Tear down resources owned by the asynchronous output subsystem.
 * \ingroup amrex_io
 */
void Finalize ();

/**
 * \brief Query whether asynchronous output is enabled at runtime.
 * \ingroup amrex_io
 */
bool UseAsyncOut ();

/**
 * \brief Compute this rank's write slot information for multi-file output.
 * \param rank MPI rank requesting metadata.
 * \return WriteInfo describing file index, slot, and slots per file.
 * \ingroup amrex_io
 */
WriteInfo GetWriteInfo (int rank);

/**
 * \brief Submit a job functor to the asynchronous writer thread.
 * \param a_f Functor that performs the output operations.
 * \ingroup amrex_io
 */
void Submit (std::function<void()>&& a_f);
/**
 * \overload
 * \ingroup amrex_io
 */
void Submit (std::function<void()> const& a_f);

/**
 * \brief Block until all submitted asynchronous jobs have completed.
 * \ingroup amrex_io
 */
void Finish ();

/**
 * \brief Within a submitted job, wait for this rank's turn to access the file.
 *
 * This is not for waiting for the job to finish; use Finish() for that.
 *
 * \ingroup amrex_io
 */
void Wait ();
/**
 * \brief Within a submitted job, notify the next rank sharing the file that it may proceed.
 * \ingroup amrex_io
 */
void Notify ();

}

#endif
