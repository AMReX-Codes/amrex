#ifndef AMREX_PARSTREAM_H_
#define AMREX_PARSTREAM_H_

#include <iostream>
#include <string>
namespace amrex
{

/**
 * \brief Use this in place of std::cout for program output.
 *
 * Can be used to replace std::cout.  In serial this just
 * returns std::cout.  In parallel, this creates a separate file for each
 * proc called \<basename\>.n where n is the procID and \<basename\> defaults to
 * "pout" but can be set by calling setPoutBaseName().  Output is then
 * directed to these files.  This keeps the output from different processors
 * from getting all jumbled up.
 If    you want fewer files, you can use ParmParse parameter amrex.pout_int=nproc
 and it will only output every nproc processors pout.n files
 (where n%nproc == 0).
*/
  std::ostream& pout();

/**
 * \brief Changes the base part of the filename for pout() files.
 *
 * When in parallel, changes the base name of the pout() files.  If pout()
 * has already been called, it closes the current output file and opens a
 * new one (unless the name is the same, in which case it does nothing).
 * In serial, ignores the argument and does nothing.
 */
  void setPoutBaseName(const std::string & a_Name);

/**
 * \brief Accesses the filename for the local pout() file.
 *
 * Returns the name used for the local pout() file.  In parallel this is
 * "\<pout_basename\>.\<procID\>", where \<pout_basename\> defaults to "pout" and
 * can be modified by calling setPoutBaseName(), and \<procID\> is the local
 * proc number.  In serial, this always returns the string "cout".  It is an
 * error (exit code 111) to call this in parallel before MPI_Initialize().
 */
  const std::string & poutFileName();

}
#endif // _PARSTREAM_H_
