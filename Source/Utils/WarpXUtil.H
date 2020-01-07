#ifndef WARPX_UTILS_H_
#define WARPX_UTILS_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <WarpXParser.H>

#include <string>

void ReadBoostedFrameParameters(amrex::Real& gamma_boost, amrex::Real& beta_boost,
                                amrex::Vector<int>& boost_direction);

void ConvertLabParamsToBoost();

void NullifyMF(amrex::MultiFab& mf, int lev, amrex::Real zmin,
               amrex::Real zmax);

/**
 * \brief Parse a string (typically a mathematical expression) from the
 * input file and store it into a variable.
 *
 * \param parmparse_string String used to initialize ParmParse
 * \param query_string ParmParse.query will look for this string
 * \param stored_string variable in which the string to parse is stored
 */
void Store_parserString(std::string query_string,
                        std::string& stored_string);

namespace WarpXUtilIO{
    /**
     * A helper function to write binary data on disk.
     * @param[in] filename where to write
     * @param[in] data Vector containing binary data to write on disk
     * return true if it succeeds, false otherwise
     */
    bool WriteBinaryDataOnFile(std::string filename, const amrex::Vector<char>& data);
}

/**
* \brief Initialize a WarpXParser object from a string containing a math expression
*
* \param parse_function String to read to initialize the parser.
*/
WarpXParser makeParser (std::string const& parse_function);

#endif //WARPX_UTILS_H_
