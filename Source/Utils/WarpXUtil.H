/* Copyright 2019-2020 Andrew Myers, Luca Fedeli, Maxence Thevenet
 * Revathi Jambunathan, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_UTILS_H_
#define WARPX_UTILS_H_

#include "Parser/WarpXParser.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

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
 * \param ParmParse pp used to read the query_string pp.<function>=string
 * \param parmparse_string String used to initialize ParmParse
 * \param query_string ParmParse.query will look for this string
 * \param stored_string variable in which the string to parse is stored
 */
void Store_parserString(amrex::ParmParse &pp, std::string query_string,
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

namespace WarpXUtilAlgo{

/** \brief Returns a pointer to the first element in the range [first, last) that is greater than val
 *
 * A re-implementation of the upper_bound algorithm suitable for GPU kernels.
 *
 * @param first: pointer to left limit of the range to consider
 * @param last: pointer to right limit of the range to consider
 * @param val: value to compare the elements of [first, last) to
 */
template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
const T* upper_bound(const T* first, const T* last, const T& val)
{
    const T* it;
    size_t count, step;
    count = last-first;
    while(count>0){
        it = first;
        step = count/2;
        it += step;
        if (!(val<*it)){
            first = ++it;
                count -= step + 1;
        }
        else{
            count = step;
        }
    }
        return first;
}

/** \brief Performs a linear interpolation
 *
 * Performs a linear interpolation at x given the 2 points
 * (x0, f0) and (x1, f1)
 */
template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T linear_interp(T x0, T x1, T f0, T f1, T x)
{
    return ((x1-x)*f0 + (x-x0)*f1)/(x1-x0);
}

/** \brief Performs a bilinear interpolation
 *
 * Performs a bilinear interpolation at (x,y) given the 4 points
 * (x0, y0, f00), (x0, y1, f01), (x1, y0, f10), (x1, y1, f11).
 */
template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T bilinear_interp(T x0, T x1, T y0, T y1, T f00, T f01, T f10, T f11, T x, T y)
{
    const T fx0 = linear_interp(x0, x1, f00, f10, x);
    const T fx1 = linear_interp(x0, x1, f01, f11, x);
    return linear_interp(y0, y1, fx0, fx1, y);
}

/** \brief Performs a trilinear interpolation
 *
 * Performs a trilinear interpolation at (x,y,z) given the 8 points
 * (x0, y0, z0, f000), (x0, y0, z1, f001), (x0, y1, z0, f010), (x0, y1, z1, f011),
 * (x1, y0, z0, f100), (x1, y0, z1, f101), (x1, y1, z0, f110), (x1, y1, z1, f111)
 */
template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T trilinear_interp(T x0, T x1,T y0, T y1, T z0, T z1,
    T f000, T f001, T f010, T f011, T f100, T f101, T f110, T f111,
    T x, T y, T z)
{
    const T fxy0 = bilinear_interp(
        x0, x1, y0, y1, f000, f010, f100, f110, x, y);
    const T fxy1 = bilinear_interp(
        x0, x1, y0, y1, f001, f011, f101, f111, x, y);
    return linear_interp(z0, z1, fxy0, fxy1, z);
}

}

/**
* \brief Initialize a WarpXParser object from a string containing a math expression
*
* \param parse_function String to read to initialize the parser.
*/
WarpXParser makeParser (std::string const& parse_function, std::vector<std::string> const& varnames);

namespace WarpXUtilMsg{

/** \brief If is_expression_true is false, this function prints msg and calls amrex::abort()
 *
 * @param[in] is_expression_true
 * @param[in] msg the string to be printed if is_expression_true is false (default value is "ERROR!")
 */
void AlwaysAssert(bool is_expression_true, const std::string& msg);

}

namespace WarpXUtilStr
{
    /** Return true if elem is in vect, false otherwise
     * @param[in] vect vector of strings, typically names
     * @param[in] elem single string
     * @return true if elem is in vect, false otherwise
     */
    bool is_in(const std::vector<std::string>& vect,
               const std::string& elem);

    /** Return true if any element in elems is in vect, false otherwise
     * @param[in] vect vector of strings, typically names
     * @param[in] elems vector of string
     * @return true if any element in elems is in vect, false otherwise
     */
    bool is_in(const std::vector<std::string>& vect,
               const std::vector<std::string>& elems);
}

#endif //WARPX_UTILS_H_
