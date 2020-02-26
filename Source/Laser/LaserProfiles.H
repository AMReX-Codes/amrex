/* Copyright 2019-2020 Luca Fedeli, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LaserProfiles_H_
#define WARPX_LaserProfiles_H_

#include "Parser/WarpXParser.H"

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

#include <map>
#include <string>
#include <memory>
#include <functional>
#include <limits>
#include <utility>


namespace WarpXLaserProfiles {

/** Common laser profile parameters
 *
 * Parameters for each laser profile as shared among all laser profile classes.
 */
struct CommonLaserParameters
{
    amrex::Real wavelength; //! central wavelength
    amrex::Real e_max;  //! maximum electric field at peak
    amrex::Vector<amrex::Real> p_X;// ! Polarization
    amrex::Vector<amrex::Real> nvec; //! Normal of the plane of the antenna
};


/** Abstract interface for laser profile classes
 *
 * Each new laser profile should inherit this interface and implement three
 * methods: init, update and fill_amplitude (described below).
 *
 * The declaration of a LaserProfile class should be placed in this file,
 * while the implementation of the methods should be in a dedicated file in
 * LaserProfilesImpl folder. LaserProfile classes should appear in
 * laser_profiles_dictionary to be used by LaserParticleContainer.
 */
class ILaserProfile
{
public:
    /** Initialize Laser Profile
     *
     * Reads the section of the inputfile relative to the laser beam
     * (e.g. laser_name.profile_t_peak, laser_name.profile_duration...)
     * and the "my_constants" section. It also receives some common
     * laser profile parameters. It uses these data to initialize the
     * member variables of the laser profile class.
     *
     * @param[in] ppl should be amrex::ParmParse(laser_name)
     * @param[in] ppc should be amrex::ParmParse("my_constants")
     * @param[in] params common laser profile parameters
     */
    virtual void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) = 0;

    /** Update Laser Profile
     *
     * Some laser profiles might need to perform an "update" operation per
     * time step.
     *
     * @param[in] t Current physical time in the simulation (seconds)
     */
    virtual void
    update (
        amrex::Real t) = 0;

    /** Fill Electric Field Amplitude for each particle of the antenna.
     *
     * Xp, Yp and amplitude must be arrays with the same length
     *
     * @param[in] np number of antenna particles
     * @param[in] Xp X coordinate of the particles of the antenna
     * @param[in] Yp Y coordinate of the particles of the antenna
     * @param[in] t time (seconds)
     * @param[out] amplitude of the electric field (V/m)
     */
    virtual void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const = 0;

    virtual ~ILaserProfile(){};
};

/**
 * Gaussian laser profile
 */
class GaussianLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    //No update needed
    void
    update (amrex::Real /*t */) override final {}

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const override final;

private:
    struct {
        amrex::Real waist          = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real duration       = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real t_peak         = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real focal_distance = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real zeta           = 0;
        amrex::Real beta           = 0;
        amrex::Real phi2           = 0;

        amrex::Vector<amrex::Real> stc_direction; //! Direction of the spatio-temporal couplings
        amrex::Real theta_stc; //! Angle between polarization (p_X) and direction of spatiotemporal coupling (stc_direction)
    } m_params;

    CommonLaserParameters m_common_params;
};

/**
 * Harris laser profile
 */
class HarrisLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    //No update needed
    void
    update (amrex::Real /*t */) override final {}

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const override final;

private:
    struct {
        amrex::Real waist          = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real duration       = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real focal_distance = std::numeric_limits<amrex::Real>::quiet_NaN();
    } m_params;

    CommonLaserParameters m_common_params;
};

/**
 * Laser profile defined by the used with an analytical expression
 */
class FieldFunctionLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    //No update needed
    void
    update (amrex::Real /*t */) override final {}

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const override final;

private:
    struct{
        std::string field_function;
    } m_params;

     WarpXParser m_parser;
};

/**
 * Laser profile read from an 'TXYE' file
 * The binary file must contain:
 * - 3 unsigned integers (4 bytes): nt (points along t), nx (points along x) and ny (points along y)
 * - nt*nx*ny doubles (8 bytes) in row major order : field amplitude
 */
class FromTXYEFileLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    /** \brief Reads new field data chunk from file if needed
    *
    * @param[in] t simulation time (seconds)
    */
    void
    update (amrex::Real t) override final;

    /** \brief compute field amplitude at particles' position for a laser beam
    * loaded from an E(x,y,t) file.
    *
    * Both Xp and Yp are given in laser plane coordinate.
    * For each particle with position Xp and Yp, this routine computes the
    * amplitude of the laser electric field, stored in array amplitude.
    *
    * \param np: number of laser particles
    * \param Xp: pointer to first component of positions of laser particles
    * \param Yp: pointer to second component of positions of laser particles
    * \param t: Current physical time
    * \param amplitude: pointer to array of field amplitude.
    */
    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const override final;

    /** \brief Function to fill the amplitude in case of a uniform grid.
    * This function cannot be private due to restrictions related to
    * the use of extended __device__ lambda
    *
    * \param idx_t_left index of the last time coordinate < t
    * \param np: number of laser particles
    * \param Xp: pointer to first component of positions of laser particles
    * \param Yp: pointer to second component of positions of laser particles
    * \param t: Current physical time
    * \param amplitude: pointer to array of field amplitude.
    */
    void internal_fill_amplitude_uniform(
        const int idx_t_left,
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const;

    /** \brief Function to fill the amplitude in case of a non-uniform grid.
    * This function cannot be private due to restrictions related to
    * the use of extended __device__ lambda
    *
    * \param idx_t_left index of the last time coordinate < t
    * \param np: number of laser particles
    * \param Xp: pointer to first component of positions of laser particles
    * \param Yp: pointer to second component of positions of laser particles
    * \param t: Current physical time
    * \param amplitude: pointer to array of field amplitude.
    */
    void internal_fill_amplitude_nonuniform(
        const int idx_t_left,
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const;

private:
    /** \brief parse a field file in the binary 'txye' format (whose details are given below).
    *
    * A 'txye' file should be a binary file with the following format:
    * -flag to indicate if the grid is uniform or not (1 byte, 0 means non-uniform, !=0 means uniform)
    * -np, number of timesteps (uint32_t, must be >=2)
    * -nx, number of points along x (uint32_t, must be >=2)
    * -ny, number of points along y (uint32_t, must be 1 for 2D simulations and >=2 for 3D simulations)
    * -timesteps (double[2] if grid is uniform, double[np] otherwise)
    * -x_coords (double[2] if grid is uniform, double[nx] otherwise)
    * -y_coords (double[1] if 2D, double[2] if 3D & uniform grid, double[ny] if 3D & non uniform grid)
    * -field_data (double[nt * nx * ny], with nt being the slowest coordinate).
    * The spatiotemporal grid must be rectangular.
    *
    * \param txye_file_name: name of the file to parse
    */
    void parse_txye_file(std::string txye_file_name);

    /** \brief Finds left and right time indices corresponding to time t
    *
    *
    * \param t: simulation time
    */
    std::pair<int,int> find_left_right_time_indices(amrex::Real t) const;

    /** \brief Load field data within the temporal range [t_begin, t_end)
    *
    * Must be called after having parsed a data file with parse_txye_file.
    *
    * \param t_begin: left limit of the timestep range to read
    * \param t_end: right limit of the timestep range to read (t_end is not read)
    */
    void read_data_t_chuck(int t_begin, int t_end);

    /**
     * \brief m_params contains all the internal parameters
     * used by this laser profile
     */
    struct{
        /** Name of the file containing the data */
        std::string txye_file_name;
        /** Flag to store if the grid is uniform */
        bool is_grid_uniform = false;
        /** Dimensions of E_data. nt, nx must be >=2.
         * If DIM=3, ny must be >=2 as well.
         * If DIM=2, ny must be 1 */
        int nt, nx, ny;
        /** Vector of temporal coordinates. For a non-uniform grid, it contains
         * all values of time. For a uniform grid, it contains only the start and stop
         * times and intermediate times are obtained with nt */
        amrex::Gpu::ManagedVector<amrex::Real> t_coords;
        /** Vector or x coordinates. For a non-uniform grid, it contains all values
         * of space dimension x. For a uniform grid, it contains only the min and max
         * x coordinates, and intermediate positions are obtained with nx */
        amrex::Gpu::ManagedVector<amrex::Real> x_coords;
        /** Vector or y coordinates. For a non-uniform grid, it contains all values
         * of space dimension y. For a uniform grid, it contains only the min and max
         * y coordinates, and intermediate positions are obtained with ny */
        amrex::Gpu::ManagedVector<amrex::Real> y_coords;
        /** Size of the timestep range to load */
        int time_chunk_size;
        /** Index of the first timestep in memory */
        int first_time_index;
        /** Index of the last timestep in memory */
        int last_time_index;
        /** Field data */
        amrex::Gpu::ManagedVector<amrex::Real> E_data;
    } m_params;

    CommonLaserParameters m_common_params;
};

/**
 * Maps laser profile names to lambdas returing unique pointers
 * to the corresponding laser profile objects.
 */
const
std::map<
std::string,
std::function<std::unique_ptr<ILaserProfile>()>
>
laser_profiles_dictionary =
{
    {"gaussian",
        [] () {return std::make_unique<GaussianLaserProfile>();} },
    {"harris",
        [] () {return std::make_unique<HarrisLaserProfile>();} },
    {"parse_field_function",
        [] () {return std::make_unique<FieldFunctionLaserProfile>();} },
    {"from_txye_file",
        [] () {return std::make_unique<FromTXYEFileLaserProfile>();} }
};

} //WarpXLaserProfiles

#endif //WARPX_LaserProfiles_H_
