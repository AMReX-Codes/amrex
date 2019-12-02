#include <LaserProfiles.H>

#include <WarpX_Complex.H>
#include <WarpXConst.H>
#include <WarpXUtil.H>

#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <algorithm>

using namespace amrex;
using namespace WarpXLaserProfiles;

void
FromTXYEFileLaserProfile::init (
    const amrex::ParmParse& ppl,
    const amrex::ParmParse& /* ppc */,
    CommonLaserParameters params)
{
    if (!std::numeric_limits< double >::is_iec559)
    {
        Print() << R"(Warning: double does not comply with IEEE 754: bad
            things will happen parsing the X, Y and T profiles for the laser!)";
    }

    // Parse the TXYE file
    ppl.get("txye_file_name", m_params.txye_file_name);
    if(m_params.txye_file_name.empty())
    {
        Abort("txye_file_name must be provided for txye_file laser profile!");
    }
    parse_txye_file(m_params.txye_file_name);

    //Set time_chunk_size
    m_params.time_chunk_size = m_params.nt;
    int temp = 1;
    if(ppl.query("time_chunk_size", temp)){
        m_params.time_chunk_size = min(
            static_cast<size_t>(temp), m_params.time_chunk_size);
    }
    if(m_params.time_chunk_size < 2){
        Abort("Error! time_chunk_size must be >= 2!");
    }

    //Allocate memory for E_data Vector
    size_t data_size = m_params.time_chunk_size*
            m_params.nx*m_params.ny;
    m_params.E_data = Gpu::ManagedVector<amrex::Real>(data_size);

    //Read first time chunck
    read_data_t_chuck(0, m_params.time_chunk_size);

    //Copy common params
    m_common_params = params;
}

/* \brief compute field amplitude at particles' position for a laser beam
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
FromTXYEFileLaserProfile::fill_amplitude (
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
{
    //Amplitude is 0 if time is out of range
    if(t < m_params.t_coords.front() ||  t > m_params.t_coords.back()){
        amrex::ParallelFor(np,
            [=] AMREX_GPU_DEVICE (int i) {
                amplitude[i] = 0.0_rt;});
        return;
    }

    //Find left and right time indices
    size_t idx_t_right;
    if(m_params.is_grid_uniform){
        const auto t_min = m_params.t_coords.front();
        const auto t_max = m_params.t_coords.back();
        idx_t_right = static_cast<size_t>(
            ceil( (m_params.nt-1)*(t-t_min)/(t_max-t_min)));
        idx_t_right = max(min(idx_t_right, m_params.nt-1),static_cast<size_t>(1));
    }
    else{
        idx_t_right = std::distance(m_params.t_coords.begin(),
        std::upper_bound(m_params.t_coords.begin(),
            m_params.t_coords.end(), t));
    }
    const size_t idx_t_left = idx_t_right-1;
    if(idx_t_left <  m_params.first_time_index){
        Abort("Something bad has happened with the simulation time");
    }

    //Load data chunck if needed
    #pragma omp critical
    {
        if(idx_t_right >  m_params.last_time_index){
            read_data_t_chuck(idx_t_left, idx_t_left+m_params.time_chunk_size);
        }
    }

    if(m_params.is_grid_uniform){
        internal_fill_amplitude_uniform(
            idx_t_left, np, Xp, Yp, t, amplitude);
    }
    else{
        internal_fill_amplitude_nonuniform(
            idx_t_left, np, Xp, Yp, t, amplitude);
    }
}

/* \brief parse a field file in the binary 'txye' format (whose details are given below).
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
void
FromTXYEFileLaserProfile::parse_txye_file(std::string txye_file_name)
{
    if(ParallelDescriptor::IOProcessor()){
        std::ifstream inp(txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");

        //Uniform grid flag
        char flag;
        inp.read(&flag, 1);
        if(!inp) Abort("Failed to read sizes from txye file");
        m_params.is_grid_uniform=flag;

        //Grid points along t, x and y
        auto const three_uint32_size = sizeof(uint32_t)*3;
        char buf[three_uint32_size];
        inp.read(buf, three_uint32_size);
        if(!inp) Abort("Failed to read sizes from txye file");
        m_params.nt = reinterpret_cast<uint32_t*>(buf)[0];
        m_params.nx = reinterpret_cast<uint32_t*>(buf)[1];
        m_params.ny = reinterpret_cast<uint32_t*>(buf)[2];
        if(m_params.nt <= 1) Abort("nt in txye file must be >=2");
        if(m_params.nx <= 1) Abort("nx in txye file must be >=2");
#if (AMREX_SPACEDIM == 3)
        if(m_params.ny <= 1) Abort("ny in txye file must be >=2 in 3D");
#elif(AMREX_SPACEDIM == 2)
        if(m_params.ny != 1) Abort("ny in txye file must be 1 in 2D");
#endif

        //Coordinates
        Vector<double> buf_t, buf_x, buf_y;
        if(m_params.is_grid_uniform){
            buf_t.resize(2);
            buf_x.resize(2);
#if (AMREX_SPACEDIM == 3)
            buf_y.resize(2);
#elif(AMREX_SPACEDIM == 2)
            buf_y.resize(1);
#endif
        }
        else{
            buf_t.resize(m_params.nt);
            buf_x.resize(m_params.nx);
            buf_y.resize(m_params.ny);
        }
        inp.read(reinterpret_cast<char*>(buf_t.dataPtr()),
            buf_t.size()*sizeof(double));
        if(!inp)
            Abort("Failed to read coords from txye file");
        if (!std::is_sorted(buf_t.begin(), buf_t.end()))
            Abort("Coordinates are not sorted  in txye file");
        inp.read(reinterpret_cast<char*>(buf_x.dataPtr()),
            buf_x.size()*sizeof(double));
        if(!inp)
            Abort("Failed to read coords from txye file");
        if (!std::is_sorted(buf_x.begin(), buf_x.end()))
            Abort("Coordinates are not sorted  in txye file");
        inp.read(reinterpret_cast<char*>(buf_y.dataPtr()),
            buf_y.size()*sizeof(double));
        if(!inp)
            Abort("Failed to read coords from txye file");
        if (!std::is_sorted(buf_y.begin(), buf_y.end()))
            Abort("Coordinates are not sorted in txye file");
        m_params.t_coords = Gpu::ManagedVector<amrex::Real>(buf_t.size());
        m_params.x_coords = Gpu::ManagedVector<amrex::Real>(buf_x.size());
        m_params.y_coords = Gpu::ManagedVector<amrex::Real>(buf_y.size());
        std::transform(buf_t.begin(), buf_t.end(), m_params.t_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
        std::transform(buf_x.begin(), buf_x.end(), m_params.x_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
        std::transform(buf_y.begin(), buf_y.end(), m_params.y_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
    }

    //Broadcast grid uniformity
    char is_grid_uniform = m_params.is_grid_uniform;
    ParallelDescriptor::Bcast(&is_grid_uniform, 1,
        ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();
    m_params.is_grid_uniform = is_grid_uniform;

    //Broadcast grid size and coordinate sizes
    size_t t_sizes[6] = {m_params.nt, m_params.nx, m_params.ny,
        m_params.t_coords.size(),
        m_params.x_coords.size(),
        m_params.y_coords.size()};
    ParallelDescriptor::Bcast(t_sizes, 6,
        ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();
    m_params.nt = t_sizes[0]; m_params.nx = t_sizes[1]; m_params.ny = t_sizes[2];

    //Broadcast coordinates
    if(!ParallelDescriptor::IOProcessor()){
        m_params.t_coords = Gpu::ManagedVector<amrex::Real>(t_sizes[3]);
        m_params.x_coords = Gpu::ManagedVector<amrex::Real>(t_sizes[4]);
        m_params.y_coords = Gpu::ManagedVector<amrex::Real>(t_sizes[5]);
    }
    ParallelDescriptor::Bcast(m_params.t_coords.dataPtr(),
        m_params.t_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(m_params.x_coords.dataPtr(),
        m_params.x_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(m_params.y_coords.dataPtr(),
        m_params.y_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();
}

/* \brief Load field data within the temporal range [t_begin, t_end)
 *
 * Must be called after having parsed a data file with parse_txye_file.
 *
 * \param t_begin: left limit of the timestep range to read
 * \param t_end: right limit of the timestep range to read (t_end is not read)
 */
void
FromTXYEFileLaserProfile::read_data_t_chuck(size_t t_begin, size_t t_end)
{
    amrex::Print() <<
        "Reading [" << t_begin << ", " << t_end <<
        ") data chunk from " << m_params.txye_file_name << "\n";

    //Indices of the first and last timestep to read
    auto i_first = max(static_cast<size_t>(0), t_begin);
    auto i_last = min(t_end-1, m_params.nt-1);
    if(i_last-i_first+1 > m_params.E_data.size())
        Abort("Data chunk to read from file is too large");

    if(ParallelDescriptor::IOProcessor()){
        //Read data chunk
        std::ifstream inp(m_params.txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");
        size_t skip_amount = 1 +
            3*sizeof(uint32_t) +
            m_params.t_coords.size()*sizeof(double) +
            m_params.x_coords.size()*sizeof(double) +
            m_params.y_coords.size()*sizeof(double) +
            t_begin*m_params.nx*m_params.ny*sizeof(double);
        inp.ignore(skip_amount);
        if(!inp) Abort("Failed to read field data from txye file");
        const size_t read_size = (i_last - i_first + 1)*
            m_params.nx*m_params.ny;
        Vector<double> buf_e(read_size);
        inp.read(reinterpret_cast<char*>(buf_e.dataPtr()), read_size*sizeof(double));
        if(!inp) Abort("Failed to read field data from txye file");
        std::transform(buf_e.begin(), buf_e.end(), m_params.E_data.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
    }

    //Broadcast E_data
    ParallelDescriptor::Bcast(m_params.E_data.dataPtr(),
        m_params.E_data.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();

    //Update first and last indices
    m_params.first_time_index = i_first;
    m_params.last_time_index = i_last;
}

/* \brief Private function to fill the amplitude in case of a uniform grid
 *
 * \param idx_t_left index of the last time coordinate < t
 * \param np: number of laser particles
 * \param Xp: pointer to first component of positions of laser particles
 * \param Yp: pointer to second component of positions of laser particles
 * \param t: Current physical time
 * \param amplitude: pointer to array of field amplitude.
 */
void
FromTXYEFileLaserProfile::internal_fill_amplitude_uniform(
    const size_t idx_t_left,
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
{
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU.
    const auto tmp_e_max = m_common_params.e_max;
    const auto tmp_x_min = m_params.x_coords.front();
    const auto tmp_x_max = m_params.x_coords.back();
    const auto tmp_y_min = m_params.y_coords.front();
    const auto tmp_y_max = m_params.y_coords.back();
    const auto tmp_nx = m_params.nx;
    const auto tmp_ny = m_params.ny;
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const size_t idx_t_right = idx_t_left+1;
    const auto t_left = idx_t_left*
        (m_params.t_coords.back()-m_params.t_coords.front())/(m_params.nt-1) +
        m_params.t_coords.front();
    const auto t_right = idx_t_right*
        (m_params.t_coords.back()-m_params.t_coords.front())/(m_params.nt-1) +
        m_params.t_coords.front();

    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
    np,
    [=] AMREX_GPU_DEVICE (int i) {
        //Amplitude is zero if we are out of bounds
        if (Xp[i] <= tmp_x_min || Xp[i] >= tmp_x_max){
            amplitude[i] = 0.0_rt;
            return;
        }
#if (AMREX_SPACEDIM == 3)
        if (Yp[i] <= tmp_y_min || Yp[i] >= tmp_y_max){
            amplitude[i] = 0.0_rt;
            return;
        }
#endif
        //Find indices and coordinates along x
        const size_t temp_idx_x_right = static_cast<size_t>(
            ceil((tmp_nx-1)*(Xp[i]- tmp_x_min)/(tmp_x_max-tmp_x_min)));
        const size_t idx_x_right =
            max(min(temp_idx_x_right,tmp_nx-1),static_cast<size_t>(1));
        const size_t idx_x_left = idx_x_right - 1;
        const auto x_0 =
            idx_x_left*(tmp_x_max-tmp_x_min)/(tmp_nx-1) + tmp_x_min;
        const auto x_1 =
            idx_x_right*(tmp_x_max-tmp_x_min)/(tmp_nx-1) + tmp_x_min;

#if (AMREX_SPACEDIM == 2)
        //Interpolate amplitude
        const auto idx = [=](int i, int j){
            return (i-tmp_idx_first_time) * tmp_nx + j;
        };
        amplitude[i] = WarpXUtilAlgo::bilinear_interp(
            t_left, t_right,
            x_0, x_1,
            p_E_data[idx(idx_t_left, idx_x_left)],
            p_E_data[idx(idx_t_left, idx_x_right)],
            p_E_data[idx(idx_t_right, idx_x_left)],
            p_E_data[idx(idx_t_right, idx_x_right)],
            t, Xp[i])*tmp_e_max;

#elif (AMREX_SPACEDIM == 3)
        //Find indices and coordinates along y
        const size_t temp_idx_y_right = static_cast<size_t>(
            ceil((tmp_ny-1)*(Yp[i]- tmp_y_min)/(tmp_y_max-tmp_y_min)));
        const size_t idx_y_right =
            max(min(temp_idx_y_right,tmp_ny-1),static_cast<size_t>(1));
        const size_t idx_y_left = idx_y_right - 1;
        const auto y_0 =
            idx_y_left*(tmp_y_max-tmp_y_min)/(tmp_ny-1) + tmp_y_min;
        const auto y_1 =
            idx_y_right*(tmp_y_max-tmp_y_min)/(tmp_ny-1) + tmp_y_min;

        //Interpolate amplitude
        const auto idx = [=](int i, int j, int k){
            return
                (i-tmp_idx_first_time)*tmp_nx*tmp_ny+
                j*tmp_ny + k;
        };
        amplitude[i] = WarpXUtilAlgo::trilinear_interp(
            t_left, t_right,
            x_0, x_1,
            y_0, y_1,
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_right)],
            t, Xp[i], Yp[i])*tmp_e_max;
#endif
        }
    );
}

/* \brief Private function to fill the amplitude in case of a non-uniform grid
 *
 * \param idx_t_left index of the last time coordinate < t
 * \param np: number of laser particles
 * \param Xp: pointer to first component of positions of laser particles
 * \param Yp: pointer to second component of positions of laser particles
 * \param t: Current physical time
 * \param amplitude: pointer to array of field amplitude.
 */
void
FromTXYEFileLaserProfile::internal_fill_amplitude_nonuniform(
    const size_t idx_t_left,
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
{
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU.
    const auto tmp_e_max = m_common_params.e_max;
    const auto tmp_x_min = m_params.x_coords.front();
    const auto tmp_x_max = m_params.x_coords.back();
    const auto tmp_y_min = m_params.y_coords.front();
    const auto tmp_y_max = m_params.y_coords.back();
    const auto p_x_coords = m_params.x_coords.dataPtr();
    const auto tmp_x_coords_size = m_params.x_coords.size();
    const auto p_y_coords = m_params.y_coords.dataPtr();
    const auto tmp_y_coords_size = m_params.y_coords.size();
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const size_t idx_t_right = idx_t_left+1;
    const auto t_left = m_params.t_coords[idx_t_left];
    const auto t_right = m_params.t_coords[idx_t_right];

    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
    np,
    [=] AMREX_GPU_DEVICE (int i) {
        //Amplitude is zero if we are out of bounds
        if (Xp[i] <= tmp_x_min || Xp[i] >= tmp_x_max){
            amplitude[i] = 0.0_rt;
            return;
        }
#if (AMREX_SPACEDIM == 3)
        if (Yp[i] <= tmp_y_min || Yp[i] >= tmp_y_max){
            amplitude[i] = 0.0_rt;
            return;
        }
#endif

        //Find indices along x
        auto const p_x_right = WarpXUtilAlgo::upper_bound(
                p_x_coords, p_x_coords+tmp_x_coords_size, Xp[i]);
        const size_t idx_x_right = p_x_right - p_x_coords;
        const size_t idx_x_left = idx_x_right - 1;

#if (AMREX_SPACEDIM == 2)
        //Interpolate amplitude
        const auto idx = [=](int i, int j){
            return (i-tmp_idx_first_time) * tmp_x_coords_size + j;
        };
        amplitude[i] = WarpXUtilAlgo::bilinear_interp(
            t_left, t_right,
            p_x_coords[idx_x_left], p_x_coords[idx_x_right],
            p_E_data[idx(idx_t_left, idx_x_left)],
            p_E_data[idx(idx_t_left, idx_x_right)],
            p_E_data[idx(idx_t_right, idx_x_left)],
            p_E_data[idx(idx_t_right, idx_x_right)],
            t, Xp[i])*tmp_e_max;

#elif (AMREX_SPACEDIM == 3)
        //Find indices along y
        auto const p_y_right = WarpXUtilAlgo::upper_bound(
            p_y_coords, p_y_coords+tmp_y_coords_size, Yp[i]);
        const size_t idx_y_right = p_y_right - p_y_coords;
        const size_t idx_y_left = idx_y_right - 1;

        //Interpolate amplitude
        const auto idx = [=](int i, int j, int k){
            return
                (i-tmp_idx_first_time)*tmp_x_coords_size*tmp_y_coords_size+
                j*tmp_y_coords_size + k;
        };
        amplitude[i] = WarpXUtilAlgo::trilinear_interp(
            t_left, t_right,
            p_x_coords[idx_x_left], p_x_coords[idx_x_right],
            p_y_coords[idx_y_left], p_y_coords[idx_y_right],
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_right)],
            t, Xp[i], Yp[i])*tmp_e_max;
#endif
        }
    );
}
