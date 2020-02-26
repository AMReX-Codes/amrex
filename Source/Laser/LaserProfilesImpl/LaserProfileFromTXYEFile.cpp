/* Copyright 2019-2020 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"
#include "Utils/WarpX_Complex.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"

#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <algorithm>


using namespace amrex;

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::init (
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
            temp, m_params.time_chunk_size);
    }
    if(m_params.time_chunk_size < 2){
        Abort("Error! time_chunk_size must be >= 2!");
    }

    //Allocate memory for E_data Vector
    const int data_size = m_params.time_chunk_size*
            m_params.nx*m_params.ny;
    m_params.E_data = Gpu::ManagedVector<amrex::Real>(data_size);

    //Read first time chunck
    read_data_t_chuck(0, m_params.time_chunk_size);

    //Copy common params
    m_common_params = params;
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::update (amrex::Real t)
{
    if(t >= m_params.t_coords.back())
        return;

    const auto idx_times = find_left_right_time_indices(t);
    const auto idx_t_left = idx_times.first;
    const auto idx_t_right = idx_times.second;

    //Load data chunck if needed
    if(idx_t_right >  m_params.last_time_index){
        read_data_t_chuck(idx_t_left, idx_t_left+m_params.time_chunk_size);
    }
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::fill_amplitude (
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    //Amplitude is 0 if time is out of range
    if(t < m_params.t_coords.front() ||  t > m_params.t_coords.back()){
        amrex::ParallelFor(np,
            [=] AMREX_GPU_DEVICE (int i) {
                amplitude[i] = 0.0_rt;});
        return;
    }

    //Find left and right time indices
    int idx_t_left, idx_t_right;
    std::tie(idx_t_left, idx_t_right) = find_left_right_time_indices(t);

    if(idx_t_left <  m_params.first_time_index){
        Abort("Something bad has happened with the simulation time");
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

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::parse_txye_file(std::string txye_file_name)
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
        // Convert from double to amrex::Real
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
    //When a non-uniform grid is used, nt, nx and ny are identical
    //to t_coords.size(), x_coords.size() and y_coords.size().
    //When a uniform grid is used, nt,nx and ny store the number of points
    //used for the mesh, while t_coords, x_coords and y_coords store the
    //extrems in each direaction. Thus t_coords and x_coords in this case
    //have size 2 and y_coords has size 1 in 2D and size 2 in 3D.
    int t_sizes[6] = {m_params.nt, m_params.nx, m_params.ny,
        static_cast<int>(m_params.t_coords.size()),
        static_cast<int>(m_params.x_coords.size()),
        static_cast<int>(m_params.y_coords.size())};
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

std::pair<int,int>
WarpXLaserProfiles::FromTXYEFileLaserProfile::find_left_right_time_indices(amrex::Real t) const
{
    int idx_t_right;
    if(m_params.is_grid_uniform){
        const auto t_min = m_params.t_coords.front();
        const auto t_max = m_params.t_coords.back();
        const auto temp_idx_t_right = static_cast<int>(
            ceil( (m_params.nt-1)*(t-t_min)/(t_max-t_min)));
        idx_t_right = max(min(temp_idx_t_right, m_params.nt-1),1);
    }
    else{
        idx_t_right = std::distance(m_params.t_coords.begin(),
        std::upper_bound(m_params.t_coords.begin(),
            m_params.t_coords.end(), t));
    }
    return std::make_pair(idx_t_right-1, idx_t_right);
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::read_data_t_chuck(int t_begin, int t_end)
{
    amrex::Print() <<
        "Reading [" << t_begin << ", " << t_end <<
        ") data chunk from " << m_params.txye_file_name << "\n";

    //Indices of the first and last timestep to read
    auto i_first = max(0, t_begin);
    auto i_last = min(t_end-1, m_params.nt-1);
    if(i_last-i_first+1 > m_params.E_data.size())
        Abort("Data chunk to read from file is too large");

    if(ParallelDescriptor::IOProcessor()){
        //Read data chunk
        std::ifstream inp(m_params.txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");
        auto skip_amount = 1 +
            3*sizeof(uint32_t) +
            m_params.t_coords.size()*sizeof(double) +
            m_params.x_coords.size()*sizeof(double) +
            m_params.y_coords.size()*sizeof(double) +
            sizeof(double)*t_begin*m_params.nx*m_params.ny;
        inp.ignore(skip_amount);
        if(!inp) Abort("Failed to read field data from txye file");
        const int read_size = (i_last - i_first + 1)*
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

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::internal_fill_amplitude_uniform(
    const int idx_t_left,
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU.
    const auto tmp_e_max = m_common_params.e_max;
    const auto tmp_x_min = m_params.x_coords.front();
    const auto tmp_x_max = m_params.x_coords.back();
    const auto tmp_y_min = m_params.y_coords.front();
    const auto tmp_y_max = m_params.y_coords.back();
    const auto tmp_nx = m_params.nx;
#if (AMREX_SPACEDIM == 3)
    const auto tmp_ny = m_params.ny;
#endif
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const int idx_t_right = idx_t_left+1;
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
        const int temp_idx_x_right = static_cast<int>(
            ceil((tmp_nx-1)*(Xp[i]- tmp_x_min)/(tmp_x_max-tmp_x_min)));
        const int idx_x_right =
            max(min(temp_idx_x_right,tmp_nx-1),static_cast<int>(1));
        const int idx_x_left = idx_x_right - 1;
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
        const int temp_idx_y_right = static_cast<int>(
            ceil((tmp_ny-1)*(Yp[i]- tmp_y_min)/(tmp_y_max-tmp_y_min)));
        const int idx_y_right =
            max(min(temp_idx_y_right,tmp_ny-1),static_cast<int>(1));
        const int idx_y_left = idx_y_right - 1;
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

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::internal_fill_amplitude_nonuniform(
    const int idx_t_left,
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU.
    const auto tmp_e_max = m_common_params.e_max;
    const auto tmp_x_min = m_params.x_coords.front();
    const auto tmp_x_max = m_params.x_coords.back();
    const auto tmp_y_min = m_params.y_coords.front();
    const auto tmp_y_max = m_params.y_coords.back();
    const auto p_x_coords = m_params.x_coords.dataPtr();
    const int tmp_x_coords_size = static_cast<int>(m_params.x_coords.size());
    const auto p_y_coords = m_params.y_coords.dataPtr();
    const int tmp_y_coords_size = static_cast<int>(m_params.y_coords.size());
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const int idx_t_right = idx_t_left+1;
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
        const int idx_x_right = p_x_right - p_x_coords;
        const int idx_x_left = idx_x_right - 1;

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
        const int idx_y_right = p_y_right - p_y_coords;
        const int idx_y_left = idx_y_right - 1;

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
