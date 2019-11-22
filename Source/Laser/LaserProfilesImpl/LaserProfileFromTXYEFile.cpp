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

    // Parse the name of the XYTE file
    ppl.get("txye_file_name", m_params.txye_file_name);
    if(m_params.txye_file_name.empty())
    {
        Abort("txye_file_name must be provided for txye_file laser profile!");
    }

    parse_txye_file(m_params.txye_file_name);

    m_params.time_chunk_size = m_params.t_coords.size();
    int temp = 1;
    if(ppl.query("time_chunk_size", temp)){
        m_params.time_chunk_size = min(
            static_cast<size_t>(temp), m_params.time_chunk_size);
    }
    if(m_params.time_chunk_size < 2){
        Abort("Error! time_chunk_size must be >= 2!");
    }

    //Read first time chunck
    read_data_t_chuck(0, m_params.time_chunk_size);

    //Copy common params
    m_common_params = params;
}

/* \brief compute [...]
 *
 */
void
FromTXYEFileLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
{
    if(t < m_params.t_coords.front() ||  t > m_params.t_coords.back()){
        amrex::ParallelFor(np,
            [=] AMREX_GPU_DEVICE (int i) {
                amplitude[i] = 0.0_rt;});
        return;
    }

    const size_t idx_t_right = std::distance(m_params.t_coords.begin(),
        std::upper_bound(m_params.t_coords.begin(),
            m_params.t_coords.end(), t));
    const size_t idx_t_left = idx_t_right - 1;

    if(idx_t_left <  m_params.first_time_index){
        Abort("Something bad has happened with the simulation time");
    }

    #pragma omp critical
    {
        if(idx_t_right >  m_params.last_time_index){
            read_data_t_chuck(idx_t_left, idx_t_left+m_params.time_chunk_size);
        }
    }

    //REPLACE
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU runs.
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
    const auto tmp_t_left = m_params.t_coords[idx_t_left];
    const auto tmp_t_right = m_params.t_coords[idx_t_right];
    const auto tmp_idx_first_time = m_params.first_time_index;
    // Loop through the macroparticle to calculate the proper amplitude

    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
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
            auto const p_x_right = WarpXUtilAlgo::upper_bound(
                    p_x_coords, p_x_coords+tmp_x_coords_size, Xp[i]);
            const size_t idx_x_right = p_x_right - p_x_coords;
            const size_t idx_x_left = idx_x_right - 1;

#if (AMREX_SPACEDIM == 2)
            auto idx = [=](int i, int j){
                return (i-tmp_idx_first_time) * tmp_x_coords_size + j;
            };
            amplitude[i] = WarpXUtilAlgo::bilinear_interp(
                tmp_t_left, tmp_t_right,
                p_x_coords[idx_x_left], p_x_coords[idx_x_right],
                p_E_data[idx(idx_t_left, idx_x_left)],
                p_E_data[idx(idx_t_left, idx_x_right)],
                p_E_data[idx(idx_t_right, idx_x_left)],
                p_E_data[idx(idx_t_right, idx_x_right)],
                t, Xp[i])*tmp_e_max;

#elif (AMREX_SPACEDIM == 3)
            auto const p_y_right = WarpXUtilAlgo::upper_bound(
                p_y_coords, p_y_coords+tmp_y_coords_size, Yp[i]);
            const size_t idx_y_right = p_y_right - p_y_coords;
            const size_t idx_y_left = idx_y_right - 1;
            auto idx = [=](int i, int j, int k){
                return
                    (i-tmp_idx_first_time)*tmp_x_coords_size*tmp_y_coords_size+
                    j*tmp_y_coords_size + k;
            };
            amplitude[i] = WarpXUtilAlgo::trilinear_interp(
                tmp_t_left, tmp_t_right,
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

void
FromTXYEFileLaserProfile::parse_txye_file(std::string txye_file_name)
{
    size_t nt, nx, ny;
    if(ParallelDescriptor::IOProcessor()){
        std::ifstream inp(txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");

        auto const three_uint32_size = sizeof(uint32_t)*3;
        char buf[three_uint32_size];
        inp.read(buf, three_uint32_size);
        if(!inp) Abort("Failed to read sizes from txye file");
        nt = reinterpret_cast<uint32_t*>(buf)[0];
        nx = reinterpret_cast<uint32_t*>(buf)[1];
        ny = reinterpret_cast<uint32_t*>(buf)[2];

        Vector<double> buf_t(nt);
        Vector<double> buf_x(nx);
        Vector<double> buf_y(ny);

        inp.read(reinterpret_cast<char*>(buf_t.dataPtr()), nt*sizeof(double));
        if(!inp)
            Abort("Failed to read coords from txye file");
        if (!std::is_sorted(buf_t.begin(), buf_t.end()))
            Abort("Coordinates are not sorted  in txye file");
        inp.read(reinterpret_cast<char*>(buf_x.dataPtr()), nx*sizeof(double));
        if(!inp)
            Abort("Failed to read coords from txye file");
        if (!std::is_sorted(buf_t.begin(), buf_t.end()))
            Abort("Coordinates are not sorted  in txye file");
        inp.read(reinterpret_cast<char*>(buf_y.dataPtr()), ny*sizeof(double));
        if(!inp)
            Abort("Failed to read coords from txye file");
        if (!std::is_sorted(buf_t.begin(), buf_t.end()))
            Abort("Coordinates are not sorted in txye file");

#if (AMREX_SPACEDIM == 3)
        if(ny <= 1) Abort("ny in txye file must be >=2 in 3D");
#endif

#if (AMREX_SPACEDIM == 2)
        if(ny != 1) Abort("ny in txye file must be 1 in 2D");
#endif

        m_params.t_coords = Gpu::ManagedVector<amrex::Real>(nt);
        m_params.x_coords = Gpu::ManagedVector<amrex::Real>(nx);
        m_params.y_coords = Gpu::ManagedVector<amrex::Real>(ny);

        std::transform(buf_t.begin(), buf_t.end(), m_params.t_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
        std::transform(buf_x.begin(), buf_x.end(), m_params.x_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
        std::transform(buf_y.begin(), buf_y.end(), m_params.y_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
    }

    size_t t_sizes[3] = {nt, nx, ny};
    ParallelDescriptor::Bcast(t_sizes, 3,
        ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();
    nt = t_sizes[0]; nx = t_sizes[1]; ny = t_sizes[2];

    if(!ParallelDescriptor::IOProcessor()){
        m_params.t_coords = Gpu::ManagedVector<amrex::Real>(nt);
        m_params.x_coords = Gpu::ManagedVector<amrex::Real>(nx);
        m_params.y_coords = Gpu::ManagedVector<amrex::Real>(ny);
    }

    ParallelDescriptor::Bcast(m_params.t_coords.dataPtr(),
        m_params.t_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(m_params.x_coords.dataPtr(),
        m_params.x_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(m_params.y_coords.dataPtr(),
        m_params.y_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();
}

void
FromTXYEFileLaserProfile::read_data_t_chuck(size_t t_begin, size_t t_end)
{
    auto i_first = max(static_cast<size_t>(0), t_begin);
    auto i_last = min(t_end-1, m_params.t_coords.size()-1);

    size_t data_size = (i_last - i_first + 1)*
            m_params.x_coords.size()*m_params.y_coords.size();

    m_params.E_data = Gpu::ManagedVector<amrex::Real>(data_size);

    if(ParallelDescriptor::IOProcessor()){
        std::ifstream inp(m_params.txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");

        size_t skip_amount = 3*sizeof(uint32_t) +
            m_params.t_coords.size()*sizeof(double) +
            m_params.x_coords.size()*sizeof(double) +
            m_params.y_coords.size()*sizeof(double) +
            t_begin*m_params.x_coords.size()*m_params.y_coords.size()*sizeof(double);

        inp.ignore(skip_amount);
        if(!inp) Abort("Failed to read field data from txye file");

        Vector<double> buf_e(data_size);

        inp.read(reinterpret_cast<char*>(buf_e.dataPtr()), data_size*sizeof(double));
        if(!inp) Abort("Failed to read field data from txye file");

        std::transform(buf_e.begin(), buf_e.end(), m_params.E_data.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
    }

    ParallelDescriptor::Bcast(m_params.E_data.dataPtr(),
        m_params.E_data.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Barrier();

    m_params.first_time_index = i_first;
    m_params.last_time_index = i_last;
}
