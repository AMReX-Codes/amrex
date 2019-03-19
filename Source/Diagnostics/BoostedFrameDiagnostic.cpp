#include <AMReX_MultiFabUtil.H>

#include "BoostedFrameDiagnostic.H"
#include "WarpX_f.H"
#include "WarpX.H"

using namespace amrex;

#ifdef WARPX_USE_HDF5

#include <hdf5.h>

/*
  Helper functions for doing the HDF5 IO.

 */
namespace
{
    const std::vector<std::string> particle_field_names = {"w", "x", "y", "z", "ux", "uy", "uz"};
    const std::vector<std::string> mesh_field_names =
        {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz", "rho"};
        
    /*
      Creates the HDF5 file in truncate mode and closes it.
      Should be run only by the root process.
    */
    void output_create(const std::string& file_path) {
        BL_PROFILE("output_create");
        hid_t file = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file < 0) {
            amrex::Abort("Error: could not create file at " + file_path);
        }
        H5Fclose(file);
    }

    /*
      Writes a single string attribute to the given group. 
      Should only be called by the root process.
    */
    void write_string_attribute(hid_t& group, const std::string& key, const std::string& val)
    {
        hid_t str_type     = H5Tcopy(H5T_C_S1);
        hid_t scalar_space = H5Screate(H5S_SCALAR);

        // Fix the str_type length for the format string.
        H5Tset_size(str_type, strlen(val.c_str()));

        hid_t attr = H5Acreate(group, key.c_str(), str_type, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_type, val.c_str());

        H5Aclose(attr);
        H5Sclose(scalar_space);
        H5Tclose(str_type);
    }

    /*
      Writes a single double attribute to the given group. 
      Should only be called by the root process.
    */
    void write_double_attribute(hid_t& group, const std::string& key, const double val)
    {
        hid_t scalar_space = H5Screate(H5S_SCALAR);

        hid_t attr = H5Acreate(group, key.c_str(), H5T_IEEE_F32LE, scalar_space,
                               H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, &val);

        H5Aclose(attr);
        H5Sclose(scalar_space);
    }

    /*
      Opens the output file and writes all of metadata attributes.
      Should be run only by the root process.
    */
    void output_write_metadata(const std::string& file_path,
                               const int istep, const Real time, const Real dt)
    {
        BL_PROFILE("output_write_metadata");
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

        write_string_attribute(file, "software", "warpx");
        write_string_attribute(file, "softwareVersion", "0.0.0");
        write_string_attribute(file, "meshesPath", "fields/");
        write_string_attribute(file, "iterationEncoding", "fileBased");
        write_string_attribute(file, "iterationFormat", "data%T.h5");
        write_string_attribute(file, "openPMD", "1.1.0");
        write_string_attribute(file, "basePath", "/data/%T/");

        hid_t group = H5Gcreate(file, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        group = H5Gcreate(group, std::to_string(istep).c_str(),
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        write_double_attribute(group, "time", time);
        write_double_attribute(group, "timeUnitSI", 1.0);
        write_double_attribute(group, "dt", dt);

        // Field groups
        group = H5Gcreate(group, "fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Close all resources.
        H5Gclose(group);
        H5Fclose(file);
        H5close();
    }

    /*
      Creates a dataset with the given cell dimensions, at the path
      "/native_fields/(field_name)".
      Should be run only by the master rank.
    */
    void output_create_field(const std::string& file_path, const std::string& field_path,
                             const unsigned nx, const unsigned ny, const unsigned nz)
    {        
        BL_PROFILE("output_create_field");

        // Open the output.
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        // Create a 3D, nx x ny x nz dataspace.
#if (AMREX_SPACEDIM == 3)
        hsize_t dims[3] = {nx, ny, nz};
#else
        hsize_t dims[3] = {nx, nz};
#endif
        hid_t grid_space = H5Screate_simple(AMREX_SPACEDIM, dims, NULL);
        
        // Create the dataset.
        hid_t dataset = H5Dcreate(file, field_path.c_str(), H5T_IEEE_F64LE,
                                  grid_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        if (dataset < 0)
        {
            amrex::Abort("Error: could not create dataset. H5 returned "
                         + std::to_string(dataset) + "\n");
        }

        // Close resources.
        H5Dclose(dataset);
        H5Sclose(grid_space);
        H5Fclose(file);
    }

    /*
      Creates a group associated with a single particle species.
      Should be run by all processes collectively.
    */    
    void output_create_species_group(const std::string& file_path, const std::string& species_name)
    {
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;
        int mpi_rank;
        MPI_Comm_rank(comm, &mpi_rank);

        // Create the file access prop list.
        hid_t pa_plist = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pa_plist, comm, info);
        
        // Open the output.
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, pa_plist);
        
        hid_t group = H5Gcreate(file, species_name.c_str(),
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group);
        H5Fclose(file);

    }

    /*
      Creates an extendible dataset, suitable for storing particle data.
      Should be run only by the master rank.
    */
    long output_resize_particle_field(const std::string& file_path, const std::string& field_path,
                                      const long num_to_add)
    {        
        BL_PROFILE("output_resize_particle_field");

        // Open the output.
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

        int rank;
        hsize_t dims[1];

        hid_t dataset = H5Dopen2 (file, field_path.c_str(), H5P_DEFAULT);
        hid_t filespace = H5Dget_space (dataset);
        rank = H5Sget_simple_extent_ndims (filespace);
        herr_t status = H5Sget_simple_extent_dims (filespace, dims, NULL);

        // set new size
        hsize_t new_size[1];
        new_size[0] = dims[0] + num_to_add;
        status = H5Dset_extent (dataset, new_size);
        
        if (status < 0)
        {
            amrex::Abort("Error: set extent filed on dataset "
                         + std::to_string(dataset) + "\n");
        }

        // Close resources.
        H5Sclose(filespace);
        H5Dclose(dataset);
        H5Fclose(file);

        return dims[0];
    }

    /*
      Creates an extendible dataset, suitable for storing particle data.
      Should be run on all ranks collectively.
    */
    void output_write_particle_field(const std::string& file_path, const std::string& field_path,
                                     const Real* data_ptr, const long count, const long index)
    {        
        BL_PROFILE("output_write_particle_field");

        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;
        int mpi_rank;
        MPI_Comm_rank(comm, &mpi_rank);

        // Create the file access prop list.
        hid_t pa_plist = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pa_plist, comm, info);
        
        // Open the output.
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, pa_plist);
        
        int RANK = 1;
        hsize_t offset[1];
        hsize_t dims[1];
        herr_t status;
        
        hid_t dataset = H5Dopen (file, field_path.c_str(), H5P_DEFAULT);

        // Make sure the dataset is there.
        if (dataset < 0)
        {
            amrex::Abort("Error on rank " + std::to_string(mpi_rank) +
                         ". Count not find dataset " + field_path + "\n");
        }
        
        hid_t filespace = H5Dget_space (dataset);
        
        offset[0] = index;
        dims[0] = count;

        // Create collective io prop list.
        hid_t collective_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(collective_plist, H5FD_MPIO_INDEPENDENT);

        if (count > 0) {
        
            /* Define memory space */
            hid_t memspace = H5Screate_simple (RANK, dims, NULL);
            
            status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
                                          dims, NULL);
            
            if (status < 0)
            {
                amrex::Abort("Error on rank " + std::to_string(ParallelDescriptor::MyProc()) +
                             " could not select hyperslab.\n");
            }
            
            /* Write the data to the extended portion of dataset  */
            status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace,
                              filespace, collective_plist, data_ptr);
        
            if (status < 0)
            {
                amrex::Abort("Error on rank " + std::to_string(ParallelDescriptor::MyProc()) +
                             " could not write hyperslab.\n");
            }

            status = H5Sclose (memspace);
        }
        
        ParallelDescriptor::Barrier();
        
        // Close resources.
        H5Pclose(collective_plist);
        H5Sclose(filespace);
        H5Dclose(dataset);
        H5Fclose(file);
        H5Pclose(pa_plist);
    }
    
    /*
      Creates an extendible dataset, suitable for storing particle data.
      Should be run on all ranks collectively.
    */
    void output_create_particle_field(const std::string& file_path, const std::string& field_path)
    {        
        BL_PROFILE("output_create_particle_field");

        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;
        int mpi_rank;
        MPI_Comm_rank(comm, &mpi_rank);

        // Create the file access prop list.
        hid_t pa_plist = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pa_plist, comm, info);
        
        // Open the output.
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, pa_plist);

        constexpr int RANK = 1;
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hsize_t      chunk_dims[2] = {4};
        
        hid_t dataspace = H5Screate_simple (RANK, dims, maxdims);

        // Enable chunking
        hid_t prop = H5Pcreate (H5P_DATASET_CREATE);
        herr_t status = H5Pset_chunk (prop, RANK, chunk_dims);

        hid_t dataset = H5Dcreate2 (file, field_path.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                                    H5P_DEFAULT, prop, H5P_DEFAULT);
        
        if (dataset < 0)
        {
            amrex::Abort("Error: could not create dataset. H5 returned "
                         + std::to_string(dataset) + "\n");
        }

        // Close resources.
        H5Dclose(dataset);
        H5Pclose(prop);
        H5Sclose(dataspace);
        H5Fclose(file);
    }
    
    /*
      Write the only component in the multifab to the dataset given by field_name.
      Uses hdf5-parallel.
    */
    void output_write_field(const std::string& file_path,
                            const std::string& field_path,
                            const MultiFab& mf, const int comp)
    {

        BL_PROFILE("output_write_field");

        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;
        int mpi_rank;
        MPI_Comm_rank(comm, &mpi_rank);

        // Create the file access prop list.
        hid_t pa_plist = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pa_plist, comm, info);

        // Open the file, and the group.
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, pa_plist);
        // Open the field dataset.
        hid_t dataset = H5Dopen(file, field_path.c_str(), H5P_DEFAULT);

        // Make sure the dataset is there.
        if (dataset < 0)
        {
            amrex::Abort("Error on rank " + std::to_string(mpi_rank) +
                         ". Count not find dataset " + field_path + "\n");
        }
        
        // Grab the dataspace of the field dataset from file.
        hid_t file_dataspace = H5Dget_space(dataset);

        // Create collective io prop list.
        hid_t collective_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(collective_plist, H5FD_MPIO_INDEPENDENT);

        // Iterate over Fabs, select matching hyperslab and write.
        hid_t status;
        // slab lo index and shape.
#if (AMREX_SPACEDIM == 3)
        hsize_t slab_offsets[3], slab_dims[3];
#else
        hsize_t slab_offsets[2], slab_dims[2];
#endif
        hid_t slab_dataspace;

        int write_count = 0;

        std::vector<Real> transposed_data;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            const int *lo_vec = box.loVect();
            const int *hi_vec = box.hiVect();
            
            transposed_data.resize(box.numPts(), 0.0);
            
            // Set slab offset and shape.
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                AMREX_ASSERT(lo_vec[idim] >= 0);
                AMREX_ASSERT(hi_vec[idim] > lo_vec[idim]);
                slab_offsets[idim] = lo_vec[idim];
                slab_dims[idim] = hi_vec[idim] - lo_vec[idim] + 1;
            }
            
            int cnt = 0;
            AMREX_D_TERM(
                         for (int i = lo_vec[0]; i <= hi_vec[0]; ++i),
                         for (int j = lo_vec[1]; j <= hi_vec[1]; ++j),
                         for (int k = lo_vec[2]; k <= hi_vec[2]; ++k))
                transposed_data[cnt++] = mf[mfi](IntVect(AMREX_D_DECL(i, j, k)), comp);

            // Create the slab space.
            slab_dataspace = H5Screate_simple(AMREX_SPACEDIM, slab_dims, NULL);

            // Select the hyperslab matching this fab.
            status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET,
                                         slab_offsets, NULL, slab_dims, NULL);
            if (status < 0)
            {
                amrex::Abort("Error on rank " + std::to_string(mpi_rank) +
                             " could not select hyperslab.\n");
            }

            // Write this pencil.
            status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, slab_dataspace,
                              file_dataspace, collective_plist, transposed_data.data());
            if (status < 0)
            {
                amrex::Abort("Error on rank " + std::to_string(mpi_rank) +
                             " could not write hyperslab.\n");
            }

            H5Sclose(slab_dataspace);
            write_count++;
        }

        ParallelDescriptor::Barrier();
        
        // Close HDF5 resources.
        H5Pclose(collective_plist);
        H5Sclose(file_dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        H5Pclose(pa_plist);
    }
}
#endif

BoostedFrameDiagnostic::
BoostedFrameDiagnostic(Real zmin_lab, Real zmax_lab, Real v_window_lab,
                       Real dt_snapshots_lab, int N_snapshots, 
                       Real gamma_boost, Real t_boost, Real dt_boost, 
                       int boost_direction, const Geometry& geom)
    : gamma_boost_(gamma_boost),
      dt_snapshots_lab_(dt_snapshots_lab),
      dt_boost_(dt_boost),
      N_snapshots_(N_snapshots),
      boost_direction_(boost_direction)
{

    BL_PROFILE("BoostedFrameDiagnostic::BoostedFrameDiagnostic");

    AMREX_ALWAYS_ASSERT(WarpX::do_boosted_frame_fields or
                        WarpX::do_boosted_frame_particles);
    
    inv_gamma_boost_ = 1.0 / gamma_boost_;
    beta_boost_ = std::sqrt(1.0 - inv_gamma_boost_*inv_gamma_boost_);
    inv_beta_boost_ = 1.0 / beta_boost_;
    
    dz_lab_ = PhysConst::c * dt_boost_ * inv_beta_boost_ * inv_gamma_boost_;
    inv_dz_lab_ = 1.0 / dz_lab_;
    Nx_lab_ = geom.Domain().length(0);
#if (AMREX_SPACEDIM == 3)
    Ny_lab_ = geom.Domain().length(1);
#else
    Ny_lab_ = 1;
#endif
    Nz_lab_ = static_cast<unsigned>((zmax_lab - zmin_lab) * inv_dz_lab_);
        
    writeMetaData();

    if (WarpX::do_boosted_frame_fields) data_buffer_.resize(N_snapshots);
    if (WarpX::do_boosted_frame_particles) particles_buffer_.resize(N_snapshots);
    for (int i = 0; i < N_snapshots; ++i) {
        Real t_lab = i * dt_snapshots_lab_;
        LabSnapShot snapshot(t_lab, t_boost,
                             zmin_lab + v_window_lab * t_lab,
                             zmax_lab + v_window_lab * t_lab, i, *this);
        snapshots_.push_back(snapshot);
        buff_counter_.push_back(0);
        if (WarpX::do_boosted_frame_fields) data_buffer_[i].reset( nullptr );
    }

    AMREX_ALWAYS_ASSERT(max_box_size_ >= num_buffer_);
}

void BoostedFrameDiagnostic::Flush(const Geometry& geom)
{
    BL_PROFILE("BoostedFrameDiagnostic::Flush");
    
    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    auto & mypc = WarpX::GetInstance().GetPartContainer();
    const std::vector<std::string> species_names = mypc.GetSpeciesNames();
    
    for (int i = 0; i < N_snapshots_; ++i) {

        int i_lab = (snapshots_[i].current_z_lab - snapshots_[i].zmin_lab) / dz_lab_;
        
        if (buff_counter_[i] != 0) {
            if (WarpX::do_boosted_frame_fields) {
                const BoxArray& ba = data_buffer_[i]->boxArray();
                const int hi = ba[0].bigEnd(boost_direction_);
                const int lo = hi - buff_counter_[i] + 1;
                
                Box buff_box = geom.Domain();
                buff_box.setSmall(boost_direction_, lo);
                buff_box.setBig(boost_direction_, hi);
                
                BoxArray buff_ba(buff_box);
                buff_ba.maxSize(max_box_size_);
                DistributionMapping buff_dm(buff_ba);
                
                const int ncomp = data_buffer_[i]->nComp();
                
                MultiFab tmp(buff_ba, buff_dm, ncomp, 0);
                
                tmp.copy(*data_buffer_[i], 0, 0, ncomp);

#ifdef WARPX_USE_HDF5
                for (int comp = 0; comp < ncomp; ++comp)
                    output_write_field(snapshots_[i].file_name, mesh_field_names[comp], tmp, comp);
#else                
                std::stringstream ss;
                ss << snapshots_[i].file_name << "/Level_0/" << Concatenate("buffer", i_lab, 5);
                VisMF::Write(tmp, ss.str());
#endif
            }
            
            if (WarpX::do_boosted_frame_particles) {
                for (int j = 0; j < mypc.nSpecies(); ++j) {
#ifdef WARPX_USE_HDF5
                    writeParticleDataHDF5(particles_buffer_[i][j],
                                          snapshots_[i].file_name,
                                          species_names[j]);
#else
                    std::stringstream part_ss;
                    part_ss << snapshots_[i].file_name + "/" + species_names[j] + "/";
                    writeParticleData(particles_buffer_[i][j], part_ss.str(), i_lab);
#endif
                }
                particles_buffer_[i].clear();
            }
            buff_counter_[i] = 0;
        }
    }

    VisMF::SetHeaderVersion(current_version);
}

void
BoostedFrameDiagnostic::
writeLabFrameData(const MultiFab* cell_centered_data,
                  const MultiParticleContainer& mypc,
                  const Geometry& geom, const Real t_boost, const Real dt) {

    BL_PROFILE("BoostedFrameDiagnostic::writeLabFrameData");

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    const RealBox& domain_z_boost = geom.ProbDomain();
    const Real zlo_boost = domain_z_boost.lo(boost_direction_);
    const Real zhi_boost = domain_z_boost.hi(boost_direction_);

    const std::vector<std::string> species_names = mypc.GetSpeciesNames();
    
    for (int i = 0; i < N_snapshots_; ++i) {
        const Real old_z_boost = snapshots_[i].current_z_boost;
        snapshots_[i].updateCurrentZPositions(t_boost,
                                              inv_gamma_boost_,
                                              inv_beta_boost_);
        
        if ( (snapshots_[i].current_z_boost < zlo_boost) or
             (snapshots_[i].current_z_boost > zhi_boost) or
             (snapshots_[i].current_z_lab < snapshots_[i].zmin_lab) or
             (snapshots_[i].current_z_lab > snapshots_[i].zmax_lab) ) continue;

        int i_lab = (snapshots_[i].current_z_lab - snapshots_[i].zmin_lab) / dz_lab_;

        if (buff_counter_[i] == 0) {
            if (WarpX::do_boosted_frame_fields) {
                const int ncomp = cell_centered_data->nComp();
                Box buff_box = geom.Domain();
                buff_box.setSmall(boost_direction_, i_lab - num_buffer_ + 1);
                buff_box.setBig(boost_direction_, i_lab);
                BoxArray buff_ba(buff_box);
                buff_ba.maxSize(max_box_size_);
                DistributionMapping buff_dm(buff_ba);
                data_buffer_[i].reset( new MultiFab(buff_ba, buff_dm, ncomp, 0) );
            }
            if (WarpX::do_boosted_frame_particles) particles_buffer_[i].resize(mypc.nSpecies());
        }

        if (WarpX::do_boosted_frame_fields) {
            const int ncomp = cell_centered_data->nComp();
            const int start_comp = 0;
            const bool interpolate = true;
            std::unique_ptr<MultiFab> slice = amrex::get_slice_data(boost_direction_,
                                                                    snapshots_[i].current_z_boost,
                                                                    *cell_centered_data, geom,
                                                                    start_comp, ncomp, interpolate);
            
            // transform it to the lab frame
            for (MFIter mfi(*slice); mfi.isValid(); ++mfi) {
                const Box& tile_box = mfi.tilebox();
                WRPX_LORENTZ_TRANSFORM_Z(BL_TO_FORTRAN_ANYD((*slice)[mfi]),
                                         BL_TO_FORTRAN_BOX(tile_box),
                                         &gamma_boost_, &beta_boost_);
            }
            
            Real dx = geom.CellSize(boost_direction_);
            int i_boost = (snapshots_[i].current_z_boost - geom.ProbLo(boost_direction_))/dx;
            Box slice_box = geom.Domain();
            slice_box.setSmall(boost_direction_, i_boost);
            slice_box.setBig(boost_direction_, i_boost);
            BoxArray slice_ba(slice_box);
            slice_ba.maxSize(max_box_size_);
            MultiFab tmp(slice_ba, data_buffer_[i]->DistributionMap(), ncomp, 0);
            tmp.copy(*slice, 0, 0, ncomp);
            
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(tmp, true); mfi.isValid(); ++mfi) {
                const Box& tile_box  = mfi.tilebox();
                WRPX_COPY_SLICE(BL_TO_FORTRAN_BOX(tile_box),
                                BL_TO_FORTRAN_ANYD(tmp[mfi]),
                                BL_TO_FORTRAN_ANYD((*data_buffer_[i])[mfi]),
                                &ncomp, &i_boost, &i_lab);
            }
        }
        
        if (WarpX::do_boosted_frame_particles) {
            mypc.GetLabFrameData(snapshots_[i].file_name, i_lab, boost_direction_,
                                 old_z_boost, snapshots_[i].current_z_boost,
                                 t_boost, snapshots_[i].t_lab, dt, particles_buffer_[i]);
        }


        ++buff_counter_[i];
        
        if (buff_counter_[i] == num_buffer_) {

            if (WarpX::do_boosted_frame_fields) {
#ifdef WARPX_USE_HDF5
                for (int comp = 0; comp < data_buffer_[i]->nComp(); ++comp)
                    output_write_field(snapshots_[i].file_name, mesh_field_names[comp],
                                       *data_buffer_[i], comp);
#else
                std::stringstream mesh_ss;
                mesh_ss << snapshots_[i].file_name << "/Level_0/" << Concatenate("buffer", i_lab, 5);
                VisMF::Write(*data_buffer_[i], mesh_ss.str());
#endif
            }
            
            if (WarpX::do_boosted_frame_particles) {
                for (int j = 0; j < mypc.nSpecies(); ++j) {
#ifdef WARPX_USE_HDF5
                    writeParticleDataHDF5(particles_buffer_[i][j],
                                          snapshots_[i].file_name,
                                          species_names[j]);
#else
                    std::stringstream part_ss;
                    part_ss << snapshots_[i].file_name + "/" + species_names[j] + "/";
                    writeParticleData(particles_buffer_[i][j], part_ss.str(), i_lab);
#endif
                }            
                particles_buffer_[i].clear();
            }
            buff_counter_[i] = 0;
        }
    }
        
    VisMF::SetHeaderVersion(current_version);    
}

#ifdef WARPX_USE_HDF5
void
BoostedFrameDiagnostic::
writeParticleDataHDF5(const WarpXParticleContainer::DiagnosticParticleData& pdata,
                      const std::string& name, const std::string& species_name)
{
    auto np = pdata.GetRealData(DiagIdx::w).size();
    
    Vector<long> particle_counts(ParallelDescriptor::NProcs(), 0);
    Vector<long> particle_offsets(ParallelDescriptor::NProcs(), 0);
    
    ParallelAllGather::AllGather(np, particle_counts.data(), ParallelContext::CommunicatorAll());
    
    long total_np = 0;
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
        particle_offsets[i] = total_np;
        total_np += particle_counts[i];
    }

    if (total_np == 0) return;
    
    long old_np = 0;
    if (ParallelDescriptor::IOProcessor())
    {
        for (int k = 0; k < static_cast<int>(particle_field_names.size()); ++k)
        {
            std::string field_path = species_name + "/" + particle_field_names[k];
            old_np = output_resize_particle_field(name, field_path, total_np);
        }
    }

    ParallelDescriptor::ReduceLongMax(old_np);
    
    // Write data here
    for (int k = 0; k < static_cast<int>(particle_field_names.size()); ++k)
    {
        std::string field_path = species_name + "/" + particle_field_names[k];
        output_write_particle_field(name, field_path,
                                    pdata.GetRealData(k).data(),
                                    particle_counts[ParallelDescriptor::MyProc()],
                                    particle_offsets[ParallelDescriptor::MyProc()] + old_np);
    }    
}
#endif

void
BoostedFrameDiagnostic::
writeParticleData(const WarpXParticleContainer::DiagnosticParticleData& pdata,
                  const std::string& name, const int i_lab)
{
    BL_PROFILE("BoostedFrameDiagnostic::writeParticleData");
    
    std::string field_name;
    std::ofstream ofs;

    const int MyProc = ParallelDescriptor::MyProc();
    auto np = pdata.GetRealData(DiagIdx::w).size();

    if (np == 0) return;

    field_name = name + Concatenate("w_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::w).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("x_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::x).data(), np, ofs);
    ofs.close();    

    field_name = name + Concatenate("y_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::y).data(), np, ofs);
    ofs.close();    

    field_name = name + Concatenate("z_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::z).data(), np, ofs);
    ofs.close();    
    
    field_name = name + Concatenate("ux_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::ux).data(), np, ofs);
    ofs.close();    

    field_name = name + Concatenate("uy_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::uy).data(), np, ofs);
    ofs.close();    

    field_name = name + Concatenate("uz_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeRealData(pdata.GetRealData(DiagIdx::uz).data(), np, ofs);
    ofs.close();
}

void
BoostedFrameDiagnostic::
writeMetaData () 
{
    BL_PROFILE("BoostedFrameDiagnostic::writeMetaData");

    if (ParallelDescriptor::IOProcessor()) {
        std::string DiagnosticDirectory = "lab_frame_data";
        
        if (!UtilCreateDirectory(DiagnosticDirectory, 0755))
            CreateDirectoryFailed(DiagnosticDirectory);

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(DiagnosticDirectory + "/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if(!HeaderFile.good())
            FileOpenFailed(HeaderFileName);
        
        HeaderFile.precision(17);
        
        HeaderFile << N_snapshots_ << "\n";
        HeaderFile << dt_snapshots_lab_ << "\n";    
        HeaderFile << gamma_boost_ << "\n";
        HeaderFile << beta_boost_ << "\n";
        HeaderFile << dz_lab_ << "\n";
        HeaderFile << Nz_lab_ << "\n";
    }
}

BoostedFrameDiagnostic::LabSnapShot::
LabSnapShot(Real t_lab_in, Real t_boost, Real zmin_lab_in, 
            Real zmax_lab_in, int file_num_in, const BoostedFrameDiagnostic& bfd)
    : t_lab(t_lab_in),
      zmin_lab(zmin_lab_in),
      zmax_lab(zmax_lab_in),
      file_num(file_num_in),
      my_bfd(bfd)
{
    current_z_lab = 0.0;
    current_z_boost = 0.0;
    updateCurrentZPositions(t_boost, my_bfd.inv_gamma_boost_, my_bfd.inv_beta_boost_);
    initial_i = (current_z_lab - zmin_lab) / my_bfd.dz_lab_;
    file_name = Concatenate("lab_frame_data/snapshot", file_num, 5);

#ifdef WARPX_USE_HDF5
    if (ParallelDescriptor::IOProcessor())
    {
        output_create(file_name);
    }

    ParallelDescriptor::Barrier();
    
    if (ParallelDescriptor::IOProcessor())
    {
        if (WarpX::do_boosted_frame_fields)
        {
            for (int comp = 0; comp < static_cast<int>(mesh_field_names.size()); ++comp) {
                output_create_field(file_name, mesh_field_names[comp],
                                    my_bfd.Nx_lab_,
                                    my_bfd.Ny_lab_,
                                    my_bfd.Nz_lab_+1);
            }
        }
    }

    ParallelDescriptor::Barrier();
    
    if (WarpX::do_boosted_frame_particles)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::vector<std::string> species_names = mypc.GetSpeciesNames();
        for (int j = 0; j < mypc.nSpecies(); ++j)
        {
            output_create_species_group(file_name, species_names[j]);
            for (int k = 0; k < static_cast<int>(particle_field_names.size()); ++k)
            {
                std::string field_path = species_names[j] + "/" + particle_field_names[k];
                output_create_particle_field(file_name, field_path);
            }
        }
    }    
#else    
    if (ParallelDescriptor::IOProcessor()) {
        
        if (!UtilCreateDirectory(file_name, 0755))
            CreateDirectoryFailed(file_name);

        const int nlevels = 1;
        for(int i = 0; i < nlevels; ++i) {
            const std::string &fullpath = LevelFullPath(i, file_name);
            if (!UtilCreateDirectory(fullpath, 0755))
                CreateDirectoryFailed(fullpath);
        }

        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::vector<std::string> species_names = mypc.GetSpeciesNames();        
        int nspecies = mypc.nSpecies();
        
        const std::string particles_prefix = "particle";
        for(int i = 0; i < nspecies; ++i) {
            const std::string fullpath = file_name + "/" + species_names[i];
            if (!UtilCreateDirectory(fullpath, 0755))
                CreateDirectoryFailed(fullpath);
        }
    }
#endif
    ParallelDescriptor::Barrier();

    writeSnapShotHeader();
}

void
BoostedFrameDiagnostic::LabSnapShot::
updateCurrentZPositions(Real t_boost, Real inv_gamma, Real inv_beta) {
    current_z_boost = (t_lab*inv_gamma - t_boost)*PhysConst::c*inv_beta;
    current_z_lab =   (t_lab - t_boost*inv_gamma)*PhysConst::c*inv_beta;
}

void
BoostedFrameDiagnostic::LabSnapShot::
writeSnapShotHeader() {
#ifndef WARPX_USE_HDF5
    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(file_name + "/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if(!HeaderFile.good())
            FileOpenFailed(HeaderFileName);
        
        HeaderFile.precision(17);
        
        HeaderFile << t_lab << "\n";
        HeaderFile << zmin_lab << "\n";
        HeaderFile << zmax_lab << "\n";
    }
#endif
}
