/* Copyright 2019 Andrew Myers, Axel Huebl, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_C.H>

#include "BackTransformedDiagnostic.H"
#include "SliceDiagnostic.H"
#include "WarpX.H"

using namespace amrex;

#ifdef WARPX_USE_HDF5

#include <hdf5.h>

/*
  Helper functions for doing the HDF5 IO.

 */
namespace
{

    const std::vector<std::string> particle_field_names = {"w", "x", "y",
                                                           "z", "ux", "uy", "uz"};

    /*
      Creates the HDF5 file in truncate mode and closes it.
      Should be run only by the root process.
    */
    void output_create(const std::string& file_path) {
        WARPX_PROFILE("output_create");
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
        WARPX_PROFILE("output_write_metadata");
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
        WARPX_PROFILE("output_create_field");

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
      Resize an extendible dataset, suitable for storing particle data.
      Should be run only by the master rank.
    */
    long output_resize_particle_field(const std::string& file_path, const std::string& field_path,
                                      const long num_to_add)
    {
        WARPX_PROFILE("output_resize_particle_field");

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
      Writes to a dataset that has been extended to the proper size. Suitable for writing particle data.
      Should be run on all ranks collectively.
    */
    void output_write_particle_field(const std::string& file_path, const std::string& field_path,
                                     const Real* data_ptr, const long count, const long index)
    {
        WARPX_PROFILE("output_write_particle_field");

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
        WARPX_PROFILE("output_create_particle_field");

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
                            const MultiFab& mf, const int comp,
                            const int lo_x, const int lo_y, const int lo_z)
    {

        WARPX_PROFILE("output_write_field");

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
        int shift[3];
        shift[0] = lo_x;
        shift[1] = lo_y;
        shift[2] = lo_z;
#else
        hsize_t slab_offsets[2], slab_dims[2];
        int shift[2];
        shift[0] = lo_x;
        shift[1] = lo_z;
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
                slab_offsets[idim] = lo_vec[idim] - shift[idim];
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

bool compare_tlab_uptr(const std::unique_ptr<LabFrameDiag>&a,
                       const std::unique_ptr<LabFrameDiag>&b)
{
    return a->m_t_lab < b->m_t_lab;
}

namespace
{
void
LorentzTransformZ(MultiFab& data, Real gamma_boost, Real beta_boost)
{
    // Loop over tiles/boxes and in-place convert each slice from boosted
    // frame to back-transformed lab frame.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(data, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tile_box = mfi.tilebox();
        Array4< Real > arr = data[mfi].array();
        // arr(x,y,z,comp) where 0->9 comps are
        // Ex Ey Ez Bx By Bz jx jy jz rho
        Real clight = PhysConst::c;
        ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Transform the transverse E and B fields. Note that ez and bz are not
                // changed by the tranform.
                Real e_lab, b_lab, j_lab, r_lab;
                e_lab = gamma_boost * (arr(i, j, k, 0) +
                                       beta_boost*clight*arr(i, j, k, 4));
                b_lab = gamma_boost * (arr(i, j, k, 4) +
                                       beta_boost*arr(i, j, k, 0)/clight);

                arr(i, j, k, 0) = e_lab;
                arr(i, j, k, 4) = b_lab;

                e_lab = gamma_boost * (arr(i, j, k, 1) -
                                       beta_boost*clight*arr(i, j, k, 3));
                b_lab = gamma_boost * (arr(i, j, k, 3) -
                                       beta_boost*arr(i, j, k, 1)/clight);

                arr(i, j, k, 1) = e_lab;
                arr(i, j, k, 3) = b_lab;

                // Transform the charge and current density. Only the z component of j is affected.
                j_lab = gamma_boost*(arr(i, j, k, 8) +
                                     beta_boost*clight*arr(i, j, k, 9));
                r_lab = gamma_boost*(arr(i, j, k, 9) +
                                     beta_boost*arr(i, j, k, 8)/clight);

                arr(i, j, k, 8) = j_lab;
                arr(i, j, k, 9) = r_lab;
            }
        );
    }
}
}

BackTransformedDiagnostic::
BackTransformedDiagnostic(Real zmin_lab, Real zmax_lab, Real v_window_lab,
                       Real dt_snapshots_lab, int N_snapshots,
                       Real dt_slice_snapshots_lab, int N_slice_snapshots,
                       Real gamma_boost, Real t_boost, Real dt_boost,
                       int boost_direction, const Geometry& geom,
                       amrex::RealBox& slice_realbox,
                       amrex::Real particle_slice_width_lab)
    : m_gamma_boost_(gamma_boost),
      m_dt_snapshots_lab_(dt_snapshots_lab),
      m_dt_boost_(dt_boost),
      m_N_snapshots_(N_snapshots),
      m_boost_direction_(boost_direction),
      m_N_slice_snapshots_(N_slice_snapshots),
      m_dt_slice_snapshots_lab_(dt_slice_snapshots_lab),
      m_particle_slice_width_lab_(particle_slice_width_lab)
{


    WARPX_PROFILE("BackTransformedDiagnostic::BackTransformedDiagnostic");

    AMREX_ALWAYS_ASSERT(WarpX::do_back_transformed_fields or
                        WarpX::do_back_transformed_particles);

    m_inv_gamma_boost_ = 1.0 / m_gamma_boost_;
    m_beta_boost_ = std::sqrt(1.0 - m_inv_gamma_boost_*m_inv_gamma_boost_);
    m_inv_beta_boost_ = 1.0 / m_beta_boost_;

    m_dz_lab_ = PhysConst::c * m_dt_boost_ * m_inv_beta_boost_ * m_inv_gamma_boost_;
    m_inv_dz_lab_ = 1.0 / m_dz_lab_;
    int Nz_lab = static_cast<unsigned>((zmax_lab - zmin_lab) * m_inv_dz_lab_);
    int Nx_lab = geom.Domain().length(0);
#if (AMREX_SPACEDIM == 3)
    int Ny_lab = geom.Domain().length(1);
    IntVect prob_ncells_lab = {Nx_lab, Ny_lab, Nz_lab};
#else
    // Ny_lab = 1;
    IntVect prob_ncells_lab = {Nx_lab, Nz_lab};
#endif
    writeMetaData();

    // Query fields to dump
    std::vector<std::string> user_fields_to_dump;
    ParmParse pp("warpx");
    bool do_user_fields;
    do_user_fields = pp.queryarr("back_transformed_diag_fields",
                                 user_fields_to_dump);
    // If user specifies fields to dump, overwrite ncomp_to_dump,
    // map_actual_fields_to_dump and mesh_field_names.
    for (int i = 0; i < 10; ++i)
        map_actual_fields_to_dump.push_back(i);

    if (do_user_fields){
        m_ncomp_to_dump = user_fields_to_dump.size();
        map_actual_fields_to_dump.resize(m_ncomp_to_dump);
        m_mesh_field_names.resize(m_ncomp_to_dump);
        for (int i=0; i<m_ncomp_to_dump; i++){
            std::string fieldstr = user_fields_to_dump[i];
            m_mesh_field_names[i] = fieldstr;
            map_actual_fields_to_dump[i] = m_possible_fields_to_dump[fieldstr];
        }
    }

    // allocating array with total number of lab frame diags (snapshots+slices)
    m_LabFrameDiags_.resize(N_snapshots+N_slice_snapshots);

    for (int i = 0; i < N_snapshots; ++i) {
        Real t_lab = i * m_dt_snapshots_lab_;
        // Get simulation domain physical coordinates (in boosted frame).
        RealBox prob_domain_lab = geom.ProbDomain();
        // Replace z bounds by lab-frame coordinates
        // x and y bounds are the same for back-transformed lab frame and boosted frame
        prob_domain_lab.setLo(AMREX_SPACEDIM-1, zmin_lab + v_window_lab * t_lab);
        prob_domain_lab.setHi(AMREX_SPACEDIM-1, zmax_lab + v_window_lab * t_lab);
        Box diag_box = geom.Domain();
        m_LabFrameDiags_[i].reset(new LabFrameSnapShot(t_lab, t_boost,
                                m_inv_gamma_boost_, m_inv_beta_boost_, m_dz_lab_,
                                prob_domain_lab, prob_ncells_lab,
                                m_ncomp_to_dump, m_mesh_field_names, prob_domain_lab,
                                diag_box, i));
    }


    for (int i = 0; i < N_slice_snapshots; ++i) {

        IntVect slice_ncells_lab ;

        // To construct LabFrameSlice(), the location of lo() and hi() of the
        // reduced diag is computed using the user-defined values of the
        // reduced diag (1D, 2D, or 3D). For visualization of the diagnostics,
        // the number of cells in each dimension is required and
        // is computed below for the reduced back-transformed lab-frame diag,
        // similar to the full-diag.
        const amrex::Real* current_slice_lo = slice_realbox.lo();
        const amrex::Real* current_slice_hi = slice_realbox.hi();

        const amrex::Real zmin_slice_lab = current_slice_lo[AMREX_SPACEDIM-1] /
                                          ( (1.+m_beta_boost_)*m_gamma_boost_);
        const amrex::Real zmax_slice_lab = current_slice_hi[AMREX_SPACEDIM-1] /
                                          ( (1.+m_beta_boost_)*m_gamma_boost_);
        int Nz_slice_lab = static_cast<unsigned>((
                           zmax_slice_lab - zmin_slice_lab) * m_inv_dz_lab_);
        int Nx_slice_lab = ( current_slice_hi[0] - current_slice_lo[0] ) /
                           geom.CellSize(0);
        if (Nx_slice_lab == 0 ) Nx_slice_lab = 1;
        // if the x-dimension is reduced, increase total_cells by 1
        // to be consistent with the number of cells created for the output.
        if (Nx_lab != Nx_slice_lab) Nx_slice_lab++;
#if (AMREX_SPACEDIM == 3)
        int Ny_slice_lab = ( current_slice_hi[1] - current_slice_lo[1]) /
                             geom.CellSize(1);
        if (Ny_slice_lab == 0 ) Ny_slice_lab = 1;
        // if the y-dimension is reduced, increase total_cells by 1
        // to be consistent with the number of cells created for the output.
        if (Ny_lab != Ny_slice_lab) Ny_slice_lab++;
        slice_ncells_lab = {Nx_slice_lab, Ny_slice_lab, Nz_slice_lab};
#else
        slice_ncells_lab = {Nx_slice_lab, Nz_slice_lab};
#endif

        IntVect slice_lo(AMREX_D_DECL(0,0,0));
        IntVect slice_hi(AMREX_D_DECL(1,1,1));

        for ( int i_dim=0; i_dim<AMREX_SPACEDIM; ++i_dim)
        {
           slice_lo[i_dim] = (slice_realbox.lo(i_dim) - geom.ProbLo(i_dim) -
                              0.5*geom.CellSize(i_dim))/geom.CellSize(i_dim);
           slice_hi[i_dim] = (slice_realbox.hi(i_dim) - geom.ProbLo(i_dim) -
                              0.5*geom.CellSize(i_dim))/geom.CellSize(i_dim);
           if (slice_lo[i_dim] == slice_hi[i_dim])
           {
              slice_hi[i_dim] = slice_lo[i_dim] + 1;
           }
        }
        Box stmp(slice_lo,slice_hi);
        Box slicediag_box = stmp;

        Real t_slice_lab = i * m_dt_slice_snapshots_lab_ ;
        RealBox prob_domain_lab = geom.ProbDomain();
        // replace z bounds by lab-frame coordinates
        prob_domain_lab.setLo(AMREX_SPACEDIM-1, zmin_lab + v_window_lab * t_slice_lab);
        prob_domain_lab.setHi(AMREX_SPACEDIM-1, zmax_lab + v_window_lab * t_slice_lab);
        RealBox slice_dom_lab = slice_realbox;
        // replace z bounds of slice in lab-frame coordinates
        // note : x and y bounds are the same for lab and boosted frames
        // initial lab slice extent //
        slice_dom_lab.setLo(AMREX_SPACEDIM-1, zmin_slice_lab + v_window_lab * t_slice_lab );
        slice_dom_lab.setHi(AMREX_SPACEDIM-1, zmax_slice_lab +
                                         v_window_lab * t_slice_lab );

        // construct labframeslice
        m_LabFrameDiags_[i+N_snapshots].reset(new LabFrameSlice(t_slice_lab, t_boost,
                                m_inv_gamma_boost_, m_inv_beta_boost_, m_dz_lab_,
                                prob_domain_lab, slice_ncells_lab,
                                m_ncomp_to_dump, m_mesh_field_names, slice_dom_lab,
                                slicediag_box, i, m_particle_slice_width_lab_));
    }
    // sort diags based on their respective t_lab
    std::stable_sort(m_LabFrameDiags_.begin(), m_LabFrameDiags_.end(), compare_tlab_uptr);

    AMREX_ALWAYS_ASSERT(m_max_box_size_ >= m_num_buffer_);
}

void BackTransformedDiagnostic::Flush(const Geometry& /*geom*/)
{
    WARPX_PROFILE("BackTransformedDiagnostic::Flush");

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    auto & mypc = WarpX::GetInstance().GetPartContainer();
    const std::vector<std::string> species_names = mypc.GetSpeciesNames();

    // Loop over BFD snapshots
    for (int i = 0; i < m_LabFrameDiags_.size(); ++i) {

        Real zmin_lab = m_LabFrameDiags_[i]->m_prob_domain_lab_.lo(AMREX_SPACEDIM-1);
        int i_lab = (m_LabFrameDiags_[i]->m_current_z_lab - zmin_lab) / m_dz_lab_;

        if (m_LabFrameDiags_[i]->m_buff_counter_ != 0) {
            if (WarpX::do_back_transformed_fields) {
                const BoxArray& ba = m_LabFrameDiags_[i]->m_data_buffer_->boxArray();
                const int hi = ba[0].bigEnd(m_boost_direction_);
                const int lo = hi - m_LabFrameDiags_[i]->m_buff_counter_ + 1;

                //Box buff_box = geom.Domain();
                Box buff_box = m_LabFrameDiags_[i]->m_buff_box_;
                buff_box.setSmall(m_boost_direction_, lo);
                buff_box.setBig(m_boost_direction_, hi);

                BoxArray buff_ba(buff_box);
                buff_ba.maxSize(m_max_box_size_);
                DistributionMapping buff_dm(buff_ba);

                const int ncomp = m_LabFrameDiags_[i]->m_data_buffer_->nComp();

                MultiFab tmp(buff_ba, buff_dm, ncomp, 0);

                tmp.copy(*m_LabFrameDiags_[i]->m_data_buffer_, 0, 0, ncomp);

#ifdef WARPX_USE_HDF5
                for (int comp = 0; comp < ncomp; ++comp) {
                    output_write_field(m_LabFrameDiags_[i]->m_file_name,
                                       m_mesh_field_names[comp], tmp, comp,
                                       lbound(buff_box).x, lbound(buff_box).y,
                                       lbound(buff_box).z);
                }
#else
                std::stringstream ss;
                ss << m_LabFrameDiags_[i]->m_file_name << "/Level_0/"
                   << Concatenate("buffer", i_lab, 5);
                VisMF::Write(tmp, ss.str());
#endif
            }

            if (WarpX::do_back_transformed_particles) {
                // Loop over species to be dumped to BFD
                for (int j = 0; j < mypc.nSpeciesBackTransformedDiagnostics(); ++j) {
                    // Get species name
                    std::string species_name =
                        species_names[mypc.mapSpeciesBackTransformedDiagnostics(j)];
#ifdef WARPX_USE_HDF5
                    // Dump species data
                    writeParticleDataHDF5(m_LabFrameDiags_[i]->m_particles_buffer_[j],
                                          m_LabFrameDiags_[i]->m_file_name,
                                          species_name);
#else
                    std::stringstream part_ss;
                    part_ss << m_LabFrameDiags_[i]->m_file_name + "/" +
                               species_name + "/";
                    // Dump species data
                    writeParticleData(m_LabFrameDiags_[i]->m_particles_buffer_[j],
                                      part_ss.str(), i_lab);
#endif
                }
                m_LabFrameDiags_[i]->m_particles_buffer_.clear();
            }
            m_LabFrameDiags_[i]->m_buff_counter_ = 0;
        }
    }

    VisMF::SetHeaderVersion(current_version);
}





void
BackTransformedDiagnostic::
writeLabFrameData(const MultiFab* cell_centered_data,
                  const MultiParticleContainer& mypc,
                  const Geometry& geom, const Real t_boost, const Real dt) {

    WARPX_PROFILE("BackTransformedDiagnostic::writeLabFrameData");
    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    const RealBox& domain_z_boost = geom.ProbDomain();
    const Real zlo_boost = domain_z_boost.lo(m_boost_direction_);
    const Real zhi_boost = domain_z_boost.hi(m_boost_direction_);

    const std::vector<std::string> species_names = mypc.GetSpeciesNames();
    Real prev_t_lab = -dt;
    std::unique_ptr<amrex::MultiFab> tmp_slice_ptr;
    std::unique_ptr<amrex::MultiFab> slice;
    amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer;

    // Loop over snapshots
    for (int i = 0; i < m_LabFrameDiags_.size(); ++i) {
        // Get updated z position of snapshot
        const Real old_z_boost = m_LabFrameDiags_[i]->m_current_z_boost;
        m_LabFrameDiags_[i]->updateCurrentZPositions(t_boost,
                                              m_inv_gamma_boost_,
                                              m_inv_beta_boost_);

        Real diag_zmin_lab = m_LabFrameDiags_[i]->m_diag_domain_lab_.lo(AMREX_SPACEDIM-1);
        Real diag_zmax_lab = m_LabFrameDiags_[i]->m_diag_domain_lab_.hi(AMREX_SPACEDIM-1);

        if ( ( m_LabFrameDiags_[i]->m_current_z_boost < zlo_boost) or
             ( m_LabFrameDiags_[i]->m_current_z_boost > zhi_boost) or
             ( m_LabFrameDiags_[i]->m_current_z_lab < diag_zmin_lab) or
             ( m_LabFrameDiags_[i]->m_current_z_lab > diag_zmax_lab) ) continue;

        // Get z index of data_buffer_ (i.e. in the lab frame) where
        // simulation domain (t', [zmin',zmax']), back-transformed to lab
        // frame, intersects with snapshot.
        Real dom_zmin_lab = m_LabFrameDiags_[i]->m_prob_domain_lab_.lo(AMREX_SPACEDIM-1);
        int i_lab = ( m_LabFrameDiags_[i]->m_current_z_lab - dom_zmin_lab) / m_dz_lab_;
        // If buffer of snapshot i is empty...
        if ( m_LabFrameDiags_[i]->m_buff_counter_ == 0) {
            // ... reset fields buffer data_buffer_
            if (WarpX::do_back_transformed_fields) {
                m_LabFrameDiags_[i]->m_buff_box_.setSmall(m_boost_direction_,
                                             i_lab - m_num_buffer_ + 1);
                m_LabFrameDiags_[i]->m_buff_box_.setBig(m_boost_direction_, i_lab);

                BoxArray buff_ba(m_LabFrameDiags_[i]->m_buff_box_);
                buff_ba.maxSize(m_max_box_size_);
                DistributionMapping buff_dm(buff_ba);
                m_LabFrameDiags_[i]->m_data_buffer_.reset( new MultiFab(buff_ba,
                                                buff_dm, m_ncomp_to_dump, 0) );
            }
            // ... reset particle buffer particles_buffer_[i]
            if (WarpX::do_back_transformed_particles)
                m_LabFrameDiags_[i]->m_particles_buffer_.resize(
                                   mypc.nSpeciesBackTransformedDiagnostics());
        }

        if (WarpX::do_back_transformed_fields) {
            const int ncomp = cell_centered_data->nComp();
            const int start_comp = 0;
            const bool interpolate = true;
            // slice containing back-transformed data is generated only if t_lab != prev_t_lab and is re-used if multiple diags have the same z_lab,t_lab.
            if (m_LabFrameDiags_[i]->m_t_lab != prev_t_lab ) {
               if (slice)
               {
                 slice.reset(new MultiFab);
                 slice.reset(nullptr);
               }
               slice = amrex::get_slice_data(m_boost_direction_,
                                             m_LabFrameDiags_[i]->m_current_z_boost,
                                             *cell_centered_data, geom,
                                             start_comp, ncomp,
                                             interpolate);
               // Back-transform data to the lab-frame
               LorentzTransformZ(*slice, m_gamma_boost_, m_beta_boost_);
             }
             // Create a 2D box for the slice in the boosted frame
             Real dx = geom.CellSize(m_boost_direction_);
             int i_boost = ( m_LabFrameDiags_[i]->m_current_z_boost -
                             geom.ProbLo(m_boost_direction_))/dx;
             //Box slice_box = geom.Domain();
             Box slice_box = m_LabFrameDiags_[i]->m_buff_box_;
             slice_box.setSmall(m_boost_direction_, i_boost);
             slice_box.setBig(m_boost_direction_, i_boost);

             // Make it a BoxArray slice_ba
             BoxArray slice_ba(slice_box);
             slice_ba.maxSize(m_max_box_size_);
             tmp_slice_ptr = std::unique_ptr<MultiFab>(new MultiFab(slice_ba,
                             m_LabFrameDiags_[i]->m_data_buffer_->DistributionMap(),
                             ncomp, 0));

             // slice is re-used if the t_lab of a diag is equal to
             // that of the previous diag.
             // Back-transformed data is copied from slice
             // which has the dmap of the domain to
             // tmp_slice_ptr which has the dmap of the
             // data_buffer that stores the back-transformed data.
             tmp_slice_ptr->copy(*slice, 0, 0, ncomp);
             m_LabFrameDiags_[i]->AddDataToBuffer(*tmp_slice_ptr, i_lab,
                                               map_actual_fields_to_dump);
             tmp_slice_ptr.reset(new MultiFab);
             tmp_slice_ptr.reset(nullptr);
        }

        if (WarpX::do_back_transformed_particles) {

            if (m_LabFrameDiags_[i]->m_t_lab != prev_t_lab ) {
               if (tmp_particle_buffer.size()>0)
               {
                  tmp_particle_buffer.clear();
                  tmp_particle_buffer.shrink_to_fit();
               }
               tmp_particle_buffer.resize(mypc.nSpeciesBackTransformedDiagnostics());
               mypc.GetLabFrameData(m_LabFrameDiags_[i]->m_file_name, i_lab,
                                    m_boost_direction_, old_z_boost,
                                    m_LabFrameDiags_[i]->m_current_z_boost,
                                    t_boost, m_LabFrameDiags_[i]->m_t_lab, dt,
                                    tmp_particle_buffer);
            }
            m_LabFrameDiags_[i]->AddPartDataToParticleBuffer(tmp_particle_buffer,
                               mypc.nSpeciesBackTransformedDiagnostics());
        }

        ++m_LabFrameDiags_[i]->m_buff_counter_;
        prev_t_lab = m_LabFrameDiags_[i]->m_t_lab;
        // If buffer full, write to disk.
        if (m_LabFrameDiags_[i]->m_buff_counter_ == m_num_buffer_) {

            if (WarpX::do_back_transformed_fields) {
#ifdef WARPX_USE_HDF5

                Box buff_box = m_LabFrameDiags_[i]->m_buff_box_;
                for (int comp = 0; comp < m_LabFrameDiags_[i]->m_data_buffer_->nComp(); ++comp)
                    output_write_field(m_LabFrameDiags_[i]->m_file_name,
                                       m_mesh_field_names[comp],
                                       *m_LabFrameDiags_[i]->m_data_buffer_, comp,
                                       lbound(buff_box).x, lbound(buff_box).y,
                                       lbound(buff_box).z);
#else
                std::stringstream mesh_ss;
                mesh_ss << m_LabFrameDiags_[i]->m_file_name << "/Level_0/" <<
                           Concatenate("buffer", i_lab, 5);
                VisMF::Write( (*m_LabFrameDiags_[i]->m_data_buffer_), mesh_ss.str());
#endif
            }

            if (WarpX::do_back_transformed_particles) {
                // Loop over species to be dumped to BFD
                for (int j = 0; j < mypc.nSpeciesBackTransformedDiagnostics(); ++j) {
                    // Get species name
                    const std::string species_name = species_names[
                                      mypc.mapSpeciesBackTransformedDiagnostics(j)];
#ifdef WARPX_USE_HDF5
                    // Write data to disk (HDF5)
                    writeParticleDataHDF5(m_LabFrameDiags_[i]->m_particles_buffer_[j],
                                          m_LabFrameDiags_[i]->m_file_name,
                                          species_name);
#else
                    std::stringstream part_ss;

                    part_ss << m_LabFrameDiags_[i]->m_file_name + "/" +
                               species_name + "/";

                    // Write data to disk (custom)
                    writeParticleData(m_LabFrameDiags_[i]->m_particles_buffer_[j],
                                      part_ss.str(), i_lab);
#endif
                }
                m_LabFrameDiags_[i]->m_particles_buffer_.clear();
            }
            m_LabFrameDiags_[i]->m_buff_counter_ = 0;
        }
    }

    VisMF::SetHeaderVersion(current_version);
}

#ifdef WARPX_USE_HDF5
void
BackTransformedDiagnostic::
writeParticleDataHDF5(const WarpXParticleContainer::DiagnosticParticleData& pdata,
                      const std::string& name, const std::string& species_name)
{
    auto np = pdata.GetRealData(DiagIdx::w).size();

    Vector<long> particle_counts(ParallelDescriptor::NProcs(), 0);
    Vector<long> particle_offsets(ParallelDescriptor::NProcs(), 0);

    ParallelAllGather::AllGather(np, particle_counts.data(),
                                 ParallelContext::CommunicatorAll());

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

    // Note, this has the effect of an MPI Barrier between the above resize operation
    // and the below write.
    ParallelDescriptor::ReduceLongMax(old_np);

    // Write data here
    for (int k = 0; k < static_cast<int>(particle_field_names.size()); ++k)
    {
        std::string field_path = species_name + "/" + particle_field_names[k];
        output_write_particle_field(name, field_path,
                                    pdata.GetRealData(k).data(),
                                    particle_counts[ParallelDescriptor::MyProc()],
                                    particle_offsets[ParallelDescriptor::MyProc()]
                                    + old_np);
    }
}
#endif

void
BackTransformedDiagnostic::
writeParticleData(const WarpXParticleContainer::DiagnosticParticleData& pdata,
                  const std::string& name, const int i_lab)
{
    WARPX_PROFILE("BackTransformedDiagnostic::writeParticleData");

    std::string field_name;
    std::ofstream ofs;

    const int MyProc = ParallelDescriptor::MyProc();
    auto np = pdata.GetRealData(DiagIdx::w).size();
    if (np == 0) return;

    field_name = name + Concatenate("w_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::w).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("x_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::x).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("y_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::y).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("z_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::z).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("ux_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::ux).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("uy_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::uy).data(), np, ofs);
    ofs.close();

    field_name = name + Concatenate("uz_", i_lab, 5) + "_" + std::to_string(MyProc);
    ofs.open(field_name.c_str(), std::ios::out|std::ios::binary);
    writeData(pdata.GetRealData(DiagIdx::uz).data(), np, ofs);
    ofs.close();
}

void
BackTransformedDiagnostic::
writeMetaData ()
{
    WARPX_PROFILE("BackTransformedDiagnostic::writeMetaData");

    if (ParallelDescriptor::IOProcessor()) {
        const std::string fullpath = WarpX::lab_data_directory + "/snapshots";
        if (!UtilCreateDirectory(fullpath, 0755))
            CreateDirectoryFailed(fullpath);

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(WarpX::lab_data_directory + "/snapshots/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if(!HeaderFile.good())
            FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << m_N_snapshots_ << "\n";
        HeaderFile << m_dt_snapshots_lab_ << "\n";
        HeaderFile << m_gamma_boost_ << "\n";
        HeaderFile << m_beta_boost_ << "\n";

        if (m_N_slice_snapshots_ > 0) {
           const std::string fullpath_slice = WarpX::lab_data_directory + "/slices";
           if (!UtilCreateDirectory(fullpath_slice, 0755))
               CreateDirectoryFailed(fullpath_slice);

           VisMF::IO_Buffer io_buffer_slice(VisMF::IO_Buffer_Size);
           std::ofstream HeaderFile_slice;
           HeaderFile_slice.rdbuf()->pubsetbuf(io_buffer_slice.dataPtr(),
                                               io_buffer_slice.size());
           std::string HeaderFileName_slice(WarpX::lab_data_directory+
                                            "/slices/Header");
           HeaderFile_slice.open(HeaderFileName_slice.c_str(),
                                                std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);

           if (!HeaderFile_slice.good())
               FileOpenFailed(HeaderFileName_slice);

           HeaderFile_slice.precision(17);

           HeaderFile_slice << m_N_slice_snapshots_ << "\n";
           HeaderFile_slice << m_dt_slice_snapshots_lab_ << "\n";
           HeaderFile_slice << m_gamma_boost_ << "\n";
           HeaderFile_slice << m_beta_boost_ << "\n";

        }

    }


}

LabFrameSnapShot::
LabFrameSnapShot(Real t_lab_in, Real t_boost, Real inv_gamma_boost_in,
                 Real inv_beta_boost_in, Real dz_lab_in, RealBox prob_domain_lab,
                 IntVect prob_ncells_lab, int ncomp_to_dump,
                 std::vector<std::string> mesh_field_names,
                 amrex::RealBox diag_domain_lab, Box diag_box, int file_num_in)
{
   m_t_lab = t_lab_in;
   m_dz_lab_ = dz_lab_in;
   m_inv_gamma_boost_ = inv_gamma_boost_in;
   m_inv_beta_boost_ = inv_beta_boost_in;
   m_prob_domain_lab_ = prob_domain_lab;
   m_prob_ncells_lab_ = prob_ncells_lab;
   m_diag_domain_lab_ = diag_domain_lab;
   m_buff_box_ = diag_box;
   m_ncomp_to_dump_ = ncomp_to_dump;
   m_mesh_field_names_ = mesh_field_names;
   m_file_num = file_num_in;
   m_current_z_lab = 0.0;
   m_current_z_boost = 0.0;
   updateCurrentZPositions(t_boost, m_inv_gamma_boost_, m_inv_beta_boost_);
   Real zmin_lab = m_prob_domain_lab_.lo(AMREX_SPACEDIM-1);
   m_initial_i = (m_current_z_lab - zmin_lab) / m_dz_lab_ ;
   m_file_name = Concatenate(WarpX::lab_data_directory + "/snapshots/snapshot",
                           m_file_num, 5);
   createLabFrameDirectories();
   m_buff_counter_ = 0;
   if (WarpX::do_back_transformed_fields) m_data_buffer_.reset(nullptr);
}

void
LabFrameDiag::
updateCurrentZPositions(Real t_boost, Real inv_gamma, Real inv_beta)
{
    m_current_z_boost = (m_t_lab*inv_gamma - t_boost)*PhysConst::c*inv_beta;
    m_current_z_lab   = (m_t_lab - t_boost*inv_gamma)*PhysConst::c*inv_beta;
}

void
LabFrameDiag::
createLabFrameDirectories() {
#ifdef WARPX_USE_HDF5
    if (ParallelDescriptor::IOProcessor())
    {
        output_create(m_file_name);
    }

    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        if (WarpX::do_back_transformed_fields)
        {
            const auto lo = lbound(m_buff_box_);
            for (int comp = 0; comp < m_ncomp_to_dump_; ++comp) {
                output_create_field(m_file_name, m_mesh_field_names_[comp],
                                    m_prob_ncells_lab_[0],
#if ( AMREX_SPACEDIM == 3 )
                                    m_prob_ncells_lab_[1],
#else
                                    1,
#endif
                                    m_prob_ncells_lab_[AMREX_SPACEDIM-1]+1);
            }
        }
    }

    ParallelDescriptor::Barrier();

    if (WarpX::do_back_transformed_particles){
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::vector<std::string> species_names = mypc.GetSpeciesNames();
        // Loop over species to be dumped to BFD
        for (int j = 0; j < mypc.nSpeciesBackTransformedDiagnostics(); ++j)
        {
            // Loop over species to be dumped to BFD
            std::string species_name =
                species_names[mypc.mapSpeciesBackTransformedDiagnostics(j)];
            output_create_species_group(m_file_name, species_name);
            for (int k = 0; k < static_cast<int>(particle_field_names.size()); ++k)
            {
                std::string field_path = species_name + "/" + particle_field_names[k];
                output_create_particle_field(m_file_name, field_path);
            }
        }
    }
#else
    if (ParallelDescriptor::IOProcessor()) {

        if (!UtilCreateDirectory(m_file_name, 0755))
            CreateDirectoryFailed(m_file_name);

        const int nlevels = 1;
        for(int i = 0; i < nlevels; ++i) {
            const std::string &fullpath = LevelFullPath(i, m_file_name);
            if (!UtilCreateDirectory(fullpath, 0755))
                CreateDirectoryFailed(fullpath);
        }

        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::vector<std::string> species_names = mypc.GetSpeciesNames();

        const std::string particles_prefix = "particle";
        // Loop over species to be dumped to BFD
        for(int i = 0; i < mypc.nSpeciesBackTransformedDiagnostics(); ++i) {
            // Get species name
            std::string species_name =
                species_names[mypc.mapSpeciesBackTransformedDiagnostics(i)];
            const std::string fullpath = m_file_name + "/" + species_name;
            if (!UtilCreateDirectory(fullpath, 0755))
                CreateDirectoryFailed(fullpath);
        }
    }
#endif
    ParallelDescriptor::Barrier();

    writeLabFrameHeader();
}

void
LabFrameDiag::
writeLabFrameHeader() {
#ifndef WARPX_USE_HDF5
    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(m_file_name + "/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if(!HeaderFile.good())
            FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << m_t_lab << "\n";
        // Write domain number of cells
        HeaderFile << m_prob_ncells_lab_[0] << ' '
#if ( AMREX_SPACEDIM==3 )
                   << m_prob_ncells_lab_[1] << ' '
#endif
                   << m_prob_ncells_lab_[AMREX_SPACEDIM-1] <<'\n';
        // Write domain physical boundaries
        // domain lower bound
        HeaderFile << m_diag_domain_lab_.lo(0) << ' '
#if ( AMREX_SPACEDIM==3 )
                   << m_diag_domain_lab_.lo(1) << ' '
#endif
                   << m_diag_domain_lab_.lo(AMREX_SPACEDIM-1) <<'\n';
        // domain higher bound
        HeaderFile << m_diag_domain_lab_.hi(0) << ' '
#if ( AMREX_SPACEDIM==3 )
                   << m_diag_domain_lab_.hi(1) << ' '
#endif
                   << m_diag_domain_lab_.hi(AMREX_SPACEDIM-1) <<'\n';
        // List of fields dumped to file
        for (int i=0; i<m_ncomp_to_dump_; i++)
        {
            HeaderFile << m_mesh_field_names_[i] << ' ';
        }
        HeaderFile << "\n";
    }
#endif

}


LabFrameSlice::
LabFrameSlice(Real t_lab_in, Real t_boost, Real inv_gamma_boost_in,
                 Real inv_beta_boost_in, Real dz_lab_in, RealBox prob_domain_lab,
                 IntVect prob_ncells_lab, int ncomp_to_dump,
                 std::vector<std::string> mesh_field_names,
                 RealBox diag_domain_lab, Box diag_box, int file_num_in,
                 amrex::Real particle_slice_dx_lab)
{
    m_t_lab = t_lab_in;
    m_dz_lab_ = dz_lab_in;
    m_inv_gamma_boost_ = inv_gamma_boost_in;
    m_inv_beta_boost_ = inv_beta_boost_in;
    m_prob_domain_lab_ = prob_domain_lab;
    m_prob_ncells_lab_ = prob_ncells_lab;
    m_diag_domain_lab_ = diag_domain_lab;
    m_buff_box_ = diag_box;
    m_ncomp_to_dump_ = ncomp_to_dump;
    m_mesh_field_names_ = mesh_field_names;
    m_file_num = file_num_in;
    m_current_z_lab = 0.0;
    m_current_z_boost = 0.0;
    updateCurrentZPositions(t_boost, m_inv_gamma_boost_, m_inv_beta_boost_);
    Real zmin_lab = m_prob_domain_lab_.lo(AMREX_SPACEDIM-1);
    m_initial_i = (m_current_z_lab - zmin_lab)/m_dz_lab_;
    m_file_name = Concatenate(WarpX::lab_data_directory+"/slices/slice",m_file_num,5);
    createLabFrameDirectories();
    m_buff_counter_ = 0;
    m_particle_slice_dx_lab_ = particle_slice_dx_lab;

    if (WarpX::do_back_transformed_fields) m_data_buffer_.reset(nullptr);
}

void
LabFrameSnapShot::
AddDataToBuffer( MultiFab& tmp, int k_lab,
                 amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump)
{
    const int ncomp_to_dump = map_actual_fields_to_dump.size();
    MultiFab& buf = *m_data_buffer_;
    for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Array4<Real> tmp_arr = tmp[mfi].array();
         Array4<Real> buf_arr = buf[mfi].array();
         // For 3D runs, tmp is a 2D (x,y) multifab that contains only
         // slice to write to file
         const Box& bx = mfi.tilebox();
         const auto field_map_ptr = map_actual_fields_to_dump.dataPtr();
         ParallelFor(bx, ncomp_to_dump,
             [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
             {
                 const int icomp = field_map_ptr[n];
#if (AMREX_SPACEDIM == 3)
                 buf_arr(i,j,k_lab,n) = tmp_arr(i,j,k,icomp);
#else
                 buf_arr(i,k_lab,k,n) = tmp_arr(i,j,k,icomp);
#endif
             }
         );
    }
}


void
LabFrameSlice::
AddDataToBuffer( MultiFab& tmp, int k_lab,
                 amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump)
{
    const int ncomp_to_dump = map_actual_fields_to_dump.size();
    MultiFab& buf = *m_data_buffer_;
    for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       Box& bx = m_buff_box_;
       const Box& bx_bf = mfi.tilebox();
       bx.setSmall(AMREX_SPACEDIM-1,bx_bf.smallEnd(AMREX_SPACEDIM-1));
       bx.setBig(AMREX_SPACEDIM-1,bx_bf.bigEnd(AMREX_SPACEDIM-1));
       Array4<Real> tmp_arr = tmp[mfi].array();
       Array4<Real> buf_arr = buf[mfi].array();
       const auto field_map_ptr = map_actual_fields_to_dump.dataPtr();
       ParallelFor(bx, ncomp_to_dump,
           [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
           {
              const int icomp = field_map_ptr[n];
#if (AMREX_SPACEDIM == 3)
              buf_arr(i,j,k_lab,n) = tmp_arr(i,j,k,icomp);
#else
              buf_arr(i,k_lab,k,n) = tmp_arr(i,j,k,icomp);
#endif
           });
    }

}


void
LabFrameSnapShot::
AddPartDataToParticleBuffer(
    Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer,
    int nspeciesBoostedFrame) {
    for (int isp = 0; isp < nspeciesBoostedFrame; ++isp) {
        auto np = tmp_particle_buffer[isp].GetRealData(DiagIdx::w).size();
        if (np == 0) continue;

        // allocate size of particle buffer array to np
        // This is a growing array. Each time we add np elements
        // to the existing array which has size = init_size
        const int init_size = m_particles_buffer_[isp].GetRealData(DiagIdx::w).size();
        const int total_size = init_size + np;
        m_particles_buffer_[isp].resize(total_size);

        // Data pointers to particle attributes //
        ParticleReal* const AMREX_RESTRICT wp_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::w).data();
        ParticleReal* const AMREX_RESTRICT x_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::x).data();
        ParticleReal* const AMREX_RESTRICT y_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::y).data();
        ParticleReal* const AMREX_RESTRICT z_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::z).data();
        ParticleReal* const AMREX_RESTRICT ux_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::ux).data();
        ParticleReal* const AMREX_RESTRICT uy_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::uy).data();
        ParticleReal* const AMREX_RESTRICT uz_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::uz).data();

        ParticleReal const* const AMREX_RESTRICT wp_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::w).data();
        ParticleReal const* const AMREX_RESTRICT x_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::x).data();
        ParticleReal const* const AMREX_RESTRICT y_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::y).data();
        ParticleReal const* const AMREX_RESTRICT z_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::z).data();
        ParticleReal const* const AMREX_RESTRICT ux_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::ux).data();
        ParticleReal const* const AMREX_RESTRICT uy_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::uy).data();
        ParticleReal const* const AMREX_RESTRICT uz_temp =
                    tmp_particle_buffer[isp].GetRealData(DiagIdx::uz).data();

        // copy all the particles from tmp to buffer
        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE(int i)
        {
            wp_buff[init_size + i] = wp_temp[i];
            x_buff[init_size + i]  = x_temp[i];
            y_buff[init_size + i]  = y_temp[i];
            z_buff[init_size + i]  = z_temp[i];
            ux_buff[init_size + i] = ux_temp[i];
            uy_buff[init_size + i] = uy_temp[i];
            uz_buff[init_size + i] = uz_temp[i];
        });
    }
}

void
LabFrameSlice::
AddPartDataToParticleBuffer(
    Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer,
    int nSpeciesBackTransformedDiagnostics) {


    for (int isp = 0; isp < nSpeciesBackTransformedDiagnostics; ++isp) {
        auto np = tmp_particle_buffer[isp].GetRealData(DiagIdx::w).size();

        if (np == 0) continue;

        ParticleReal const* const AMREX_RESTRICT wp_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::w).data();
        ParticleReal const* const AMREX_RESTRICT x_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::x).data();
        ParticleReal const* const AMREX_RESTRICT y_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::y).data();
        ParticleReal const* const AMREX_RESTRICT z_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::z).data();
        ParticleReal const* const AMREX_RESTRICT ux_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::ux).data();
        ParticleReal const* const AMREX_RESTRICT uy_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::uy).data();
        ParticleReal const* const AMREX_RESTRICT uz_temp =
             tmp_particle_buffer[isp].GetRealData(DiagIdx::uz).data();

        // temporary arrays to store copy_flag and copy_index
        // for particles that cross the reduced domain for diagnostics.
        amrex::Gpu::ManagedDeviceVector<int> FlagForPartCopy(np);
        amrex::Gpu::ManagedDeviceVector<int> IndexForPartCopy(np);

        int* const AMREX_RESTRICT Flag = FlagForPartCopy.dataPtr();
        int* const AMREX_RESTRICT IndexLocation = IndexForPartCopy.dataPtr();

        // Compute extent of the reduced domain +/- user-defined physical width
        Real const xmin = m_diag_domain_lab_.lo(0)-m_particle_slice_dx_lab_;
        Real const xmax = m_diag_domain_lab_.hi(0)+m_particle_slice_dx_lab_;
#if (AMREX_SPACEDIM == 3)
        Real const ymin = m_diag_domain_lab_.lo(1)-m_particle_slice_dx_lab_;
        Real const ymax = m_diag_domain_lab_.hi(1)+m_particle_slice_dx_lab_;
#endif

        //Flag particles that need to be copied if they are
        // within the reduced slice +/- user-defined physical width
        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE(int i)
        {
            Flag[i] = 0;
            if ( x_temp[i] >= (xmin) &&
                 x_temp[i] <= (xmax) ) {
#if (AMREX_SPACEDIM == 3)
               if (y_temp[i] >= (ymin) &&
                   y_temp[i] <= (ymax) )
#endif
               {
                   Flag[i] = 1;
               }
            }
        });

        // Call exclusive scan to obtain location indices using
        // flag values. These location indices are used to copy data
        // from src to dst when the copy-flag is set to 1.
        amrex::Gpu::exclusive_scan(Flag,Flag+np,IndexLocation);
        const int copy_size = IndexLocation[np-1] + Flag[np-1];
        const int init_size = m_particles_buffer_[isp].GetRealData(DiagIdx::w).size();
        const int total_reducedDiag_size = copy_size + init_size;

        // allocate array size for reduced diagnostic buffer array
        m_particles_buffer_[isp].resize(total_reducedDiag_size);

        // Data pointers to particle attributes //
        ParticleReal* const AMREX_RESTRICT wp_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::w).data();
        ParticleReal* const AMREX_RESTRICT x_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::x).data();
        ParticleReal* const AMREX_RESTRICT y_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::y).data();
        ParticleReal* const AMREX_RESTRICT z_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::z).data();
        ParticleReal* const AMREX_RESTRICT ux_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::ux).data();
        ParticleReal* const AMREX_RESTRICT uy_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::uy).data();
        ParticleReal* const AMREX_RESTRICT uz_buff =
              m_particles_buffer_[isp].GetRealData(DiagIdx::uz).data();

        // Selective copy of particle data from tmp array to reduced buffer
        // array on the GPU using the flag value and index location.
        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE(int i)
        {
            if (Flag[i] == 1)
            {
               const int loc = IndexLocation[i] + init_size;
               wp_buff[loc] = wp_temp[i];
               x_buff[loc]  = x_temp[i];
               y_buff[loc]  = y_temp[i];
               z_buff[loc]  = z_temp[i];
               ux_buff[loc] = ux_temp[i];
               uy_buff[loc] = uy_temp[i];
               uz_buff[loc] = uz_temp[i];
            }
        });

    }
}
