/* Copyright 2019 Andrew Myers, Axel Huebl, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_BackTransformedDiagnostic_H_
#define WARPX_BackTransformedDiagnostic_H_

#include "Particles/MultiParticleContainer.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <map>
#include <vector>


/** \brief
 * The capability for back-transformed lab-frame data is implemented to generate
 * the full diagnostic snapshot for the entire domain and reduced diagnostic
 * (1D, 2D or 3D 'slices') for a sub-domain.
 * LabFrameDiag class defines the parameters required to backtrasform data from
 * boosted frame at (z_boost,t_boost) to lab-frame at (z_lab, t_lab) using Lorentz
 * transformation. This Lorentz transformation picks out one slice corresponding
 * to both of those times, at position current_z_boost and current_z_lab in the
 * boosted and lab frames, respectively.
 * Two derived classes, namely, LabFrameSnapShot and LabFrameSlice are defined to
 * store the full back-transformed diagnostic snapshot of the entire domain and
 * reduced back-transformed diagnostic for a sub-domain, respectively.
 * The approach here is to define an array of LabFrameDiag  which would include
 * both, full domain snapshots and reduced domain 'slices', sorted based on their
 * respective t_lab. This is done to re-use the backtransformed data stored in
 * the slice multifab at (z_lab,t_lab)
 * for the full domain snapshot and sub-domain slices that have the same t_lab,
 * instead of re-generating the backtransformed slice data at z_lab for a given
 * t_lab for each diagnostic.
 */
class LabFrameDiag {
    public:
    std::string m_file_name;
    amrex::Real m_t_lab;
    amrex::RealBox m_prob_domain_lab_;
    amrex::IntVect m_prob_ncells_lab_;
    amrex::RealBox m_diag_domain_lab_;
    amrex::Box m_buff_box_;

    amrex::Real m_current_z_lab;
    amrex::Real m_current_z_boost;
    amrex::Real m_inv_gamma_boost_;
    amrex::Real m_inv_beta_boost_;
    amrex::Real m_dz_lab_;
    amrex::Real m_particle_slice_dx_lab_;

    int m_ncomp_to_dump_;
    std::vector<std::string> m_mesh_field_names_;

    int m_file_num;
    int m_initial_i;

    // For back-transformed diagnostics of grid fields, data_buffer_
    // stores a buffer of the fields in the lab frame (in a MultiFab, i.e.
    // with all box data etc.). When the buffer if full, dump to file.
    std::unique_ptr<amrex::MultiFab> m_data_buffer_;
    // particles_buffer_ is currently blind to refinement level.
    // particles_buffer_[j] is a WarpXParticleContainer::DiagnosticParticleData
    //  where - j is the species number for the current diag
    amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> m_particles_buffer_;
    // buff_counter_ is the number of z slices in data_buffer_
    int m_buff_counter_;
    int m_num_buffer_ = 256;
    int m_max_box_size = 256;
    void updateCurrentZPositions(amrex::Real t_boost, amrex::Real inv_gamma,
                                 amrex::Real inv_beta);

    void createLabFrameDirectories();

    void writeLabFrameHeader();

    /// Back-transformed lab-frame field data is copied from
    /// tmp_slice to data_buffer where it is stored.
    /// For the full diagnostic, all the data in the MultiFab is copied.
    /// For the reduced diagnostic, the data is selectively copied if the
    /// extent of the z_lab multifab intersects with the user-defined sub-domain
    /// for the reduced diagnostic (i.e., a 1D, 2D, or 3D region of the domain).
    virtual void AddDataToBuffer(amrex::MultiFab& tmp_slice_ptr, int i_lab,
                 amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump){};

    /// Back-transformed lab-frame particles is copied from
    /// tmp_particle_buffer to particles_buffer.
    /// For the full diagnostic, all the particles are copied,
    /// while for the reduced diagnostic, particles are selectively
    /// copied if their position in within the user-defined
    /// sub-domain +/- 1 cell size width for the reduced slice diagnostic.
    virtual void AddPartDataToParticleBuffer(
                 amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer,
                 int nSpeciesBoostedFrame) {};
};

/** \brief
 * LabFrameSnapShot stores the back-transformed lab-frame metadata
 * corresponding to a single time snapshot of the full domain.
 * The snapshot data is written to disk in the directory lab_frame_data/snapshots/.
 * zmin_lab, zmax_lab, and t_lab are all constant for a given snapshot.
 * current_z_lab and current_z_boost for each snapshot are updated as the
 * simulation time in the boosted frame advances.
 */

class LabFrameSnapShot : public LabFrameDiag {
    public:
    LabFrameSnapShot(amrex::Real t_lab_in, amrex::Real t_boost,
                     amrex::Real inv_gamma_boost_in, amrex::Real inv_beta_boost_in,
                     amrex::Real dz_lab_in, amrex::RealBox prob_domain_lab,
                     amrex::IntVect prob_ncells_lab, int ncomp_to_dump,
                     std::vector<std::string> mesh_field_names,
                     amrex::RealBox diag_domain_lab,
                     amrex::Box diag_box, int file_num_in);
    void AddDataToBuffer( amrex::MultiFab& tmp_slice, int k_lab,
         amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump) override;
    void AddPartDataToParticleBuffer(
         amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer,
         int nSpeciesBoostedFrame) override;
};


/** \brief
 * LabFrameSlice stores the back-transformed metadata corresponding to a single time at the
 * user-defined slice location. This could be a 1D line, 2D slice, or 3D box
 * (a reduced back-transformed diagnostic) within the computational
 * domain, as defined in the input file by the user. The slice is written to disk in the
 * lab_frame_data/slices. Similar to snapshots, zmin_lab, zmax_lab, and
 * t_lab are constant for a given slice. current_z_lab and current_z_boost
 * for each snapshot are updated as the sim time in boosted frame advances.
 */
class LabFrameSlice : public LabFrameDiag {
    public:
    LabFrameSlice(amrex::Real t_lab_in, amrex::Real t_boost,
                     amrex::Real inv_gamma_boost_in, amrex::Real inv_beta_boost_in,
                     amrex::Real dz_lab_in, amrex::RealBox prob_domain_lab,
                     amrex::IntVect prob_ncells_lab, int ncomp_to_dump,
                     std::vector<std::string> mesh_field_names,
                     amrex::RealBox diag_domain_lab,
                     amrex::Box diag_box, int file_num_in,
                     amrex::Real particle_slice_dx_lab);
    void AddDataToBuffer( amrex::MultiFab& tmp_slice_ptr, int i_lab,
         amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump) override;
    void AddPartDataToParticleBuffer(
         amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer,
         int nSpeciesBoostedFrame) override;
};

/** \brief
 * BackTransformedDiagnostic class handles the back-transformation of data when
 * running simulations in a boosted frame of reference to the lab-frame.
 * Because of the relativity of simultaneity, events that are synchronized
 * in the simulation boosted frame are not
 * synchronized in the lab frame. Thus, at a given t_boost, we must write
 * slices of back-transformed data to multiple output files, each one
 * corresponding to a given time in the lab frame. The member function
 * writeLabFrameData() orchestrates the operations required to
 * Lorentz-transform data from boosted-frame to lab-frame and store them
 * in the LabFrameDiag class, which writes out the field and particle data
 * to the output directory. The functions Flush() and writeLabFrameData()
 * are called at the end of the simulation and when the
 * the buffer for data storage is full, respectively. The particle data
 * is collected and written only if particle.do_back_transformed_diagnostics = 1.
 */
class BackTransformedDiagnostic {

public:

    BackTransformedDiagnostic(amrex::Real zmin_lab, amrex::Real zmax_lab,
                           amrex::Real v_window_lab, amrex::Real dt_snapshots_lab,
                           int N_snapshots, amrex::Real dt_slice_snapshots_lab,
                           int N_slice_snapshots, amrex::Real gamma_boost,
                           amrex::Real t_boost, amrex::Real dt_boost,
                           int boost_direction, const amrex::Geometry& geom,
                           amrex::RealBox& slice_realbox,
                           amrex::Real particle_slice_width_lab);

    /// Flush() is called at the end of the simulation when the buffers that contain
    /// back-transformed lab-frame data even if they are not full.
    void Flush(const amrex::Geometry& geom);

    /** \brief
     * The order of operations performed in writeLabFrameData is as follows :
     * 1. Loops over the sorted back-transformed diags and for each diag
     *    steps 2-7 are performed
     * 2. Based on t_lab and t_boost, obtain z_lab and z_boost.
     * 3. Define data_buffer multifab that will store the data in the BT diag.
     * 4. Define slice multifab at z_index that corresponds to z_boost and
     *    getslicedata using cell-centered data at z_index and its distribution map.
     * 5. Lorentz transform data stored in slice from z_boost,t_Boost to z_lab,t_lab
     *    and store in slice multifab.
     * 6. Generate a temporary slice multifab with distribution map of lab-frame
     *    data but at z_boost and
     *    ParallelCopy data from the slice multifab to the temporary slice.
     * 7. Finally, AddDataToBuffer is called where the data from temporary slice
     *    is simply copied from tmp_slice(i,j,k_boost) to
     *    LabFrameDiagSnapshot(i,j,k_lab) for full BT lab-frame diagnostic
     *    OR from tmp_slice(i,j,k_boost) to
     *    LabFrameDiagSlice(i,j,k_lab) for the reduced slice diagnostic
     * 8. Similarly, particles that crossed the z_boost plane are selected
     *    and lorentz-transformed to the lab-frame and copied to the full
     *    and reduce diagnostic and stored in particle_buffer.
     */
    void writeLabFrameData(const amrex::MultiFab* cell_centered_data,
                           const MultiParticleContainer& mypc,
                           const amrex::Geometry& geom,
                           const amrex::Real t_boost, const amrex::Real dt);
    /// The metadata containg information on t_boost, num_snapshots, and Lorentz parameters.
    void writeMetaData();

private:
    amrex::Real m_gamma_boost_;
    amrex::Real m_inv_gamma_boost_;
    amrex::Real m_beta_boost_;
    amrex::Real m_inv_beta_boost_;
    amrex::Real m_dz_lab_;
    amrex::Real m_inv_dz_lab_;
    amrex::Real m_dt_snapshots_lab_;
    amrex::Real m_dt_boost_;
    int m_N_snapshots_;
    int m_boost_direction_;
    int m_N_slice_snapshots_;
    amrex::Real m_dt_slice_snapshots_lab_;
    amrex::Real m_particle_slice_width_lab_;

    int m_num_buffer_ = 256;
    int m_max_box_size_ = 256;

    std::vector<std::unique_ptr<LabFrameDiag> > m_LabFrameDiags_;

    void writeParticleData(
         const WarpXParticleContainer::DiagnosticParticleData& pdata,
         const std::string& name, const int i_lab);

#ifdef WARPX_USE_HDF5
    void writeParticleDataHDF5(
         const WarpXParticleContainer::DiagnosticParticleData& pdata,
         const std::string& name, const std::string& species_name);
#endif
    // Map field names and component number in cell_centered_data
    std::map<std::string, int> m_possible_fields_to_dump = {
        {"Ex" , 0},
        {"Ey" , 1},
        {"Ez" , 2},
        {"Bx" , 3},
        {"By" , 4},
        {"Bz" , 5},
        {"jx" , 6},
        {"jy" , 7},
        {"jz" , 8},
        {"rho", 9} };

    // maps field index in data_buffer_[i] -> cell_centered_data for
    // snapshots i. By default, all fields in cell_centered_data are dumped.
    // Needs to be amrex::Vector because used in a ParallelFor kernel.
    amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump;
    // Name of fields to dump. By default, all fields in cell_centered_data.
    // Needed for file headers only.
    std::vector<std::string> m_mesh_field_names = {"Ex", "Ey", "Ez",
                                                 "Bx", "By", "Bz",
                                                 "jx", "jy", "jz", "rho"};
    int m_ncomp_to_dump = 10;


};

#endif
