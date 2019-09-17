#ifndef WARPX_BoostedFrameDiagnostic_H_
#define WARPX_BoostedFrameDiagnostic_H_

#include <vector>
#include <map>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Geometry.H>
#include <AMReX_CudaContainers.H>

#include "MultiParticleContainer.H"
#include "WarpXConst.H"

///
/// BoostedFrameDiagnostic is for handling IO when running in a boosted
/// frame of reference. Because of the relativity of simultaneity, events that
/// are synchronized in the simulation frame are not synchronized in the
/// lab frame. Thus, at a given t_boost, we must write slices of data to
/// multiple output files, each one corresponding to a given time in the lab frame.
///
class LabFrameDiag { 
   public:
   std::string file_name;
   amrex::Real t_lab;
   amrex::RealBox prob_domain_lab_;
   amrex::IntVect prob_ncells_lab_;
   amrex::RealBox diag_domain_lab_;
   amrex::Box buff_box_;

   amrex::Real current_z_lab;
   amrex::Real current_z_boost;
   amrex::Real inv_gamma_boost_;
   amrex::Real inv_beta_boost_;
   amrex::Real dz_lab_;
   amrex::Real dx_;
   amrex::Real dy_;

   int ncomp_to_dump_;
   std::vector<std::string> mesh_field_names_;
           
   int file_num;
   int initial_i;

   // For back-transformed diagnostics of grid fields, data_buffer_
   // stores a buffer of the fields in the lab frame (in a MultiFab, i.e.
   // with all box data etc.). When the buffer if full, dump to file.
   std::unique_ptr<amrex::MultiFab> data_buffer_;
   // particles_buffer_ is currently blind to refinement level.
   // particles_buffer_[j] is a WarpXParticleContainer::DiagnosticParticleData where
   // - j is the species number for the current diag 
   amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> particles_buffer_;
   // buff_counter_ is the number of z slices in data_buffer_
   int buff_counter_;
   int num_buffer_ = 256;
   int max_box_size = 256;
   ///
   /// This snapshot is at time t_lab, and the simulation is at time t_boost.
   /// The Lorentz transformation picks out one slice corresponding to both
   /// of those times, at position current_z_boost and current_z_lab in the 
   /// boosted and lab frames, respectively. 
   ///
   void updateCurrentZPositions(amrex::Real t_boost, amrex::Real inv_gamma,
                                amrex::Real inv_beta);

   void createLabFrameDirectories();
   ///
   /// Write some useful metadata about this snapshot.
   ///
   void writeLabFrameHeader();

   virtual int getId() {return 0;}
   virtual void AddDataToBuffer(amrex::MultiFab& tmp_slice_ptr, int i_lab,
                amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump){}; 
//   virtual void initializeBufferMF(int boost_direction, int i_lab); 
   virtual void AddPartDataToParticleBuffer( amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer, int nSpeciesBoostedFrame) {};

};

/// LabFrameSnapShot stores metadata corresponding to a single time
/// snapshot in the lab frame. The snapshot is written to disk 
/// in the directory "file_name". zmin_lab, zmax_lab, and t_lab
/// are all constant for a given snapshot. current_z_lab and 
/// current_z_boost for each snapshot are updated as the 
/// simulation time in the boosted frame advances.
///
class LabFrameSnapShot : public LabFrameDiag {
   public:
   int snapshot_id; 
   int getId() {return snapshot_id;}
   LabFrameSnapShot(amrex::Real t_lab_in, amrex::Real t_boost, 
                    amrex::Real inv_gamma_boost_, amrex::Real inv_beta_boost_, 
                    amrex::Real dz_lab_, amrex::RealBox prob_domain_lab, 
                    amrex::IntVect prob_ncells_lab, int ncomp_to_dump, 
                    std::vector<std::string> mesh_field_names, 
                    amrex::RealBox diag_domain_lab, 
                    amrex::Box diag_box, int file_num_in);
   void AddDataToBuffer( amrex::MultiFab& tmp_slice, int k_lab,
                amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump); 
   void AddPartDataToParticleBuffer( amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer, int nSpeciesBoostedFrame);
};

class LabFrameSlice : public LabFrameDiag {
   public:
   int slice_id;
   int getId() {return slice_id;}
   LabFrameSlice(amrex::Real t_lab_in, amrex::Real t_boost, 
                    amrex::Real inv_gamma_boost_, amrex::Real inv_beta_boost_, 
                    amrex::Real dz_lab_, amrex::RealBox prob_domain_lab, 
                    amrex::IntVect prob_ncells_lab, int ncomp_to_dump, 
                    std::vector<std::string> mesh_field_names, 
                    amrex::RealBox diag_domain_lab, 
                    amrex::Box diag_box, int file_num_in, 
                    amrex::Real cell_dx, amrex::Real cell_dy);
   void AddDataToBuffer( amrex::MultiFab& tmp_slice_ptr, int i_lab,
                amrex::Gpu::ManagedDeviceVector<int> map_actual_fields_to_dump); 
   void AddPartDataToParticleBuffer( amrex::Vector<WarpXParticleContainer::DiagnosticParticleData> tmp_particle_buffer, int nSpeciesBoostedFrame);
};

class BoostedFrameDiagnostic {

    amrex::Real gamma_boost_;
    amrex::Real inv_gamma_boost_;
    amrex::Real beta_boost_;
    amrex::Real inv_beta_boost_;
    amrex::Real dz_lab_;
    amrex::Real inv_dz_lab_;
    amrex::Real dt_snapshots_lab_;
    amrex::Real dt_boost_;
    int N_snapshots_;
    int boost_direction_;
    int N_slice_snapshots_;
    amrex::Real dt_slice_snapshots_lab_;

    int num_buffer_ = 256;
    int max_box_size_ = 256;

    std::vector<std::unique_ptr<LabFrameDiag> > LabFrameDiags_;
 
    void writeParticleData(const WarpXParticleContainer::DiagnosticParticleData& pdata,
                           const std::string& name, const int i_lab);

#ifdef WARPX_USE_HDF5
    void writeParticleDataHDF5(const WarpXParticleContainer::DiagnosticParticleData& pdata,
                               const std::string& name, const std::string& species_name);
#endif
public:

    BoostedFrameDiagnostic(amrex::Real zmin_lab, amrex::Real zmax_lab,
                           amrex::Real v_window_lab, amrex::Real dt_snapshots_lab,
                           int N_snapshots, amrex::Real dt_slice_snapshots_lab, 
                           int N_slice_snapshots, amrex::Real gamma_boost,
                           amrex::Real t_boost, amrex::Real dt_boost, int boost_direction,
                           const amrex::Geometry& geom, amrex::RealBox& slice_realbox);
    
    void Flush(const amrex::Geometry& geom);

    void writeLabFrameData(const amrex::MultiFab* cell_centered_data,
                           const MultiParticleContainer& mypc,
                           const amrex::Geometry& geom,
                           const amrex::Real t_boost, const amrex::Real dt);

    void writeMetaData();

private:
    // Map field names and component number in cell_centered_data
    std::map<std::string, int> possible_fields_to_dump = {
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
    std::vector<std::string> mesh_field_names = {"Ex", "Ey", "Ez",
                                                 "Bx", "By", "Bz",
                                                 "jx", "jy", "jz", "rho"};
    int ncomp_to_dump = 10;


};

#endif
