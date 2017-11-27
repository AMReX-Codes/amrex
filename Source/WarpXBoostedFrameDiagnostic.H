#ifndef WARPX_BoostedFrameDiagnostic_H_
#define WARPX_BoostedFrameDiagnostic_H_

#include <vector>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include "WarpXConst.H"

///
/// BoostedFrameDiagnostic is for handling IO when running in a boosted 
/// frame of reference. Because of the relativity of simultaneity, events that
/// are synchronized in the simulation frame are not synchronized in the 
/// lab frame. Thus, at a given t_boost, we must write slices of data to 
/// multiple output files, each one corresponding to a given time in the lab frame.
///
class BoostedFrameDiagnostic {

    ///
    /// LabSnapShot stores metadata corresponding to a single time
    /// snapshot in the lab frame. The snapshot is written to disk 
    /// in the directory "file_name". zmin_lab, zmax_lab, and t_lab
    /// are all constant for a given snapshot. current_z_lab and 
    /// current_z_boost for each snapshot are updated as the 
    /// simulation time in the boosted frame advances.
    ///
    struct LabSnapShot {
        
        std::string file_name;
        amrex::Real t_lab;
        amrex::Real zmin_lab;
        amrex::Real zmax_lab;
        amrex::Real current_z_lab;
        amrex::Real current_z_boost;
        int file_num;

        LabSnapShot(amrex::Real t_lab_in, amrex::Real zmin_lab_in, 
                    amrex::Real zmax_lab_in, int file_num_in);
        
        ///
        /// This snapshot is at time t_lab, and the simulation is at time t_boost.
        /// The Lorentz transformation picks out one slice corresponding to both
        /// of those times, at position current_z_boost and current_z_lab in the 
        /// boosted and lab frames, respectively. 
        ///
        void updateCurrentZPositions(amrex::Real t_boost, amrex::Real inv_gamma,
                                     amrex::Real inv_beta);

        ///
        /// Write some useful metadata about this snapshot.
        ///
        void writeSnapShotHeader();
    };

    amrex::Real gamma_boost_;
    amrex::Real inv_gamma_boost_;
    amrex::Real beta_boost_;
    amrex::Real inv_beta_boost_;
    amrex::Real dz_lab_;
    amrex::Real inv_dz_lab_;
    amrex::Real dt_snapshots_lab_;
    amrex::Real dt_boost_;
    int N_snapshots_;
    int Nz_lab_;    
    int boost_direction_;

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > data_buffer_;
    int num_buffer_ = 32;
    int max_box_size_ = 64;
    amrex::Vector<int> buff_counter_;

    amrex::Vector<LabSnapShot> snapshots_;

public:

    BoostedFrameDiagnostic(amrex::Real zmin_lab, amrex::Real zmax_lab, 
                           amrex::Real v_window_lab, amrex::Real dt_snapshots_lab,
                           int N_snapshots, amrex::Real gamma_boost,
                           amrex::Real dt_boost, int boost_direction);
    
    void writeLabFrameData(const amrex::MultiFab& cell_centered_data,
                           const amrex::Geometry& geom, amrex::Real t_boost);
    
    void writeMetaData();
};

#endif
