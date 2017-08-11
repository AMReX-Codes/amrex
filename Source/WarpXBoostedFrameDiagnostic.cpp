#include "WarpXBoostedFrameDiagnostic.H"

BoostedFrameDiagnostic::
BoostedFrameDiagnostic(Real zmin_lab, Real zmax_lab, Real v_window_lab,
                       Real dt_snapshots_lab, int N_snapshots, 
                       Real gamma_boost, Real dt_boost)
    : gamma_boost_(gamma_boost),
      dt_snapshots_lab_(dt_snapshots_lab),
      dt_boost_(dt_boost),
      N_snapshots_(N_snapshots)
{
    inv_gamma_boost_ = 1.0 / gamma_boost_;
    beta_boost_ = std::sqrt(1.0 - inv_gamma_boost_*inv_gamma_boost_);
    inv_beta_boost_ = 1.0 / beta_boost_;
    
    dz_lab_ = PhysConst::c * dt_boost_ * inv_beta_boost_ * inv_gamma_boost_;
    inv_dz_lab_ = 1.0 / dz_lab_;
    Nz_lab_ = static_cast<int>((zmax_lab - zmin_lab) * inv_dz_lab_);

    writeMetaData();

    for (int i = 0; i < N_snapshots; ++i) {
        Real t_lab = i * dt_snapshots_lab_;
        LabSnapShot snapshot(t_lab, zmin_lab + v_window_lab * t_lab,
                             zmax_lab + v_window_lab * t_lab, i);
        snapshots_.push_back(snapshot);
    }
}

void
BoostedFrameDiagnostic::
writeMetaData() 
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string DiagnosticDirectory = "lab_frame_data";
        
        if (!amrex::UtilCreateDirectory(DiagnosticDirectory, 0755))
            amrex::CreateDirectoryFailed(DiagnosticDirectory);
        
        std::string HeaderFileName(DiagnosticDirectory + "/Header");
        std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
                                 std::ofstream::trunc |
                                 std::ofstream::binary);
        if(!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
        
        HeaderFile.precision(17);
        
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        
        HeaderFile << N_snapshots_ << "\n";
        HeaderFile << dt_snapshots_lab_ << "\n";    
        HeaderFile << gamma_boost_ << "\n";
        HeaderFile << beta_boost_ << "\n";
        HeaderFile << dz_lab_ << "\n";
        HeaderFile << Nz_lab_ << "\n";
    }
}

BoostedFrameDiagnostic::LabSnapShot::
LabSnapShot(amrex::Real t_lab_in, amrex::Real zmin_lab_in, 
            amrex::Real zmax_lab_in, int file_num_in) 
    : t_lab(t_lab_in),
      zmin_lab(zmin_lab_in),
      zmax_lab(zmax_lab_in),
      file_num(file_num_in) 
{
    current_z_lab = 0.0;
    current_z_boost = 0.0;
    file_name = amrex::Concatenate("lab_frame_data/snapshot", file_num, 5);
    
    const int nlevels = 1;
    const std::string level_prefix = "Level_";
    amrex::PreBuildDirectorHierarchy(file_name,
                                     level_prefix, nlevels, true);
    writeSnapShotHeader();
}

void
BoostedFrameDiagnostic::LabSnapShot::
updateCurrentZPositions(amrex::Real t_boost, amrex::Real inv_gamma,
                        amrex::Real inv_beta) {
    current_z_boost = (t_lab*inv_gamma - t_boost)*PhysConst::c*inv_beta;
    current_z_lab =   (t_lab - t_boost*inv_gamma)*PhysConst::c*inv_beta;
}

void
BoostedFrameDiagnostic::LabSnapShot::
writeSnapShotHeader() {
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName(file_name + "/Header");
        std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
                                 std::ofstream::trunc |
                                 std::ofstream::binary);
        if(!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
        
        HeaderFile.precision(17);
        
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        
        HeaderFile << t_lab << "\n";
        HeaderFile << zmin_lab << "\n";
        HeaderFile << zmax_lab << "\n";
    }
}
