#include "WarpXBoostedFrameDiagnostic.H"
#include <AMReX_MultiFabUtil.H>
#include "WarpX_f.H"

using namespace amrex;

BoostedFrameDiagnostic::
BoostedFrameDiagnostic(Real zmin_lab, Real zmax_lab, Real v_window_lab,
                       Real dt_snapshots_lab, int N_snapshots, 
                       Real gamma_boost, Real dt_boost, 
                       int boost_direction)
    : gamma_boost_(gamma_boost),
      dt_snapshots_lab_(dt_snapshots_lab),
      dt_boost_(dt_boost),
      N_snapshots_(N_snapshots),
      boost_direction_(boost_direction)
{

    BL_PROFILE("BoostedFrameDiagnostic::BoostedFrameDiagnostic");

    inv_gamma_boost_ = 1.0 / gamma_boost_;
    beta_boost_ = std::sqrt(1.0 - inv_gamma_boost_*inv_gamma_boost_);
    inv_beta_boost_ = 1.0 / beta_boost_;
    
    dz_lab_ = PhysConst::c * dt_boost_ * inv_beta_boost_ * inv_gamma_boost_;
    inv_dz_lab_ = 1.0 / dz_lab_;
    Nz_lab_ = static_cast<int>((zmax_lab - zmin_lab) * inv_dz_lab_);

    writeMetaData();

    data_buffer_.resize(N_snapshots);
    for (int i = 0; i < N_snapshots; ++i) {
        Real t_lab = i * dt_snapshots_lab_;
        LabSnapShot snapshot(t_lab, zmin_lab + v_window_lab * t_lab,
                             zmax_lab + v_window_lab * t_lab, i);
        snapshots_.push_back(snapshot);
        buff_counter_.push_back(0);
        data_buffer_[i].reset( nullptr );
    }

    AMREX_ALWAYS_ASSERT(max_box_size_ >= num_buffer_);
}

void
BoostedFrameDiagnostic::
writeLabFrameData(const MultiFab& cell_centered_data, const Geometry& geom, Real t_boost) {

    BL_PROFILE("BoostedFrameDiagnostic::writeLabFrameData");

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    const RealBox& domain_z_boost = geom.ProbDomain();
    const Real zlo_boost = domain_z_boost.lo(boost_direction_);
    const Real zhi_boost = domain_z_boost.hi(boost_direction_);

    for (int i = 0; i < N_snapshots_; ++i) {
        const Real old_z_boost = snapshots_[i].current_z_boost;
        snapshots_[i].updateCurrentZPositions(t_boost,
                                              inv_gamma_boost_,
                                              inv_beta_boost_);
        
        if ( (snapshots_[i].current_z_boost < zlo_boost) or
             (snapshots_[i].current_z_boost > zhi_boost) or
             (snapshots_[i].current_z_lab < snapshots_[i].zmin_lab) or
             (snapshots_[i].current_z_lab > snapshots_[i].zmax_lab) ) continue;

        // for each z position, fill a slice with the data.
        int i_lab = (snapshots_[i].current_z_lab - snapshots_[i].zmin_lab) / dz_lab_;

        const int ncomp = cell_centered_data.nComp();
        const int start_comp = 0;
        std::unique_ptr<MultiFab> slice = amrex::get_slice_data(boost_direction_,
                                                                snapshots_[i].current_z_boost,
                                                                cell_centered_data, geom,
                                                                start_comp, ncomp);

        // transform it to the lab frame
        for (MFIter mfi(*slice); mfi.isValid(); ++mfi) {
            const Box& tile_box = mfi.tilebox();
            WRPX_LORENTZ_TRANSFORM_Z(BL_TO_FORTRAN_ANYD((*slice)[mfi]),
                                     BL_TO_FORTRAN_BOX(tile_box),
                                     &gamma_boost_, &beta_boost_);
        }

        if (buff_counter_[i] == 0) {
            Box buff_box = geom.Domain();
            buff_box.setSmall(boost_direction_, i_lab - num_buffer_ + 1);
            buff_box.setBig(boost_direction_, i_lab);
            BoxArray buff_ba(buff_box);
            buff_ba.maxSize(max_box_size_);
            DistributionMapping buff_dm(buff_ba);
            data_buffer_[i].reset( new MultiFab(buff_ba, buff_dm, ncomp, 0) );
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

        ++buff_counter_[i];

        if (buff_counter_[i] == num_buffer_) {
            std::stringstream ss;
            ss << snapshots_[i].file_name << "/Level_0/" << Concatenate("buffer", i_lab, 5);
            VisMF::Write(*data_buffer_[i], ss.str());
            buff_counter_[i] = 0;
        }
    }

    VisMF::SetHeaderVersion(current_version);
}

void
BoostedFrameDiagnostic::
writeMetaData() 
{

    BL_PROFILE("BoostedFrameDiagnostic::writeMetaData");

    if (ParallelDescriptor::IOProcessor()) {
        std::string DiagnosticDirectory = "lab_frame_data";
        
        if (!UtilCreateDirectory(DiagnosticDirectory, 0755))
            CreateDirectoryFailed(DiagnosticDirectory);
        
        std::string HeaderFileName(DiagnosticDirectory + "/Header");
        std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
                                 std::ofstream::trunc |
                                 std::ofstream::binary);
        if(!HeaderFile.good())
            FileOpenFailed(HeaderFileName);
        
        HeaderFile.precision(17);
        
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
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
LabSnapShot(Real t_lab_in, Real zmin_lab_in, 
            Real zmax_lab_in, int file_num_in) 
    : t_lab(t_lab_in),
      zmin_lab(zmin_lab_in),
      zmax_lab(zmax_lab_in),
      file_num(file_num_in) 
{
    current_z_lab = 0.0;
    current_z_boost = 0.0;
    file_name = Concatenate("lab_frame_data/snapshot", file_num, 5);
    
    if (ParallelDescriptor::IOProcessor()) {

        const int nlevels = 1;
        const std::string level_prefix = "Level_";
        
        if (!UtilCreateDirectory(file_name, 0755))
            CreateDirectoryFailed(file_name);
        for(int i(0); i < nlevels; ++i) {
            const std::string &fullpath = LevelFullPath(i, file_name);
            if (!UtilCreateDirectory(fullpath, 0755))
                CreateDirectoryFailed(fullpath);
        }
    }
    
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) writeSnapShotHeader();
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
    if (ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName(file_name + "/Header");
        std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
                                 std::ofstream::trunc |
                                 std::ofstream::binary);
        if(!HeaderFile.good())
            FileOpenFailed(HeaderFileName);
        
        HeaderFile.precision(17);
        
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        
        HeaderFile << t_lab << "\n";
        HeaderFile << zmin_lab << "\n";
        HeaderFile << zmax_lab << "\n";
    }
}
