#include "WarpXBoostedFrameDiagnostic.H"

#include "WarpX_f.H"

using namespace amrex;

namespace {

    Box 
    getIndexBox(const RealBox& real_box, const Geometry& geom) {
        IntVect slice_lo, slice_hi;
    
        D_TERM(slice_lo[0]=floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0));,
               slice_lo[1]=floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1));,
               slice_lo[2]=floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2)););
        
        D_TERM(slice_hi[0]=floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0));,
               slice_hi[1]=floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1));,
               slice_hi[2]=floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2)););
        
        return Box(slice_lo, slice_hi) & geom.Domain();
    }    

    std::unique_ptr<MultiFab> allocateSlice(const MultiFab& cell_centered_data, 
                                            const Geometry& geom, Real z_coord,
                                            Array<int>& slice_to_full_ba_map) {
        
        // Get our slice and convert to index space
        RealBox real_slice = geom.ProbDomain();
        real_slice.setLo(2, z_coord);
        real_slice.setHi(2, z_coord);
        Box slice_box = getIndexBox(real_slice, geom);
        
        // define the multifab that stores slice
        BoxArray ba = cell_centered_data.boxArray();
        const DistributionMapping& dm = cell_centered_data.DistributionMap();

        std::vector< std::pair<int, Box> > isects;
        ba.intersections(slice_box, isects, false, 0);

        Array<Box> boxes;
        Array<int> procs;
        for (int i = 0; i < isects.size(); ++i) {
            procs.push_back(dm[isects[i].first]);
            boxes.push_back(isects[i].second);
            slice_to_full_ba_map.push_back(isects[i].first);
        }

        BoxArray slice_ba(&boxes[0], boxes.size());
        DistributionMapping slice_dmap(procs);
        std::unique_ptr<MultiFab> slice(new MultiFab(slice_ba, slice_dmap, 10, 0));
        return slice;
    }
    
    std::unique_ptr<MultiFab> getSliceData(const MultiFab& cell_centered_data, 
                                           const Geometry& geom, Real z_coord) {
        
        Array<int> slice_to_full_ba_map;
        std::unique_ptr<MultiFab> slice = allocateSlice(cell_centered_data, geom, z_coord,
                                                        slice_to_full_ba_map);
        
        // Fill the slice with sampled data
        const BoxArray& ba = cell_centered_data.boxArray();
        for (MFIter mfi(*slice); mfi.isValid(); ++mfi) {
            int slice_gid = mfi.index();
            int full_gid = slice_to_full_ba_map[slice_gid];
            
            const Box& slice_box = mfi.validbox();
            const Box& full_box  = ba[full_gid];
            const Box& tile_box  = mfi.tilebox();
            
            WRPX_FILL_SLICE(cell_centered_data[full_gid].dataPtr(),
                            full_box.loVect(), full_box.hiVect(),
                            (*slice)[slice_gid].dataPtr(),
                            slice_box.loVect(), slice_box.hiVect(),
                            tile_box.loVect(), tile_box.hiVect());
        }
        
        return slice;
    }    
}

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
writeLabFrameData(const MultiFab& cell_centered_data, const Geometry& geom, Real t_boost) {
    
    for (int i = 0; i < N_snapshots_; ++i) {
        snapshots_[i].updateCurrentZPositions(t_boost, 
                                              inv_gamma_boost_,
                                              inv_beta_boost_);
        
        // for each z position, fill a slice with the data.
        int i_lab = (snapshots_[i].current_z_lab - snapshots_[i].zmin_lab) / dz_lab_;
        std::unique_ptr<MultiFab> slice = getSliceData(cell_centered_data, geom,
                                                       snapshots_[i].current_z_boost);

        // transform it to the lab frame
        for (MFIter mfi(*slice); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            const Box& tile_box = mfi.tilebox();
            WRPX_LORENTZ_TRANSFORM((*slice)[mfi].dataPtr(), box.loVect(), box.hiVect(),
                                   tile_box.loVect(), tile_box.hiVect(),
                                   &gamma_boost_, &beta_boost_);
        }
        
        // and write it to disk.
        std::stringstream ss;
        ss << snapshots_[i].file_name << "/Level_0/" << Concatenate("slice", i_lab, 5);    
        VisMF::Write(*slice, ss.str());
    }
}

void
BoostedFrameDiagnostic::
writeMetaData() 
{
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
    
    const int nlevels = 1;
    const std::string level_prefix = "Level_";
    PreBuildDirectorHierarchy(file_name,
                              level_prefix, nlevels, true);
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
