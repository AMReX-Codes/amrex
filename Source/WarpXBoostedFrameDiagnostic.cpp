#include "WarpXBoostedFrameDiagnostic.H"

using namespace amrex;

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

Box getIndexBox(const RealBox& real_box, const Geometry& geom) {
    IntVect slice_lo, slice_hi;
    
    D_TERM(slice_lo[0]=floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0));,
           slice_lo[1]=floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1));,
           slice_lo[2]=floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2)););
    
    D_TERM(slice_hi[0]=floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0));,
           slice_hi[1]=floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1));,
           slice_hi[2]=floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2)););
    
    return Box(slice_lo, slice_hi) & geom.Domain();
}

void getSliceData(const MultiFab& cell_centered_data, 
                  const Geometry& geom, Real z_coord) {

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
    }  
    procs.push_back(ParallelDescriptor::MyProc());
    BoxArray slice_ba(&boxes[0], boxes.size());
    DistributionMapping slice_dmap(procs);
    MultiFab slice(slice_ba, slice_dmap, 6, 0);
    
    // Fill the slice with sampled data
    const Real* dx      = geom.CellSize();
    for (MFIter mfi(slice); mfi.isValid(); ++mfi) {
        int slice_gid = mfi.index();
        Box grid = slice_ba[slice_gid];
        int xlo = grid.smallEnd(0), ylo = grid.smallEnd(1), zlo = grid.smallEnd(2);
        int xhi = grid.bigEnd(0), yhi = grid.bigEnd(1), zhi = grid.bigEnd(2);
    
        IntVect iv = grid.smallEnd();
        ba.intersections(Box(iv, iv), isects, true, 0);
        BL_ASSERT(!isects.empty());
        int full_gid = isects[0].first;
    
        for (int k = zlo; k <= zhi; k++) {
            for (int j = ylo; j <= yhi; j++) {
                for (int i = xlo; i <= xhi; i++) {
                    for (int comp = 0; comp < 6; comp++) {
                        Real x = geom.ProbLo(0) + i*dx[0];
                        Real y = geom.ProbLo(1) + j*dx[1];
                        Real z = z_coord;
                            
                        D_TERM(iv[0]=floor((x - geom.ProbLo(0))/geom.CellSize(0));,
                               iv[1]=floor((y - geom.ProbLo(1))/geom.CellSize(1));,
                               iv[2]=floor((z - geom.ProbLo(2))/geom.CellSize(2)););
                            
                        slice[slice_gid](IntVect(i, j, k), comp) = cell_centered_data[full_gid](iv, comp);
                    }
                }
            }
        }
    }
}

void
BoostedFrameDiagnostic::
writeLabFrameData(const MultiFab& cell_centered_data, Real t_boost)
{
    
    for (int i = 0; i < N_snapshots_; ++i) {
        snapshots_[i].updateCurrentZPositions(t_boost, 
                                              inv_gamma_boost_,
                                              inv_beta_boost_);
        
        // for each z position, fill a slice with the data.
        int i_lab = (snapshots_[i].current_z_lab - snapshots_[i].zmin_lab) / dz_lab_;
        std::cout << i_lab << " ";
        std::cout << snapshots_[i].current_z_lab << " ";
        std::cout << snapshots_[i].zmin_lab << " ";
        std::cout << dz_lab_ << " ";
        std::cout<< snapshots_[i].current_z_boost << std::endl;
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
