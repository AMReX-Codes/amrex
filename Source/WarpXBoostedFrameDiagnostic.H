#ifndef WARPX_BoostedFrameDiagnostic_H_
#define WARPX_BoostedFrameDiagnostic_H_

#include <vector>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>

#include "WarpXConst.H"

class BoostedFrameDiagnostic {

    amrex::Real gamma_boost_;
    amrex::Real inv_gamma_boost_;
    amrex::Real beta_boost_;
    amrex::Real inv_beta_boost_;
    amrex::Real dz_lab_;
    amrex::Real inv_dz_lab_;
    amrex::Real dt_snapshots_lab_;
    amrex::Real dt_boost_;
    int  N_snapshots_;
    int  Nz_lab_;
    
    struct LabSnapShot {
        std::string file_name;
        amrex::Real t_lab;
        amrex::Real zmin_lab;
        amrex::Real zmax_lab;
        amrex::Real current_z_lab;
        amrex::Real current_z_boost;
        int  file_num;
        
        LabSnapShot(amrex::Real t_lab_in, amrex::Real zmin_lab_in, 
                    amrex::Real zmax_lab_in, int file_num_in) 
            : t_lab(t_lab_in),
              zmin_lab(zmin_lab_in),
              zmax_lab(zmax_lab_in),
              file_num(file_num_in) 
        {
            current_z_lab = 0.0;
            current_z_boost = 0.0;
            file_name = amrex::Concatenate("lab_data", file_num, 5);
        }
    };
    
    std::vector<LabSnapShot> snapshots_;

public:
    BoostedFrameDiagnostic(amrex::Real zmin_lab, amrex::Real zmax_lab, 
                           amrex::Real v_window_lab, amrex::Real dt_snapshots_lab,
                           int N_snapshots, amrex::Real gamma_boost,
                           amrex::Real dt_boost);
};

#endif
