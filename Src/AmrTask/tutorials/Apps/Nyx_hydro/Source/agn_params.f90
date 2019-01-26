module agn_params_module

    use amrex_fort_module, only : rt => amrex_real

    real(rt), save :: l_merge
    logical, save :: cutoff_vel
    real(rt), save :: eps_rad, eps_coupling, T_min, bondi_boost, max_frac_removed, frac_kinetic, eps_kinetic

end module agn_params_module
