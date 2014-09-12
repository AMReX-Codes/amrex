module prob_module

    use bl_types

    ! Here we hold the specific problem parameters
    real(dp_t), parameter :: uadv = 1.d0
    real(dp_t), parameter :: vadv = 0.d0
    real(dp_t), parameter :: wadv = 0.d0
    real(dp_t), parameter ::   mu = 2.d-2

end module prob_module
