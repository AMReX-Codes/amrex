module comoving_module

    use constants_module, only : rt => type_real, M_PI

    integer,  save :: comoving_type
    real(rt), save :: comoving_OmM, comoving_OmB, comoving_OmN, comoving_h

end module comoving_module
