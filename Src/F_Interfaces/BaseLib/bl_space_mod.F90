
module bl_space_module

#if (BL_SPACEDIM == 1)
  integer, parameter :: bl_num_dims = 1
#elif (BL_SPACEDIM == 2)
  integer, parameter :: bl_num_dims = 2
#else
  integer, parameter :: bl_num_dims = 3
#endif

end module bl_space_module
