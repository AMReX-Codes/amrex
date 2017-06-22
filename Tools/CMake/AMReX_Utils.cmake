###################################################################
# Check if dir or file given by path exists and issue a warning or
#  error if not
##################################################################
function ( check_path  path  message_type )
   if ( EXISTS ${path} )
   else ()
      message(${message_type} ${path} " does not exist!")
   endif ( EXISTS ${path} )
endfunction ()


# ###################################################################
# # Set default options for reasonable beahvior
# ##################################################################
# macro ( set_default_options )

#    # Default to debug is no other build type specified
#    if ( NOT CMAKE_BUILD_TYPE )
#       set ( CMAKE_BUILD_TYPE "Debug" CACHE STRING
# 	 "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
# 	 FORCE )
#    endif ()
   
#    # Compile options:
#    option (ENABLE_FORTRAN_MPI "Enable Fortran MPI Communicator" ON)
#    option (ENABLE_FMG "Enable Fortran for MultiGrid Solver" OFF)
#    option (ENABLE_FBASELIB "Enable Fortran BaseLib" OFF)
#    option (ENABLE_CXX11 "Enable C++11" ON)
#    option (ENABLE_POSITION_INDEPENDENT_CODE "Compile with position-independent code enabled" OFF)

#    set (BL_SPACEDIM 3 CACHE INT "Dimension of AMReX build")
#    set (ENABLE_MPI 1 CACHE INT "Enable build with MPI")
#    set (ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
#    set (BL_PRECISION "DOUBLE" CACHE INT "Precision of AMReX build")
#    set (BL_USE_PARTICLES 0 CACHE INT "Include Particles classes in AMReX build")
#    set (ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
#    set (ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in AMReX build")
  


# endmacro ()
