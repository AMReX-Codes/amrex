set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND ALLHEADERS AMReX_EB2_F.H           AMReX_EB2_IF_Box.H        AMReX_EB2_IF_Lathe.H )
list ( APPEND ALLHEADERS AMReX_EB2_IF_Union.H    AMReX_EB2_GeometryShop.H  AMReX_EB2_IF_Complement.H )
list ( APPEND ALLHEADERS AMReX_EB2_IF_Plane.H    AMReX_EB2_IndexSpaceI.H   AMReX_EB2_Graph.H )
list ( APPEND ALLHEADERS AMReX_EB2_IF_Cylinder.H AMReX_EB2_IF_Scale.H      AMReX_EB2_Level.H )
list ( APPEND ALLHEADERS AMReX_EB2.H             AMReX_EB2_IF_Ellipsoid.H  AMReX_EB2_IF_Sphere.H )
list ( APPEND ALLHEADERS AMReX_EB2_MultiGFab.H   AMReX_EB2_IF_AllRegular.H AMReX_EB2_IF_Intersection.H )
list ( APPEND AMReX_EB2_IF_Translation.H )

list ( APPEND CXXSRC  AMReX_EB2.cpp  AMReX_EB2_Level.cpp  AMReX_EB2_MultiGFab.cpp )

list ( APPEND F90SRC  AMReX_eb2_${DIM}d.F90 )

#
# This is an option to display arrays in GDB. It is for debugging
# AMReX and currently not supported by cmake.# 
# if(USE_ARRAYVIEW EQUAL TRUE)
#   set(CXX_source_files ${CXX_source_files} AMReX_DatasetClient.cpp)
#   set(CXX_header_files ${CXX_header_files} AMReX_DatasetClient.H)
# endif()

# 
#  Collect sources
# 
set ( ALLSRC  ${CXXSRC} ${F90SRC} ${F77SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR} )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})
