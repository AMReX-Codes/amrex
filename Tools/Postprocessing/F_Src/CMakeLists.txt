set( FSRC plt_compare_diff_grids.f90 fextract.f90 fboxinfo.f90 )

foreach ( MAIN ${FSRC} )

   string(REPLACE ".f90" "" THIS_EXE ${MAIN} )

   add_executable( ${THIS_EXE} ${MAIN} )
   target_link_libraries(${THIS_EXE} fboxlib cboxlib fboxlib cfboxlib
      box_camrdata ${MPI_Fortran_LIBRARIES})
   add_install_binary(${THIS_EXE})

endforeach ()
