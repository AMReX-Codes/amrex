# -*- mode: cmake -*-------------------------------------------
project(GSLIB)

set(F77_source_files acorni.f backtr.f beyond.f blue.f chknam.f chktitle.f CMakeLists.txt cova3.f dlocate.f dpowint.f dsortem.f gauinv.f gcum.f getindx.f getz.f green.f hexa.f ksol.f ktsol.f locate.f nscore.f numtext.f ordrel.f picksupr.f powint.f psfill.f psline.f pstext.f rand.f red.f resc.f scal.f setrot.f setsupr.f sortem.f sqdist.f srchsupr.f strlen.f)

add_library(gslib ${F77_source_files})
set_target_properties(gslib PROPERTIES LINKER_LANGUAGE Fortran)

if (${CCSE_INSTALL_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
else ()
    install(TARGETS gslib
      LIBRARY DESTINATION ${CCSE_INSTALL_DIR}/lib${MPI_SUFFIX}${OMP_SUFFIX}
      ARCHIVE DESTINATION ${CCSE_INSTALL_DIR}/lib${MPI_SUFFIX}${OMP_SUFFIX}
      CONFIGURATIONS ${CMAKE_BUILD_TYPE}
      )
endif ()

if (BUILD_TESTS)
endif()

