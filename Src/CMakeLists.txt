#
# Set subdirectories path
# 
set (CBOXDIR Base)
set (FBOXDIR F_BaseLib)
set (CBNDRYDIR Boundary)
set (CPARTDIR Particle)
set (CAMRCOREDIR AmrCore)
set (CAMRDIR Amr)
set (CMGDIR LinearSolvers/C_CellMG)
set (CTMGDIR LinearSolvers/C_TensorMG)
set (CFMGDIR LinearSolvers/C_to_F_MG)
set (FMGDIR LinearSolvers/F_MG)
set (CAMRDATADIR Extern/amrdata)
set (FINTDIR F_Interfaces)

#
# Add all the subdir to the include path
#  
include_directories ( ${AMREX_EXTRA_Fortran_INCLUDE_PATH} ${AMREX_EXTRA_CXX_INCLUDE_PATH} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CBOXDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CBNDRYDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CPARTDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CAMRDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CAMRCOREDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CMGDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CTMGDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CFMGDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${CAMRDATADIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${FBOXDIR} )
include_directories ( ${INCLUDE_PATHS} )

list_include_dirs ( INCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/${FMGDIR} )
include_directories ( ${INCLUDE_PATHS} )

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)


add_subdirectory(${CBOXDIR})
add_subdirectory(${CBNDRYDIR})
add_subdirectory(${CPARTDIR})
add_subdirectory(${CAMRDIR})
add_subdirectory(${CAMRCOREDIR})
add_subdirectory(${CMGDIR})
add_subdirectory(${CTMGDIR})
add_subdirectory(${FBOXDIR})
add_subdirectory(${FMGDIR})
add_subdirectory(${CFMGDIR})


add_library(cboxlib $<TARGET_OBJECTS:box_c> $<TARGET_OBJECTS:box_cbndry> $<TARGET_OBJECTS:box_cpart> $<TARGET_OBJECTS:box_camr> $<TARGET_OBJECTS:box_camrcore> $<TARGET_OBJECTS:box_cmg> $<TARGET_OBJECTS:box_ctmg> $<TARGET_OBJECTS:box_cfmg>)
add_library(cfboxlib $<TARGET_OBJECTS:box_cfmg>)
add_library(fboxlib $<TARGET_OBJECTS:box_f> $<TARGET_OBJECTS:box_fmg>)
add_install_library(cboxlib)
add_install_library(cfboxlib)
add_install_library(fboxlib)

# Move Tools dir to install location
install(CODE "MESSAGE(\"-- Copying   : ${AMREX_TOOLS_DIR} --> ${CCSE_INSTALL_DIR}\")" )
install(CODE "FILE( COPY ${AMREX_TOOLS_DIR} DESTINATION ${CCSE_INSTALL_DIR} )")
SET_TARGET_PROPERTIES(cboxlib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(cfboxlib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(fboxlib PROPERTIES LINKER_LANGUAGE C)

add_subdirectory(${CAMRDATADIR})


