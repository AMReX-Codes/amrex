
set(CBOXDIR C_BaseLib)
set(FBOXDIR F_BaseLib)
set(CBNDRYDIR C_BoundaryLib)
set(CAMRDIR C_AMRLib)
set(CMGDIR LinearSolvers/C_CellMG)
set(CTMGDIR LinearSolvers/C_TensorMG)
set(CFMGDIR LinearSolvers/C_to_F_MG)
set(FMGDIR LinearSolvers/F_MG)
set(CAMRDATADIR Extern/amrdata)

# Pick one place to put modules so everyone can see them
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files CACHE PATH "Folder for fortran module files")
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

set(CBOXLIB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${CBOXDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CBNDRYDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CAMRDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CMGDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CTMGDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CFMGDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CAMRDATADIR} ${BOXLIB_EXTRA_CXX_INCLUDE_PATH})

include(PreprocessBoxLibFortran)

add_subdirectory(${CBOXDIR})
add_subdirectory(${CBNDRYDIR})
add_subdirectory(${CAMRDIR})
add_subdirectory(${CMGDIR})
add_subdirectory(${CTMGDIR})

set(FBOXLIB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${FBOXDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${FMGDIR} ${BOXLIB_EXTRA_Fortran_INCLUDE_PATH})

add_subdirectory(${FBOXDIR})
add_subdirectory(${FMGDIR})

add_subdirectory(${CFMGDIR})

add_library(cboxlib $<TARGET_OBJECTS:box_c> $<TARGET_OBJECTS:box_cbndry> $<TARGET_OBJECTS:box_camr> $<TARGET_OBJECTS:box_cmg> $<TARGET_OBJECTS:box_ctmg> $<TARGET_OBJECTS:box_cfmg>)
add_library(cfboxlib $<TARGET_OBJECTS:box_cfmg>)
add_library(fboxlib $<TARGET_OBJECTS:box_f> $<TARGET_OBJECTS:box_fmg>)
add_install_library(cboxlib)
add_install_library(cfboxlib)
add_install_library(fboxlib)
SET_TARGET_PROPERTIES(cboxlib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(cfboxlib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(fboxlib PROPERTIES LINKER_LANGUAGE C)

add_subdirectory(${CAMRDATADIR})

add_subdirectory(Extern/gslib)

