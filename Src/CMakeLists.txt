
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
set(CMAKE_FORTRAN_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${FBOXDIR} CACHE PATH "Folder for fortran module files")

set(CBOXLIB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${CBOXDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CBNDRYDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CAMRDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CMGDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CTMGDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CFMGDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${CAMRDATADIR} ${BOXLIB_EXTRA_INCLUDE_PATH})

include(PreprocessBoxLibFortran)

add_subdirectory(${CBOXDIR})
add_subdirectory(${CBNDRYDIR})
add_subdirectory(${CAMRDIR})
add_subdirectory(${CMGDIR})
add_subdirectory(${CTMGDIR})

set(FBOXLIB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${FBOXDIR} ${CMAKE_CURRENT_SOURCE_DIR}/${FMGDIR})

add_subdirectory(${FBOXDIR})
add_subdirectory(${FMGDIR})

add_subdirectory(${CFMGDIR})
add_subdirectory(${CAMRDATADIR})

add_subdirectory(Extern/gslib)

