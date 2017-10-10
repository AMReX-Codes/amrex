#
# Take care of subdir amrdata ( this can be enabled separatly or
# implicitly by enabling profparser )
# 
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND ALLHEADERS AMReX_AmrData.H   AMReX_XYPlotDataList.H   AMReX_DataServices.H   AMReX_AmrvisConstants.H )
list ( APPEND CXXSRC     AMReX_AmrData.cpp AMReX_XYPlotDataList.cpp AMReX_DataServices.cpp )

if (ENABLE_PROPARSER)
   list ( APPEND CXXSRC     AMReX_XYPlotDataList.cpp )
   list ( APPEND ALLHEADERS AMReX_XYPlotDataList.H   )
endif ()

list ( APPEND F77SRC     AMReX_FABUTIL_${DIM}D.F )

# Accumulate sources
set ( ALLSRC ${CXXSRC} ${F77SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/amrdata )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/amrdata )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/amrdata)


# ----------------------------------------------------------------------#

#
# The next part is needed only if ProfParser is enabled
# 
if ( NOT ENABLE_PROFPARSER )
   return ()
endif ()


set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


# 
# Extra pre-processing steps if any profiling is active
# 
set (PWD ${CMAKE_CURRENT_LIST_DIR}/ProfParser)

print (PWD)

add_custom_command (
   OUTPUT  ${PWD}/BLProfParser.tab.cpp ${PWD}/BLProfParser.tab.H 
   COMMAND cat   BLProfParser.y $(SED0) $(SED1) > BLProfParserNC.y
   COMMAND bison --defines=BLProfParser.tab.H --output=BLProfParser.tab.cpp BLProfParserNC.y
   COMMAND rm    BLProfParserNC.y
   DEPENDS ${PWD}/BLProfParser.y
   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ProfParser
   COMMENT "Generating BLProfParser.tab.H and BLProfParser.tab.cpp" )

add_custom_command (
   OUTPUT ${PWD}/BLProfParser.lex.yy.cpp
   COMMAND flex --outfile=BLProfParser.lex.yy.cpp BLProfParser.l
   DEPENDS ${PWD}/BLProfParser.l
   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ProfParser
   COMMENT "Generating BLProfParser.lex.yy.cpp" )


list ( APPEND ALLHEADERS AMReX_BLProfStats.H AMReX_BLProfUtilities.H )

list ( APPEND CXXSRC     AMReX_BLProfUtilities.cpp AMReX_BLWritePlotFile.cpp BLProfParser.tab.cpp )
list ( APPEND CXXSRC     BLProfParser.lex.yy.cpp AMReX_BLProfStats.cpp AMReX_CommProfStats.cpp )
list ( APPEND CXXSRC     AMReX_RegionsProfStats.cpp )

list ( APPEND F77SRC     AMReX_AVGDOWN_${DIM}D.F )

# 
# Accumulate sources
# 
set ( ALLSRC ${CXXSRC} ${F77SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/ProfParser )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/ProfParser )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
#
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/ProfParser)
