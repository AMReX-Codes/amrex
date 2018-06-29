#
# Take care of subdir amrdata ( this can be enabled separatly or
# implicitly by enabling profparser )
# 

add_sources ( amrdata/AMReX_AmrData.H        amrdata/AMReX_XYPlotDataList.H )
add_sources ( amrdata/AMReX_DataServices.H   amrdata/AMReX_AmrvisConstants.H )
add_sources ( amrdata/AMReX_AmrData.cpp      amrdata/AMReX_XYPlotDataList.cpp )
add_sources ( amrdata/AMReX_DataServices.cpp amrdata/AMReX_FABUTIL_${DIM}D.F )

if (ENABLE_PROPARSER)
   add_sources ( amrdata/AMReX_XYPlotDataList.cpp amrdata/AMReX_XYPlotDataList.H )
endif ()


# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#

#
# The next part is needed only if ProfParser is enabled
# 
if ( NOT ENABLE_PROFPARSER )
   return ()
endif ()

# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#
# ----------------------------------------------------------------------#


# 
# Extra pre-processing steps if any profiling is active
# 
set (PWD ${CMAKE_CURRENT_LIST_DIR}/ProfParser)

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

unset (PWD)
add_sources ( ProfParser/AMReX_BLProfStats.H        ProfParser/AMReX_BLProfUtilities.H )
add_sources ( ProfParser/AMReX_BLProfUtilities.cpp  ProfParser/AMReX_BLWritePlotFile.cpp )
add_sources ( ProfParser/BLProfParser.tab.cpp )
add_sources ( ProfParser/BLProfParser.lex.yy.cpp    ProfParser/AMReX_BLProfStats.cpp )
add_sources ( ProfParser/AMReX_CommProfStats.cpp    ProfParser/AMReX_RegionsProfStats.cpp )
add_sources ( ProfParser/AMReX_AVGDOWN_${DIM}D.F )
