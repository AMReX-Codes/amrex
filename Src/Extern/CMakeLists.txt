set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

# 
# For now we accumulate sources for folder amrdata and ProfParser only
#
set (CXXSRC
   amrdata/AMReX_AmrData.cpp  amrdata/AMReX_DataServices.cpp
   ProfParser/AMReX_XYPlotDataList.cpp)

set (F77SRC amrdata/AMReX_FABUTIL_${BL_SPACEDIM}D.F)

set (ALLHEADERS
   amrdata/AMReX_AmrData.H  amrdata/AMReX_AmrvisConstants.H  amrdata/AMReX_DataServices.H
   ProfParser/AMReX_XYPlotDataList.H )


# Extra steps if any profiling is active
if ( AMREX_ANY_PROFILING )

   set (PWD ${CMAKE_CURRENT_LIST_DIR}/ProfParser)
   
   # add_custom_command (
   #    OUTPUT  BLProfParser.tab.H  BLProfParser.tab.cpp
   #    COMMAND cat ${PWD}/BLProfParser.y $(SED0) $(SED1) > BLProfParserNC.y
   #    COMMAND bison --defines=BLProfParser.tab.H --output=BLProfParser.tab.cpp BLProfParserNC.y
   #    COMMAND rm  BLProfParserNC.y
   #    DEPENDS ${PWD}/BLProfParser.y
   #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   #    COMMENT "Generating BLProfParser.tab.H and BLProfParser.tab.cpp" )


   # add_custom_command (
   #    OUTPUT BLProfParser.lex.yy.cpp
   #    COMMAND flex --outfile=BLProfParser.lex.yy.cpp ${PWD}/BLProfParser.l
   #    DEPENDS ${PWD}/BLProfParser.l
   #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   #    COMMENT "Generating BLProfParser.lex.yy.cpp" )

   
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
      
   list ( APPEND CXXSRC
      ProfParser/AMReX_BLProfUtilities.cpp ProfParser/AMReX_BLWritePlotFile.cpp
      ProfParser/BLProfParser.tab.cpp     ProfParser/BLProfParser.lex.yy.cpp
      ProfParser/AMReX_BLProfStats.cpp ProfParser/AMReX_CommProfStats.cpp
      ProfParser/AMReX_RegionsProfStats.cpp ProfParser/AMReX_XYPlotDataList.cpp )

   list ( APPEND F77SRC ProfParser/AMReX_AVGDOWN_${BL_SPACEDIM}D.F)

   list ( APPEND ALLHEADERS
      ProfParser/AMReX_BLProfStats.H  ProfParser/AMReX_BLProfUtilities.H
      ProfParser/AMReX_XYPlotDataList.H ProfParser/BLProfParser.tab.H)
   
endif ()





# Accumulate sources
set ( ALLSRC ${CXXSRC} ${F77SRC} )

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

