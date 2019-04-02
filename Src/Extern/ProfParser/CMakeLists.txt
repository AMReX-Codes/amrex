target_compile_definitions(amrex PUBLIC $<BUILD_INTERFACE:AMREX_USE_PROFPARSER BL_USE_PROFPARSER>)
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# 
# Extra pre-processing steps if any profiling is active
# 
set (PWD ${CMAKE_CURRENT_LIST_DIR})

add_custom_command (
   OUTPUT  ${PWD}/BLProfParser.tab.cpp ${PWD}/BLProfParser.tab.H 
   COMMAND cat   BLProfParser.y $(SED0) $(SED1) > BLProfParserNC.y
   COMMAND bison --defines=BLProfParser.tab.H --output=BLProfParser.tab.cpp BLProfParserNC.y
   COMMAND rm    BLProfParserNC.y
   DEPENDS ${PWD}/BLProfParser.y
   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
   COMMENT "Generating BLProfParser.tab.H and BLProfParser.tab.cpp" )

add_custom_command (
   OUTPUT ${PWD}/BLProfParser.lex.yy.cpp
   COMMAND flex --outfile=BLProfParser.lex.yy.cpp BLProfParser.l
   DEPENDS ${PWD}/BLProfParser.l
   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
   COMMENT "Generating BLProfParser.lex.yy.cpp" )

unset (PWD)



target_sources( amrex
   PRIVATE 
   AMReX_BLProfStats.H
   AMReX_BLProfStats.cpp
   AMReX_BLProfUtilities.H
   AMReX_BLProfUtilities.cpp
   AMReX_BLWritePlotFile.cpp
   BLProfParser.tab.cpp
   BLProfParser.lex.yy.cpp
   AMReX_CommProfStats.cpp
   AMReX_RegionsProfStats.cpp 
   AMReX_AVGDOWN_${DIM}D.F
   )
