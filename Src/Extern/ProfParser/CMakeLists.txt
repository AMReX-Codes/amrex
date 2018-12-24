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

add_sources( AMReX_BLProfStats.H        AMReX_BLProfUtilities.H )
add_sources( AMReX_BLProfUtilities.cpp  AMReX_BLWritePlotFile.cpp )
add_sources( BLProfParser.tab.cpp )
add_sources( BLProfParser.lex.yy.cpp    AMReX_BLProfStats.cpp )
add_sources( AMReX_CommProfStats.cpp    AMReX_RegionsProfStats.cpp )
add_sources( AMReX_AVGDOWN_${DIM}D.F )
