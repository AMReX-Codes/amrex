#
# This file provides the following variables
# 
#   CCSE_MACHINES
#   NERSC_MACHINES
#   OCLF_MACHINES
#   LLNL_MACHINES
#   ALCF_MACHINES
#   NREAL_MACHINES
#   FLASH_MACHINES
#   SITE
#   MACHINE


#
#  CCSE machines
#
set ( CCSE_MACHINES )
list ( APPEND CCSE_MACHINES angilas atragon baragon battra biollante ebirah gamera gigan )
list ( APPEND CCSE_MACHINES gimantis godzilla gojira hedorah kiryu kumonga manda )
list ( APPEND CCSE_MACHINES megalon mothra rodan varan naphta orga ghidorah )

#
# NERSC machines
# 
set ( NERSC_MACHINES )
list (APPEND NERSC_MACHINES edison cori )

#
# OLCF machines
# 
set ( OLCF_MACHINES )
list ( APPEND OLCF_MACHINES titan summit summitdev )

#
# LLNL machines
# 
set ( LLNL_MACHINES )
list ( APPEND LLNL_MACHINES ray rzmanta )

#
# ALCF machines
#
set ( ACLF_MACHINES )
list ( APPEND ACLF_MACHINES mmira theta )

#
# NREL machines
#
set ( NREL_MACHINES )
list ( APPEND NREL_MACHINES merlin anymachine )

#
# FLASH machines
#
set ( FLASH_MACHINES )
list ( APPEND FLASH_MACHINES asterix )

#
# Find machine and site where AMReX is being used
#
set ( SITE     unknown )
set ( MACHINE  unknown )
cmake_host_system_information ( RESULT MACHINE QUERY HOSTNAME )

list ( FIND CCSE_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE ccse )
endif ()

list ( FIND NERSC_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE nersc )
endif ()

list ( FIND OLCF_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE olcf )
endif ()

list ( FIND LLNL_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE llnl )
endif ()

list ( FIND ALCF_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE alcf )
endif ()

list ( FIND FLASH_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE flash )
endif ()

list ( FIND NREL_MACHINES ${MACHINE} index)
if ( ${index} GREATER -1 )
   set ( SITE nrel )
endif ()

message ( STATUS "Machine name is ${MACHINE}")
message ( STATUS "Site name is ${SITE}" )




