#
# This file defines the following variables
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
list (APPEND NERSC_MACHINES cori )

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
list ( APPEND ACLF_MACHINES mira theta )

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

# function to assign SITE
function ( assign_site group machine site_name )
   foreach ( item ${${group}} )
      string ( FIND "${${machine}}" "${item}" pos ) 
      if ( ${pos} GREATER -1 )
	 set ( SITE ${site_name} PARENT_SCOPE )
	 return ()
      endif ()
   endforeach()
endfunction ()

assign_site ( CCSE_MACHINES MACHINE ccse )

assign_site ( NERSC_MACHINES MACHINE nersc )

assign_site ( OLCF_MACHINES MACHINE olcf )

assign_site ( LLNL_MACHINES MACHINE llnl )

assign_site ( ALCF_MACHINES MACHINE alcf )

assign_site ( FLASH_MACHINES MACHINE flash )

assign_site ( NREAL_MACHINES MACHINE nrel )

message ( STATUS "Machine name is ${MACHINE}")
message ( STATUS "Site name is ${SITE}" )




