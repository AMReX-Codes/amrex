set ( AMREX_TUTORIALS_DIR ${CMAKE_CURRENT_LIST_DIR} )

set ( AMREX_TUTORIALS_SUBDIRS Amr) # For now only Amr

prepend ( AMREX_TUTORIALS_SUBDIRS ${AMREX_TUTORIALS_DIR})

#
# Get all the amrex properties 
#
get_target_property ( AMREX_INCLUDE_PATHS amrex INCLUDE_DIRECTORIES )


#
# Create a target to build all tutorials
#
add_custom_target (tutorials)

#
# Here we define a macro to accumulate all the tutorial targets
#
set ( TUTORIAL_TARGETS )

macro (add_tutorial targets_list)
   add_dependencies (tutorials ${targets_list})
endmacro ()

#
# Loop over subdirs
# 
foreach (subdir ${AMREX_TUTORIALS_SUBDIRS})
   
   file ( GLOB_RECURSE  ALL_TUTORIALS "${subdir}/*CMakeLists.txt" ) 

   #
   # Loop over all the tutorials in this subdir
   # 
   foreach ( this_tutorial ${ALL_TUTORIALS})

      get_filename_component (dir ${this_tutorial} DIRECTORY )
      add_subdirectory (${dir})

   endforeach ()

   
endforeach ()
