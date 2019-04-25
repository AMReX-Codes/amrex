set ( AMREX_TUTORIALS_DIR ${CMAKE_CURRENT_LIST_DIR} )


set ( AMREX_TUTORIALS_SUBDIRS Amr Basic CVODE GPU LinearSolvers ) # For now only Amr

prepend ( AMREX_TUTORIALS_SUBDIRS ${AMREX_TUTORIALS_DIR})

#
# Create a target to build all tutorials and an indentical one for
# test_install
#
# !!!!!! ATTENTION: target "tutorial" and "test_install" are the same for
#  the time being until we decide what to do with test_install
# 
add_custom_target (tutorials)
add_custom_target (test_install)

#
# Here we define a macro to accumulate all the tutorial targets
#
set ( TUTORIAL_TARGETS )

macro (add_tutorial targets_list)
   add_dependencies (tutorials ${targets_list})
   add_dependencies (test_install ${targets_list})
endmacro ()

#
# Loop over subdirs
# 
foreach (subdir ${AMREX_TUTORIALS_SUBDIRS})
   
   file ( GLOB_RECURSE  ALL_TUTORIALS "${subdir}/*CMakeLists.txt" ) 

   foreach ( this_tutorial ${ALL_TUTORIALS})
      get_filename_component (dir ${this_tutorial} DIRECTORY )
      add_subdirectory (${dir})
   endforeach ()
   
endforeach ()


# Just write something if
add_custom_command(TARGET test_install POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan
                   "AMReX installation is working")
