macro(setup_clang_tidy)
   find_program(AMReX_CLANG_TIDY_EXE NAMES "clang-tidy-17" "clang-tidy-16"
	"clang-tidy-15" "clang-tidy-14" "clang-tidy-13" "clang-tidy-12" "clang-tidy")
   if (AMReX_CLANG_TIDY_EXE)
      set(_tmp "")

      execute_process(COMMAND ${AMReX_CLANG_TIDY_EXE} --version
                      OUTPUT_VARIABLE _tmp)
      if (_tmp MATCHES "LLVM version ([0-9\.]+)")
         message(STATUS "Found clang-tidy ${CMAKE_MATCH_1}")
         if ("${CMAKE_MATCH_1}" VERSION_GREATER_EQUAL 12.0.0)
            # Config file not supported in earlier versions
            set(AMReX_CLANG_TIDY_CONFIG_FILE_NAME ${PROJECT_SOURCE_DIR}/.clang-tidy)
         endif()
      endif()

      # Need --extra-arg to suppress warnings like clang-diagnostic-unknown-warning-option
      # when GCC is used.
      set(AMReX_CLANG_TIDY_COMMAND "${AMReX_CLANG_TIDY_EXE};--extra-arg=-Wno-unknown-warning-option")
      if (AMReX_CLANG_TIDY_CONFIG_FILE_NAME)
         set(AMReX_CLANG_TIDY_COMMAND "${AMReX_CLANG_TIDY_COMMAND}"
             "--config-file=${AMReX_CLANG_TIDY_CONFIG_FILE_NAME}")
      endif()
      if (AMReX_CLANG_TIDY_WERROR)
         set(AMReX_CLANG_TIDY_COMMAND "${AMReX_CLANG_TIDY_COMMAND}"
             "--warnings-as-errors=*")
      endif()

      foreach(D IN LISTS AMReX_SPACEDIM)
         set_target_properties(amrex_${D}d PROPERTIES CXX_CLANG_TIDY "${AMReX_CLANG_TIDY_COMMAND}")
      endforeach()

      unset(_tmp)
   else()
      message(WARNING "AMReX_CLANG_TIDY was enabled, but clang-tidy was not found.")
   endif()
endmacro(setup_clang_tidy)
