macro(setup_ccache)
    find_program(AMReX_CCACHE_EXE ccache)
    if (AMReX_CCACHE_EXE)
        message(STATUS "Found CCache: ${AMReX_CCACHE_EXE}")
        set(CMAKE_CXX_COMPILER_LAUNCHER "${AMReX_CCACHE_EXE}")
        if (AMREX_GPU_BACKEND STREQUAL CUDA)
            set(CMAKE_CUDA_COMPILER_LAUNCHER "${AMReX_CCACHE_EXE}")
        endif()
	if (MSVC AND (CMAKE_GENERATOR MATCHES "Visual Studio"))
            # https://github.com/ccache/ccache/wiki/MS-Visual-Studio#usage-with-cmake
            file(COPY_FILE
                 ${AMReX_CCACHE_EXE} ${CMAKE_BINARY_DIR}/cl.exe
                 ONLY_IF_DIFFERENT)

            # By default Visual Studio generators will use /Zi which is not compatible
            # with ccache, so tell Visual Studio to use /Z7 instead.
            message(STATUS "Setting MSVC debug information format to 'Embedded'")
            set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")

            set(CMAKE_VS_GLOBALS
              "CLToolExe=cl.exe"
              "CLToolPath=${CMAKE_BINARY_DIR}"
              "TrackFileAccess=false"
              "UseMultiToolTask=true"
              "DebugInformationFormat=OldStyle"
            )
        endif()
    else()
        message(WARNING "AMReX_CCACHE was enabled, but ccache was not found.")
    endif()
    mark_as_advanced(AMReX_CCACHE_EXE)
endmacro(setup_ccache)
