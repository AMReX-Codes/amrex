project('amrex', 'cpp', 'fortran', 'c',
        default_options: ['c_std=c11', 'cpp_std=c++11'])

subdir('Src')
incdirs = include_directories(
    'Src/Amr',
    'Src/AmrCore',
    'Src/AmrTask/Amr',
    'Src/AmrTask/graph',
    'Src/AmrTask/rts_impls/Serial',
    'Src/Base',
    'Src/Boundary',
)


# libamrex = static_library('amrex',
#     [],
#     include_directories: incdirs,
#     fortran_args: [
#         '-DAMREX_SPACEDIM=3',
#         '-DBL_FORT_USE_LOWERCASE=1',
#         '-DBL_USE_MPI=0',
#         '-cpp',
#     ],
#     cpp_args: [
#         '-DAMREX_SPACEDIM=3',
#         '-DBL_FORT_USE_LOWERCASE=1',
#     ],
#     link_with: base,
#     install: true,
# )

libamrex_dep = declare_dependency(
    include_directories : incdirs_amrex,
    link_with : base,
)
