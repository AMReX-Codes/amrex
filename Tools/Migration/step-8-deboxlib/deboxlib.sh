#!/usr/bin/env bash

# remove #include <AMReX_windstd.H>
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" \) -exec grep -Iq . {} \; \
    -exec sed -i '/#include\s*\(<\|\"\)\s*AMReX_winstd\.H\s*\(>\|\"\)/d' {} +

# AMReX_BoxLib.H --> AMReX.H
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" \) -exec grep -Iq . {} \; \
    -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)AMReX_BoxLib\.H\(\s*\(>\|\"\)\)/\1AMReX\.H\3/g' {} +

# buildInfo.H --> AMReX_buildInfo.H
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" \) -exec grep -Iq . {} \; \
    -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)buildInfo\.H\(\s*\(>\|\"\)\)/\1AMReX_buildInfo\.H\3/g' {} +

# buildInfoGetBoxlibDir --> buildInfoGetAMReXDir
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -name "*.cpp" \
    -exec sed -i 's/buildInfoGetBoxlibDir/buildInfoGetAMReXDir/g' {} +

# --boxlib_home --> --amrex_home
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -name "Make.*" \
    -exec sed -i 's/--boxlib_home/--amrex_home/g;s/buildInfo\.cpp/AMReX_buildInfo\.cpp/g;s/buildInfo\.H/AMReX_buildInfo\.H/g;s/buildInfo\.o/AMReX_buildInfo\.o/g' {} + \


