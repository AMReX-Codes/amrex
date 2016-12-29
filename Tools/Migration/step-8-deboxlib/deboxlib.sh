#!/usr/bin/env bash

# remove #include <AMReX_windstd.H>
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" \) -exec grep -Iq . {} \; -exec sed -i '/#include\s*\(<\|\"\)\s*AMReX_winstd\.H\s*\(>\|\"\)/d' {} +

# AMReX_BoxLib.H --> AMReX.H
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" \) -exec grep -Iq . {} \; -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)AMReX_BoxLib\.H\(\s*\(>\|\"\)\)/\1AMReX\.H\3/g' {} +
