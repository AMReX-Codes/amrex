#!/usr/bin/env bash

#  #include <ARMeX_Array.H>   #include "AMReX_Array.H"
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)'AMReX_Array\.H'\(\s*\(>\|\"\)\)/#include <AMReX_Vector\.H>/g' {} +

#  amrex::Array
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/amrex::Array\(\s*\)</amrex::Vector\1</g' {} +

#  ^Array<
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/^Array\(\s*\)</Vector\1</g' {} +

#  _Array<
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/\(\s\+\)Array\(\s*\)</\1Vector\2</g' {} +

#  XArray<, where X is < or ( or " or , or { or //
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/\(<\|(\|"\|,\|{|\/\/\)Array\(\s*\)</\1Vector\2</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/GetArrOfArrOf/GetVecOfVecOf/g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" -o -name "*.tex" \) -exec grep -Iq . {} \; -exec sed -i 's/GetArrOf/GetVecOf/g' {} +

