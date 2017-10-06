#!/usr/bin/env bash

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)'AMReX_Array\.H'\(\s*\(>\|\"\)\)/#include <AMReX_Vector\.H>/g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/amrex::Array</amrex::Vector</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/\(\s\+\)Array\(\s*\)</\1Vector\2</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/<Array\(\s*\)</<Vector\1</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/^Array\(\s*\)</Vector\1</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/(Array\(\s*\)</(Vector\1</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/"Array\(\s*\)</"Vector\1</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/,Array\(\s*\)</,Vector\1</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/\/\/Array\(\s*\)</\/\/Vector\1</g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/GetArrOfArrOf/GetVecOfVecOf/g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/GetArrOf/GetVecOf/g' {} +


