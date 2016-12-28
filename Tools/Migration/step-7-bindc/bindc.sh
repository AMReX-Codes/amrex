#!/usr/bin/env bash

# subroutine probinit (*) --> subroutine amrex_probinit (*) bind(c)

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.f" -o -name "*.f90" -o -name "*.F" -o -name "*.F90" \) -exec grep -Iq . {} \; -exec sed -i 's/\(\(subroutine\|SUBROUTINE\)\s*\)\(probinit\|PROBINIT\)\s*\((.*)\)/\1amrex_probinit \4 bind(c)/g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.f" -o -name "*.f90" -o -name "*.F" -o -name "*.F90" \) -exec grep -Iq . {} \; -exec sed -i 's/\(end\|END\)\s*\(subroutine\|SUBROUTINE\)\s*\(probinit\|PROBINIT\)/end subroutine amrex_probinit/g' {} +

