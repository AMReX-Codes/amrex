#!/usr/bin/env bash

# remove #inlcude<AMReX_PROB_AMR_F.H> from Fortran .F files
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -name "*.F" -exec grep -Iq . {} \; -exec sed -i '/#include\s*\(<\|\"\)\s*AMReX_PROB_AMR_F\.H/d' {} +


# subroutine probinit (*) --> subroutine amrex_probinit (*) bind(c)

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.f" -o -name "*.f90" -o -name "*.F" -o -name "*.F90" \) -exec grep -Iq . {} \; -exec sed -i 's/\(\(subroutine\|SUBROUTINE\)\s*\)\(probinit\|PROBINIT\|FORT_PROBINIT\)\s*\((.*)\)/\1amrex_probinit \4 bind(c)/g' {} +

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.f" -o -name "*.f90" -o -name "*.F" -o -name "*.F90" \) -exec grep -Iq . {} \; -exec sed -i 's/\(end\|END\)\s*\(subroutine\|SUBROUTINE\)\s*\(probinit\|PROBINIT\|FORT_PROBINIT\)/end subroutine amrex_probinit/g' {} +

