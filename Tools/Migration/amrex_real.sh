#!/usr/bin/env bash

# amrex::Real --> amrex_real in all _f.H and _F.H files

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*_f.H" -o -name "*_F.H" \) -exec grep -Iq . {} \; -exec sed -i 's/amrex::Real/amrex_real/g' {} +
