#!/usr/bin/env bash

# Starting from the current directory, this script searches and
# replaces "BOXLIB_HOME" with "AMREX_HOME" and "/BoxLib" with "/amrex"
# in files named `Make.*`, `GNUmakefile`, and `FParallelMG.mak`.

OLD="BOXLIB_HOME"
NEW="AMREX_HOME"
echo ${OLD}" --> "${NEW}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "Make\.*" -o -name "GNUmakefile" -o -name "FParallelMG\.mak" \) -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD}"'/'"${NEW}"'/g' {} +

OLD="\/BoxLib"
NEW="\/amrex"
echo ${OLD}" --> "${NEW}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "Make\.package" -o -name "GNUmakefile" \) -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD}"'/'"${NEW}"'/g' {} +

#OLD="BOXLIB_HOME"
#NEW="AMREX_HOME"
#echo ${OLD}" --> "${NEW}
#find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "GPackage\.mak" -o -name "GMakedefs\.mak" -o -name "GMakerules\.mak" \) -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD}"'/'"${NEW}"'/g' {} +
