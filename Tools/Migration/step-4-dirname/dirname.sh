#!/usr/bin/env bash

echo "Remove 'FCOMP = ' line from GNUmakefile"
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -name "GNUmakefile" -exec sed -i '/FCOMP\s*=/d' {} +

OLD_SRC_PATH="Tools\/C_mk"
NEW_SRC_PATH="Tools\/GNUMake"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_AMRLib"
NEW_SRC_PATH="Amr"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_AmrCoreLib"
NEW_SRC_PATH="AmrCore"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_BaseLib"
NEW_SRC_PATH="Base"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_BoundaryLib"
NEW_SRC_PATH="Boundary"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_ParticleLib"
NEW_SRC_PATH="Particle"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +
