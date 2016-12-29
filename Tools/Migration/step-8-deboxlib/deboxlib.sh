#!/usr/bin/env bash

find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.H" -o -name "*.h" -o -name "*.cpp" \) -exec grep -Iq . {} \; -exec sed -i '/#include\s*\(<\|\"\)\s*AMReX_winstd\.H\s*\(>\|\"\)/d' {} +

