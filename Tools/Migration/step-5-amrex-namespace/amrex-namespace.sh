#!/usr/bin/env bash

OLD="BoxLib::"
NEW="amrex::"
echo ${OLD}" --> "${NEW}
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD}"'/'"${NEW}"'/g' {} +

