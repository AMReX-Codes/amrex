#!/usr/bin/env bash

# Copyright 2020 Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This script updates release version number in all source files.
#
# updates occurences like "version = '??.??'" where ??.?? is the version number
# updates occurences like "WarpX v??.??" where ??.?? is the version number
#
# Requirements:
# - gnu grep (default grep on MacOS does not have same -P option)
# - gnu sed  (default grep on MacOS does not have same -i option)

set -e

# Get old tag number from git
git_tag=`git describe --tags`
old_release_number="${git_tag%%-*}"

# Construct new tag number from current date
d=`date +%Y.%m`
new_release_number=${d:2}

echo "Replace release number $old_release_number by $new_release_number"

# Loop over files and sed patterns with version number
pattern="\.c$|\.cpp$|\.F90$|\.h$|\.H$|\.ini$|\.md$|\.py$|"\
"\.rst$|\.sh$|\.tex$|\.txt$|\.xml$|\.yml$|"\
"CMakeLists\.txt|inputs"
for i in $(find .. \
                -not -path "../.git/*"   \
                -not -path "../.idea/*"  \
                -not -path "../Docs/source/api/*" \
                -not -path "../Docs/build/*" \
                -not -path "../Docs/doxyxml/*" \
                -not -path "*wp_parse*" \
                -not -path "../tmp_build_dir/*" \
                -type f | \
           grep -P "${pattern}")
do
    echo $i
    # update occurences like "version = '??.??'" where ??.?? is the version number
    # Note: sleep is required due to a bug in sed: without, the file
    # permissions are modified
    sed -i "s/version = "\'"$old_release_number"\'"/version = "\'"$new_release_number"\'"/g" $i && sleep .001
    # update occurences like "WarpX v??.??" where ??.?? is the version number
    sed -i "s/WarpX v$old_release_number/WarpX v$new_release_number/g" $i && sleep .001
done

echo ""
echo "WARNING: Remaining occurences of $old_release_number are listed below."
echo "         Is this expected? Or should these be updated too?"
echo ""
git grep "$old_release_number" .
