#!/usr/bin/env bash

# Copyright 2020 Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This file loops over all WarpX source files, uses git to get the list of
# contributors, and writes a copyright header at the beginning of each file
# with the list of contributors.
#
# To use it, execute from the WarpX directory.
# Warning: it'll delete and create file tmp.txt
# It uses the gnu-sed (sed on Linux, gsed on MacOS)

rm -r tmp.txt
set -e

# Loop over all source files
pattern="\.c$|\.cpp$|\.F90$|\.h$|\.H$|\.py$|"\
"\.sh$|\.tex$|\.txt$|\.yml$|"\
"CMakeLists\.txt|inputs"
for i in $(find . \
                -not -path "./.git/*" \
                -not -path "./.idea/*" \
                -not -path "./Docs/source/api/*" \
                -not -path "./Docs/build/*" \
                -not -path "./Docs/doxyxml/*" \
                -not -path "*wp_parse*" \
                -not -path "*LEGAL.txt" \
                -not -path "*LICENSE.txt" \
                -not -path "./tmp_build_dir/*" \
                -not -path "*/inputs*" \
                -not -path "*/PICMI_inputs*" \
                -not -path "./Tools/performance_tests/performance_log.txt" \
                -type f | \
               grep -P "${pattern}")
do
    echo "   ---   " $i
    # If Copyright information is present, remove it.
    # WARNING: This only works for C++ files.
    gsed -i '/^\/\* Copyright/,/\*\//{/^#/!{/^\$/!d}}' $i ; sleep 0.01
    # Get year of first modification of the file
    year_line=`git log --format=%aD $i | tail -1`
    year_first=($year_line)
    year_first=${year_first[3]}
    # Get year of last modification of the file
    year_line=`git log --format=%aD $i | head -1`
    year_last=($year_line)
    year_last=${year_last[3]}
    # Format year string, something like "2020" or "2016-2020"
    if [ $year_first == $year_last ]; then year_string=$year_first; else year_string=$year_first-$year_last; fi
    cp $i tmp.txt
    # If bash or python or txt or yml file, comment character is #
    if [ "${i: -2}" == "sh" ] || [ "${i: -2}" == "py" ] || [ "${i: -3}" == "txt" ] || [ "${i: -3}" == "yml" ]; then
        echo "sh or py or txt"
        pattern1="#"
        pattern2="#"
        pattern3="
"
    # If C or C++ file, comment characters are /*, * and */
    elif [ "${i: -1}" == "H" ] || [ "${i: -3}" == "cpp" ] || [ "${i: -1}" == "c" ] || [ "${i: -1}" == "h" ]; then
        echo "cpp"
        pattern1="/*"
        pattern2=" *"
        pattern3="
 */"
    # If Fortran file, comment character is !
    elif [ "${i: -3}" == "F90" ]; then
        echo "Fortran"
        pattern1="!"
        pattern2="!"
        pattern3="
"
    # If Latex, comment character is %
    elif [ "${i: -3}" == "tex" ]; then
        echo "tex"
        pattern1="%"
        pattern2="%"
        pattern3="
"
    else
        echo "error: unknown file type"
        exit
    fi
    # Get formatted authors list
    # sorted, unique, delete authors "Tools" (used by ax3l), and remove newlines
    authors_list=`git log --follow --pretty=format:'%aN' $i | sort | uniq | grep -v Tools | gsed 's/$/, /g'`
    # Put 2 authors per line, to avoid very long lines.
    authors_list=`echo $authors_list | gsed 's/\([^,]*,[^,]*,[^,]*\),/\1,\n'"$pattern2"'/g' | gsed s/,$//g`
    # Copy current file + Copyright to tmp.txt
    # rm -rf tmp.txt
    cp $i tmp.txt
    first_line=`head -n1 $i`
    # If a shebang is present, keep it as first line
    if [ "${first_line:0:2}" == "#!" ]; then
        echo "keeping shebang"
        echo "$first_line" > tmp.txt
        echo "" >> tmp.txt
    else
        truncate -s 0 tmp.txt
    fi
    # Write copyright
    echo "$pattern1 Copyright $year_string $authors_list
$pattern2
$pattern2 This file is part of WarpX.
$pattern2
$pattern2 License: BSD-3-Clause-LBNL$pattern3" >> tmp.txt
    # If no shebang, put first line after Copyright
    if [ "${first_line:0:2}" != "#!" ]; then
        echo "$first_line" >> tmp.txt
    fi
    # Then copy the content of the file
    tail -n +2 $i >> tmp.txt
    # Then overwrite current file with tmp.txt
    mv tmp.txt $i
done
