#!/bin/bash

# Specify the root directory for traversal
root_directory="./"

# Create a temporary file to store the list of directories
temp_file=$(mktemp)
find "$root_directory" -type d -print0 > "$temp_file"

# Initialize an empty string to store file information
files_list=""

# file counter to be used as time entry
count=0

ls -d */ | sort -z > temporary_file

# Read from the temporary file
for dir in */; do

    dir_name=$(basename "$dir")

    # Check if the folder starts with "plt" and contains a file named "Header"
    if [[ "$dir_name" == plt* && -f "$dir/Header" ]]; then
        # Extract version number from folder name
        version="${dir_name#plt}"
        echo $version

        # Create file information
        files_list+="$(printf "{ \"name\": \"plt$version\", \"time\": $count},")"
        files_list+=$'\n'

        ((count++))
    fi
done < "$temp_file"

# Remove trailing comma from the last entry
files_list="${files_list%,}"

# Create the final JSON structure
# Header line
header_line="{ \"file-series-version\": \"1.0\", \"files\": ["
# Write the files list
all_files="$(printf '%s\n' "$files_list") ] }"

file_series_data="$header_line"
file_series_data+=$'\n'
file_series_data+="$all_files"

# Write the generated JSON structure to a file named plot_files.series
echo "$file_series_data" > plot_files.series

# Remove the temporary file
rm "$temp_file"

echo "JSON structure has been written to plot_files.series"
