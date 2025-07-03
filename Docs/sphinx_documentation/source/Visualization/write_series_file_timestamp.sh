#!/bin/bash

# Specify the root directory for traversal
root_directory="./"

# Initialize an empty array to store directory information
declare -a dir_info

# Find all plt directories
for dir in plt*; do
    if [[ -d "$dir" && -f "$dir/Header" ]]; then
        # Extract time from the fifth line of the Header file
        time_value=$(sed -n '5p' "$dir/Header")

        # Store directory and time as a pair
        dir_info+=("$dir:$time_value")
    fi
done

# Sort directories numerically by their numbers
IFS=$'\n' sorted_dirs=($(for entry in "${dir_info[@]}"; do
    echo "$entry"
done | sort -t: -k1.4,1 -n))
unset IFS

# Initialize an empty string to store file information
files_list=""

# Process sorted directories
for entry in "${sorted_dirs[@]}"; do
    dir=$(echo "$entry" | cut -d: -f1)
    time_value=$(echo "$entry" | cut -d: -f2)

    # Round to 6 decimal places and format nicely
    rounded_time=$(awk -v t="$time_value" 'BEGIN{printf "%.6f", t}' | sed 's/\.0*$//;s/\.\([0-9]*[1-9]\)0*$/.\1/')

    echo "Processing $dir with time $rounded_time"

    # Create file information
    files_list+="$(printf "{ \"name\": \"$dir\", \"time\": $rounded_time},")"
    files_list+=$'\n'
done

# Remove trailing comma from the last entry
files_list="${files_list%,}"

# Create the final JSON structure
header_line="{ \"file-series-version\": \"1.0\", \"files\": ["
all_files="$(printf '%s\n' "$files_list") ] }"

file_series_data="$header_line"
file_series_data+=$'\n'
file_series_data+="$all_files"

# Write the generated JSON structure to a file named plot_files.series
echo "$file_series_data" > plot_files.series

echo "JSON structure has been written to plot_files.series"
