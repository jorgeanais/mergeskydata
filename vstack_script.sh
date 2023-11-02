#!/bin/bash

# Set the directory where your files are located
input_dir="/home/jorge/Documents/data/virac2/merged"
output_file="/home/jorge/Documents/data/virac2/merged/all.fits"

# Create an array to store the list of files
files=("$input_dir"/*.fits)

# Check if there are at least 2 files to merge
if [ ${#files[@]} -lt 2 ]; then
  echo "Not enough files to merge."
  exit 1
fi

# Initialize the output file with the first two files
temp_file="/tmp/temp_table_$(date +%s%N).fits"
python merge_tables.py -f1 "${files[0]}" -f2 "${files[1]}" -o "$temp_file"

# Iterate through the remaining files and merge them with the output file
for ((i=2; i<${#files[@]}; i++)); do
  python merge_tables.py -f1 "${files[i]}" -f2 "$temp_file" -o "$temp_file"
done

# Rename the final combined file to the desired output name
mv "$temp_file" "$output_file"

echo "Merged all files into $output_file"
