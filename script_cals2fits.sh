#!/bin/bash

source venv/bin/activate

# Set the directory where your files are located
input_dir="/run/media/jorge/TWOTERAS/DATA/phd/PSF_bulge"

# Create an array to store the list of files
files=("$input_dir"/*cals.gz)

# Iterate through all the files and run the cals_2_fits.py script
for ((i=0; i<${#files[@]}; i++)); do
  python cals_2_fits.py -f "${files[i]}"
done
