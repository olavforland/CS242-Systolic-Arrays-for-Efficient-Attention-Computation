#!/bin/bash

# Define arrays of K values and magnitude scales
K_VALUES=(2 4 8 16 32)
MAGNITUDE_SCALES=(0.01 0.1 1.0 10.0 100.0)

# Create results directory if it doesn't exist
mkdir -p results

# Remove previous results file
rm -f results/error_results.csv

# Loop over all combinations
for K in "${K_VALUES[@]}"; do
    for SCALE in "${MAGNITUDE_SCALES[@]}"; do
        echo "Running test with K=$K, Magnitude Scale=$SCALE"
        make clean
        make all K=$K MAGNITUDE_SCALE=$SCALE
        make sim
    done
done
