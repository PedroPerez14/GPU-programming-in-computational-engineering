#   Shallow Water Equations 2D Simulator (SWE2D)
#  ===========================================================================================
#   Square grid
#  ===========================================================================================
#   Version 1.0 - Mar 2025
#  ===========================================================================================
#  ===========================================================================================
#   Computational Hydraulics Group - University of Zaragoza   
#  =========================================================================================== 

#!/bin/bash

# Source files
SRC1="swe2d.c"
SRC2="lib/shallow_water.c"
OUTPUT="runSWE2D"
RESULTS_PATH="outputFiles"

# Remove old files
rm -f $OUTPUT
rm -rf $RESULTS_PATH

# Compile with nvcc (GPU compiler with C++ support)
nvcc -g -G -arch=sm_80 -o $OUTPUT $SRC1 $SRC2 -lm -x cu -std=c++17  -lcublas

# Create results directory
mkdir -p $RESULTS_PATH

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "GPU code compiled successfully. Run with ./$OUTPUT"
else
    echo "Compilation failed"
fi
