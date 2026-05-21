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

# Nombres de los archivos fuente
SRC1="swe2d.c"
SRC2="lib/shallow_water.c"
OUTPUT="runSWE2D"
resultspath="outputFiles"

rm $OUTPUT
rm -r $resultspath

# Compilar los archivos con g++
g++ -o $OUTPUT $SRC1 $SRC2 -lm -std=c++17
mkdir $resultspath

# Verificar si la compilación fue exitosa
if [ $? -eq 0 ]; then
    echo "Code compiled succesfully. Run with ./$OUTPUT"
else
    echo "Compilation fail"
fi
