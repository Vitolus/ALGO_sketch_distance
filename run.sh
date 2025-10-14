#!/bin/sh

# This script automates the process of building the project, generating sketches
# for several species, calculating their distances, and running the unit test.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Building the project... ---"
# Assuming g++ is available. The -O3 flag enables optimizations.
g++ -std=c++17 -O3 main.cpp FracMinHash.h FracMinHash.cpp phylogenerator.h -o dist
echo "Build complete. Executable 'dist' created."
echo ""
echo ""

echo "--- Creating genome sketches (first 1MB of each)... ---"
echo "Creating sketch for Saccharomyces cerevisiae..."
./script_stream.sh -s saccharomyces_cerevisiae -l 1000000 | ./dist --create-sketch --k 21 --scale 0.05 sc.sketch
echo ""
echo "Creating sketch for Homo sapiens..."
./script_stream.sh -s homo_sapiens -l 1000000 | ./dist --create-sketch --k 21 --scale 0.05 hs.sketch
echo ""
echo "Creating sketch for Mus musculus..."
./script_stream.sh -s mus_musculus -l 1000000 | ./dist --create-sketch --k 21 --scale 0.05 mm.sketch
echo ""
echo "Creating sketch for Loxodonta africana..."
./script_stream.sh -s loxodonta_africana -l 1000000 | ./dist --create-sketch --k 21 --scale 0.05 la.sketch
echo ""
echo "Creating sketch for Drosophila melanogaster..."
./script_stream.sh -s drosophila_melanogaster -l 1000000 | ./dist --create-sketch --k 21 --scale 0.05 dm.sketch
echo ""
echo "--- All sketches created successfully. ---"
echo ""
echo ""

echo "--- Computing distance matrix and generating phylogenetic trees... ---"
./dist --distance hs.sketch mm.sketch la.sketch dm.sketch sc.sketch
echo ""

echo "--- Analysis complete. ---"
