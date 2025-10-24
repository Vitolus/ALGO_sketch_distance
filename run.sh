#!/bin/sh

# This script automates the process of building the project, generating sketches
# for several species

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Building the project... ---"
# Assuming g++ is available. The -O3 flag enables optimizations.
g++ -std=c++17 -O3 main.cpp FracMinHash.h FracMinHash.cpp phylogenerator.h BloomFilter.cpp BloomFilter.h -o dist
echo "Build complete. Executable 'dist' created."
echo ""
echo ""

echo "--- Creating genome sketches... ---"
K=13
SCALE=0.01
N_BASES=10000000
echo "Creating sketch for Saccharomyces cerevisiae..."
./script_stream.sh -s saccharomyces_cerevisiae | ./dist --create-sketch --k $K --scale $SCALE sc.sketch
echo ""
echo "Creating sketch for Homo sapiens..."
./script_stream.sh -s homo_sapiens | ./dist --create-sketch --k $K --scale $SCALE hs.sketch
echo ""
echo "Creating sketch for Mus musculus..."
./script_stream.sh -s mus_musculus | ./dist --create-sketch --k $K --scale $SCALE mm.sketch
echo ""
echo "Creating sketch for Loxodonta africana..."
./script_stream.sh -s loxodonta_africana | ./dist --create-sketch --k $K --scale $SCALE la.sketch
echo ""
echo "Creating sketch for Drosophila melanogaster..."
./script_stream.sh -s drosophila_melanogaster | ./dist --create-sketch --k $K --scale $SCALE dm.sketch
echo ""
echo "--- All sketches created successfully. ---"
echo ""
echo ""

echo "--- Computing distance matrix and generating phylogenetic trees... ---"
./dist --distance hs.sketch mm.sketch la.sketch dm.sketch sc.sketch
echo ""

echo "--- Analysis complete. ---"
