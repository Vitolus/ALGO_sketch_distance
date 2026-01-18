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
SCALE=0.00001 # Reduced scale to keep sketch size < 50KB
# Genome sizes for Bloom filter dimensioning
SIZE_HS=3200000000 # Homo sapiens ~3.2Gb (Largest)

# To compare Bloom Filters, they must have the same size and hash functions.
# We size all filters based on the largest genome in the set to avoid saturation.
COMMON_GENOME_SIZE=$SIZE_HS

echo "Creating sketch for Homo sapiens..."
./script_stream.sh -s homo_sapiens | ./dist --create-sketch --k $K --scale $SCALE --genome-size $COMMON_GENOME_SIZE hs.sketch
echo ""
echo "Creating sketch for Saccharomyces cerevisiae..."
./script_stream.sh -s saccharomyces_cerevisiae | ./dist --create-sketch --k $K --scale $SCALE --genome-size $COMMON_GENOME_SIZE sc.sketch
echo ""
echo "Creating sketch for Mus musculus..."
./script_stream.sh -s mus_musculus | ./dist --create-sketch --k $K --scale $SCALE --genome-size $COMMON_GENOME_SIZE mm.sketch
echo ""
echo "Creating sketch for Loxodonta africana..."
./script_stream.sh -s loxodonta_africana | ./dist --create-sketch --k $K --scale $SCALE --genome-size $COMMON_GENOME_SIZE la.sketch
echo ""
echo "Creating sketch for Drosophila melanogaster..."
./script_stream.sh -s drosophila_melanogaster | ./dist --create-sketch --k $K --scale $SCALE --genome-size $COMMON_GENOME_SIZE dm.sketch
echo ""
echo "--- All sketches created successfully. ---"
echo ""
echo ""
echo "--- Computing distance matrix and generating phylogenetic trees... ---"
./dist --distance hs.sketch mm.sketch la.sketch dm.sketch sc.sketch
echo ""
echo "--- Analysis complete. ---"
