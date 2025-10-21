#include <algorithm>
#include <complex>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "FracMinHash.h"
#include "phylogenerator.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

void printUsage(const std::string& program_name) {
    std::cerr << "Usage: \n"
              << "  " << program_name << " --create-sketch [options] <output.sketch>\n"
              << "  " << program_name << " --distance <G1.sketch> <G2.sketch> ...\n\n"
              << "Commands:\n"
              << "  --create-sketch  Creates a sketch from a genome sequence provided via stdin.\n"
              << "  --distance       Calculates and prints the distance matrix, UPGMA tree, and NJ tree.\n\n"
              << "Options for --create-sketch:\n"
              << "  --k <int>        K-mer size. Overrides --d-max if both are provided. (default: 21)\n"
              << "  --d-max <float>  Max evolutionary distance. Used to calculate k if --k is not set.\n"
              << "  --scale <float>  Scaling factor (default: 0.001)\n"
              << "  --seed <int>     Seed for hashing (default: 1469598103934665603)\n";
}

/**
 * @brief helper function to get a sequence name from a file path
 */
string extractBaseName(const string& path){
    // Find the last slash to remove the directory part
    const size_t last_slash_pos = path.find_last_of("/\\");
    string filename = (last_slash_pos == string::npos) ? path : path.substr(last_slash_pos + 1);
    // remove file extension if present
    if(const size_t extension_pos = filename.rfind(".sketch"); extension_pos != string::npos){
        return filename.substr(0, extension_pos);
    }
    return filename;
}

void createSketch(const string& output_file, const unsigned k, const double scale, const uint64_t seed){
    // print progress messages to stderr
    cerr << "Starting sketch creation from standard input...\n";
    cerr << "   Parameters: k=" << k << ", scale=" << scale << ", seed=" << seed << endl;
    cerr << "   Output will be saved to: " << output_file << endl;
    FracMinHash sketch(scale, k, seed);
    char current_base;
    unsigned long long base_count = 0;
    // reading one character at a time from stdin
    while(std::cin.get(current_base)){
        // upstream commands already filter for ACTG
        sketch.add_char(current_base);
        base_count++;
    }
    // finalize and save the sketch
    try{
        sketch.save(output_file);
    }catch(const std::exception& e){
        cerr << "Error saving sketch: " << e.what() << endl;
        return;
    }
    cerr << "   Processed " << base_count << " bases\n";
    cerr << "   Retained " << sketch.sketch_size() << " hashes in the sketch\n";
    cerr << "   Sketch saved successfully to " << output_file << endl;
    cout << "Memory usage: " << (std::filesystem::file_size(output_file) / 1024) << " bytes\n";
}

/**
 * @brief load two sketch files and compute their distance.
 */
double computeDistanceBetweenSketches(const string& sketch_file1, const string& sketch_file2){
    double dist = 0.0;
    if(sketch_file1 == sketch_file2){
        return dist; // distance to self is always 0
    }
    try{
        const FracMinHash sketch1 = FracMinHash::load(sketch_file1);
        const FracMinHash sketch2 = FracMinHash::load(sketch_file2);
        dist = sketch1.distance(sketch2);
    }catch(const std::exception& e){
        cerr << "Error reading sketches: " << e.what() << endl;
        return -1.0; // indicate error with a negative distance
    }
    return dist;
}

/**
 * @brief Compute and print the n x n distance matrix in Phylip format.
 */
std::pair<std::vector<std::string>, std::vector<std::vector<double>>> computeDistance(const vector<string>& sketch_files){
    const size_t n = sketch_files.size();
    // store pairs of <base_name, original_path>
    vector<std::pair<string, string>> sorted_sketches;
    sorted_sketches.reserve(n);
    for(const auto& path : sketch_files){
        sorted_sketches.emplace_back(extractBaseName(path), path);
    }
    // Sort the vector of pairs
    std::sort(sorted_sketches.begin(), sorted_sketches.end());
    // compute and store the n x n distance matrix
    std::vector<std::string> names(n);
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    for(size_t i = 0; i < n; i++){
        names[i] = sorted_sketches[i].first;
        for(size_t j = i + 1; j < n; j++){ // matrix is symmetric (compute only upper triangle)
            // use the original paths for the calculation
            const double dist = computeDistanceBetweenSketches(sorted_sketches[i].second, sorted_sketches[j].second);
            matrix[i][j] = matrix[j][i] = dist;
        }
    }
    return {names, matrix};
}

/**
 * @brief print distance matrix in Phylip format
 */
void printDistanceMatrix(const std::vector<std::string>& names, const std::vector<std::vector<double>>& matrix) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << names.size() << "\n";
    for (size_t i = 0; i < names.size(); ++i) {
        std::cout << names[i];
        for (size_t j = 0; j < names.size(); ++j) {
            std::cout << " " << matrix[i][j];
        }
        std::cout << "\n";
    }
}

void unitTest(){
    // create three small sketches from hardcoded sequences
    cout << "--- Running Unit Test ---\n";
    FracMinHash sketch1(0.1, 5);
    std::string seq1_str = "CTACTACGCCGATTCTGCTG";
    for(char c : seq1_str) sketch1.add_char(c);
    FracMinHash sketch2(0.1, 5);
    std::string seq2_str = "CTACTACGCCAATTCTGCTG";
    for(char c : seq2_str) sketch2.add_char(c);
    FracMinHash sketch3(0.1, 5);
    std::string seq3_str = "ATACTACGCCGATTCTGCTG";
    for(char c : seq3_str) sketch3.add_char(c);
    // compute distances
    const double dist12 = sketch1.distance(sketch2);
    const double dist13 = sketch1.distance(sketch3);
    const double dist23 = sketch2.distance(sketch3);
    cout << "Distance seq1-seq2: " << dist12 << " (Expected ~0.4762)\n";
    cout << "Distance seq1-seq3: " << dist13 << " (Expected ~0.1176)\n";
    cout << "Distance seq2-seq3: " << dist23 << " (Expected ~0.5455)\n\n";
    // prepare data for tree building
    std::vector<std::string> names = {"seq1", "seq2", "seq3"};
    std::vector<std::vector<double>> matrix = {
        {0.0, dist12, dist13},
        {dist12, 0.0, dist23},
        {dist13, dist23, 0.0}
    };
    cout << "--- Phylip Distance Matrix ---\n";
    printDistanceMatrix(names, matrix);
    cout << "\n--- UPGMA Tree (Newick Format) ---\n";
    buildUPGMATree(names, matrix);
    cout << "\n--- Neighbor-Joining Tree (Newick Format) ---\n";
    buildNJTree(names, matrix);
    cout << "--- Unit Test Finished ---\n";
}

int main(int argc, char* argv[]){
    // performance with large I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    if(argc < 2){
        printUsage(argv[0]);
        return 1;
    }
    if(const string command = argv[1]; command == "--create-sketch"){
        if (argc < 3) {
            cerr << "Error: --create-sketch requires exactly one output file name.\n";
            printUsage(argv[0]);
            return 1;
        }
        unsigned k = 0;
        double d_max = 0.0;
        double scale = 0.001;
        uint64_t seed = 1469598103934665603ULL;
        string output_file;
        for(int i = 2; i < argc; i++){
            if(string arg = argv[i]; arg == "--k"){
                if(i + 1 >= argc){
                    cerr << "Error: --k requires an integer argument.\n";
                    printUsage(argv[0]);
                    return 1;
                }
                try{
                    k = std::stoi(argv[++i]);
                }catch(const std::exception& e){
                    cerr << "Error: --k argument must be an integer: " << e.what() << endl;
                }
            }else if(arg == "--d-max"){
                if(i + 1 >= argc){
                    cerr << "Error: --d-max requires a float argument.\n";
                    printUsage(argv[0]);
                    return 1;
                }
                try{
                    d_max = std::stod(argv[++i]);
                }catch(const std::exception& e){
                    cerr << "Error: --d-max argument must be a float: " << e.what() << endl;
                }
            }else if(arg == "--scale"){
                if(i + 1 >= argc){
                    cerr << "Error: --scale requires a double argument.\n";
                    printUsage(argv[0]);
                    return 1;
                }
                try{
                    scale = std::stod(argv[++i]);
                }catch(const std::exception& e){
                    cerr << "Error: --scale argument must be a double: " << e.what() << endl;
                }
            }else if(arg == "--seed"){
                if(i + 1 >= argc){
                    cerr << "Error: --seed requires an long argument.\n";
                    printUsage(argv[0]);
                    return 1;
                }
                try{
                    seed = std::stoll(argv[++i]);
                }catch(const std::exception& e){
                    cerr << "Error: --seed argument must be an integer: " << e.what() << endl;
                }
            }else{
                // assume it's the output file name
                if(!output_file.empty()){
                    cerr << "Error: only one output file name is allowed.\n";
                    printUsage(argv[0]);
                    return 1;
                }
                output_file = argv[i];
            }
        }
        if(output_file.empty()){
            cerr << "Error: output file name is required.\n";
            printUsage(argv[0]);
            return 1;
        }
        // Determine final k value based on user input
        if (k > 0) { // --k was provided, it takes precedence
            if (d_max > 0.0) {
                cerr << "Warning: Both --k and --d-max were provided. Using --k " << k << ".\n";
            }
        } else if (d_max > 0.0) { // --d-max was provided, but --k was not
            if (d_max <= 0.0 || d_max >= 1.0) {
                cerr << "Error: --d-max must be between 0.0 and 1.0 (exclusive).\n";
                return 1;
            }
            // Calculate k using the formula k = log(P_chance) / log(1 - d_max)
            const double p_chance = 0.01;
            k = static_cast<unsigned>(std::round(std::log(p_chance) / std::log(1.0 - d_max)));
            cerr << "Using automatically determined k=" << k << " from d_max=" << d_max << endl;
        } else { // Neither was provided, use default k
            k = 21;
        }
        // Clamp k to the valid range [1, 31] supported by FracMinHash
        k = std::max(1u, std::min(31u, k));
        auto start_time = std::chrono::high_resolution_clock::now();
        createSketch(output_file, k, scale, seed);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        cout << "Total time: " << elapsed.count() << "s\n";
    }else if(command == "--distance"){
        if(argc < 4){
            cerr << "Error: --distance requires at least two sketch files.\n";
            printUsage(argv[0]);
            return 1;
        }
        vector<string> files;
        for(int i = 2; i < argc; i++){
            files.emplace_back(argv[i]);
        }
        auto start_time = std::chrono::high_resolution_clock::now();
        auto [names, matrix] = computeDistance(files);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        cout << "Distance matrix computed in " << elapsed.count() << "s\n";
        cout << "--- Phylip Distance Matrix ---\n";
        printDistanceMatrix(names, matrix);
        cout << "\n--- UPGMA Tree (Newick Format) ---\n";
        buildUPGMATree(names, matrix);
        cout << "\n--- Neighbor-Joining Tree (Newick Format) ---\n";
        buildNJTree(names, matrix);
    }else if(command == "--test"){
        unitTest();
    }else{
        cerr << "Error: Unknown command '" << command << "'.\n";
        printUsage(argv[0]);
        return 1;
    }
    return 0;
}