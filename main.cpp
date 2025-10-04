#include <algorithm>
#include <complex>
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

void printUsage(const string& program_name){
    cerr << "Usage: \n"
              << program_name << " --create-sketch [options] <output.sketch>\n"
              << program_name << " --distance <G1.sketch> <G2.sketch>\n\n"
              << "Options for --create-sketch:\n"
              << "--k <int>    K-mer size (default: 21)\n"
              << "--scale <float>    Scaling factor for sketch size (default: 0.001)\n"
              << "--seed <int>    Seed for hashing (default: 1469598103934665603)\n";
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
vector<vector<double>> computeDistance(const vector<string>& sketch_files){
    const size_t n = sketch_files.size();
    // store pairs of <base_name, original_path>
    vector<std::pair<string, string>> sorted_sketches;
    sorted_sketches.reserve(n);
    for(const auto& path : sketch_files){
        sorted_sketches.emplace_back(extractBaseName(path), path);
    }
    // Sort the vector of pairs
    std::sort(sorted_sketches.begin(), sorted_sketches.end(), [](const auto& a, const auto& b){
        return a.first < b.first;
    });
    // compute and store the n x n distance matrix
    vector<vector<double>> distance_matrix(n, vector<double>(n, 0.0));
    for(size_t i = 0; i < n; ++i){
        distance_matrix[i][i] = 0.0;
        for(size_t j = i + 1; j < n; ++j){ // matrix is symmetric (compute only upper triangle)
            // use the original paths for the calculation
            const string& file1 = sorted_sketches[i].second;
            const string& file2 = sorted_sketches[j].second;
            const double dist = computeDistanceBetweenSketches(file1, file2);
            distance_matrix[i][j] = dist;
            distance_matrix[j][i] = dist; // assign symmetric value
        }
    }
    // stream the results to standard output in Phylip format
    cout << std::fixed << std::setprecision(6);
    // number of sequences
    cout << n << endl;
    for(size_t i = 0; i < n; ++i){
        // sequence name
        cout << sorted_sketches[i].first;
        // distances for that row
        for(size_t j = 0; j < n; ++j){
            cout << " " << distance_matrix[i][j];
        }
        cout << endl;
    }
    return distance_matrix;
}

int main(int argc, char* argv[]){
    // performance with large I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if(argc < 2){
        printUsage(argv[0]);
        return 1;
    }
    vector<vector<double>> distance_matrix;
    if(const string command = argv[1]; command == "--create-sketch"){
        if (argc < 3) {
            cerr << "Error: --create-sketch requires exactly one output file name.\n";
            printUsage(argv[0]);
            return 1;
        }
        unsigned k = 21;
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
        createSketch(output_file, k, scale, seed);
    }else if(command == "--distance"){
        if(argc < 4){
            cerr << "Error: --distance requires at least two sketch files.\n";
            printUsage(argv[0]);
            return 1;
        }
        vector<string> files;
        for(int i = 2; i < argc; ++i){
            files.emplace_back(argv[i]);
        }
        distance_matrix = computeDistance(files);
    }else{
        cerr << "Error: Unknown command '" << command << "'.\n";
        printUsage(argv[0]);
        return 1;
    }
    return 0;

}