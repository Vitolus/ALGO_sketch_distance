#include <iostream>
#include <string>
#include <vector>
#include <algorithm> // For std::sort
#include <utility>   // For std::pair
#include <iomanip>   // For std::fixed and std::setprecision

void printUsage(const std::string& program_name){
    std::cerr << "Usage: \n"
              << "  " << program_name << " --create-sketch <output.sketch>\n"
              << "  " << program_name << " --distance <G1.sketch> <G2.sketch> ...\n";
}

// helper function to get a clean sequence name from a file path
std::string extractBaseName(const std::string& path){
    // Find the last slash to remove the directory part
    size_t last_slash_pos = path.find_last_of("/\\");
    std::string filename = (last_slash_pos == std::string::npos) ? path : path.substr(last_slash_pos + 1);
    // remove file extension if present
    if(size_t extension_pos = filename.rfind(".sketch"); extension_pos != std::string::npos){
        return filename.substr(0, extension_pos);
    }
    return filename;
}

void createSketch(const std::string& output_file){
    // print progress messages to stderr
    std::cerr << "Starting sketch creation from standard input...\n";
    std::cerr << "   Output will be saved to: " << output_file << "\n";
    // TODO: initialize sketching data structure
    char current_base;
    unsigned long long base_count = 0;
    // reading one character at a time from stdin
    while(std::cin.get(current_base)){
        // upstream commands already filter for ACTG
        // TODO: add current_base to the sketch
        base_count++;
    }
    // finalize and save the sketch
    std::cerr << "   Processed " << base_count << " bases.\n";
    std::cerr << "   Sketch saved successfully to " << output_file << ".\n";
}

/**
 * @brief load two sketch files and compute their distance.
 */
double computeDistanceBetweenSketches(const std::string& sketch_file1, const std::string& sketch_file2){
    double dist = 0.0;
    // The distance from a sequence to itself is always 0.
    if(sketch_file1 == sketch_file2){
        return dist;
    }
    // TODO: add distance computation logic
    return dist;
}

/**
 * @brief Compute and print the n x n distance matrix in Phylip format.
 */
void computeDistance(const std::vector<std::string>& sketch_files){
    const size_t n = sketch_files.size();
    // store pairs of <base_name, original_path>
    std::vector<std::pair<std::string, std::string>> sorted_sketches;
    for(const auto& path : sketch_files){
        sorted_sketches.emplace_back(extractBaseName(path), path);
    }
    // Sort the vector of pairs
    std::sort(sorted_sketches.begin(), sorted_sketches.end());
    // compute and store the n x n distance matrix
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));
    for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
            // use the original paths for the calculation
            const std::string& file1 = sorted_sketches[i].second;
            const std::string& file2 = sorted_sketches[j].second;
            distance_matrix[i][j] = computeDistanceBetweenSketches(file1, file2);
        }
    }
    // stream the results to standard output in Phylip format
    std::cout << std::fixed << std::setprecision(6);
    // number of sequences
    std::cout << n << "\n";
    for(size_t i = 0; i < n; ++i){
        // sequence name
        std::cout << sorted_sketches[i].first;
        // distances for that row
        for(size_t j = 0; j < n; ++j){
            std::cout << " " << distance_matrix[i][j];
        }
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]){
    // good practice for performance with large I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if(argc < 2){
        printUsage(argv[0]);
        return 1;
    }
    if(std::string command = argv[1]; command == "--create-sketch"){
        if (argc != 3) {
            std::cerr << "Error: --create-sketch requires exactly one output file name.\n";
            printUsage(argv[0]);
            return 1;
        }
        createSketch(argv[2]);
    }else if(command == "--distance"){
        if(argc < 4){
            std::cerr << "Error: --distance requires at least two sketch files.\n";
            printUsage(argv[0]);
            return 1;
        }
        std::vector<std::string> files;
        for(int i = 2; i < argc; ++i){
            files.emplace_back(argv[i]);
        }
        computeDistance(files);
    }else{
        std::cerr << "Error: Unknown command '" << command << "'.\n";
        printUsage(argv[0]);
        return 1;
    }
    return 0;
}