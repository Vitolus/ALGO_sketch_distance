#include <iostream>
#include <string>
#include <vector>
#include <algorithm> // For std::sort
#include <utility>   // For std::pair
#include <iomanip>   // For std::fixed and std::setprecision

void printUsage(const std::string& program_name) {
    std::cerr << "Usage: \n"
              << "  " << program_name << " --create-sketch <output.sketch>\n"
              << "  " << program_name << " --distance <G1.sketch> <G2.sketch> ...\n";
}

// --- Helper function to get a clean sequence name from a file path ---
// E.g., "data/human.sketch" -> "human"
std::string extractBaseName(const std::string& path) {
    // Find the last slash to remove the directory part
    size_t last_slash_pos = path.find_last_of("/\\");
    std::string filename = (last_slash_pos == std::string::npos) ? path : path.substr(last_slash_pos + 1);

    // Find the .sketch extension to remove it
    size_t extension_pos = filename.rfind(".sketch");
    if (extension_pos != std::string::npos) {
        return filename.substr(0, extension_pos);
    }
    return filename;
}

void createSketch(const std::string& output_file) {
    // Print progress messages to standard error (stderr)
    std::cerr << "Starting sketch creation from standard input...\n";
    std::cerr << "   Output will be saved to: " << output_file << "\n";

    // TODO: YOUR SKETCHING LOGIC GOES HERE
    // 1. Initialize your sketching data structures.
    //    Example: MinHash sketcher(k, sketch_size);

    char current_base;
    unsigned long long base_count = 0;

    // 2. Loop, reading one character at a time from standard input (stdin).
    //    std::cin.get() is efficient and reads until the input stream is closed.
    while (std::cin.get(current_base)) {
        // The upstream commands already filter for ACTG, so we can
        // directly process the base.
        // Example: sketcher.add(current_base);
        base_count++;
    }

    // 3. Once the stream ends, finalize and save the sketch.
    //    Example: sketcher.save_to_file(output_file);

    std::cerr << "   Processed " << base_count << " bases.\n";
    std::cerr << "   Sketch saved successfully to " << output_file << ".\n";
}

// --- Placeholder for your actual distance calculation logic ---
// This is where you would load two sketch files and compute their distance.
// For this example, it just returns a dummy value.
double computeDistanceBetweenSketches(const std::string& sketch_file1, const std::string& sketch_file2) {
    // The distance from a sequence to itself is always 0.
    if (sketch_file1 == sketch_file2) {
        return 0.0;
    }
    // TODO: REPLACE THIS WITH YOUR REAL DISTANCE LOGIC!
    // This dummy logic creates a predictable, non-zero value for demonstration.
    return (double)(sketch_file1.length() % 10 + sketch_file2.length() % 10) / 10.0 + 0.1;
}

void computeDistance(const std::vector<std::string>& sketch_files) {
    const size_t n = sketch_files.size();

    // 1. EXTRACT and SORT sequence names
    // We store pairs of <base_name, original_path> to sort by name
    // while keeping track of the full path needed for calculations.
    std::vector<std::pair<std::string, std::string>> sorted_sketches;
    for (const auto& path : sketch_files) {
        sorted_sketches.push_back({extractBaseName(path), path});
    }

    // Sort the vector of pairs. By default, it sorts based on the first
    // element of the pair (the base name), which is exactly what we need.
    std::sort(sorted_sketches.begin(), sorted_sketches.end());

    // 2. CALCULATE and STORE the n x n distance matrix
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Use the original full paths for the calculation
            const std::string& file1 = sorted_sketches[i].second;
            const std::string& file2 = sorted_sketches[j].second;
            distance_matrix[i][j] = computeDistanceBetweenSketches(file1, file2);
        }
    }

    // 3. STREAM the results to standard output in Phylip format
    std::cout << std::fixed << std::setprecision(6); // Set floating point precision

    // First line is the number of sequences
    std::cout << n << "\n";

    // Subsequent n lines are the matrix rows
    for (size_t i = 0; i < n; ++i) {
        // First element is the (sorted) sequence name
        std::cout << sorted_sketches[i].first;

        // Then, the distances for that row
        for (size_t j = 0; j < n; ++j) {
            std::cout << " " << distance_matrix[i][j];
        }
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]) {
    // This is a good practice for performance with large I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }
    std::string command = argv[1];
    if (command == "--create-sketch") {
        if (argc != 3) {
            std::cerr << "Error: --create-sketch requires exactly one output file name.\n";
            printUsage(argv[0]);
            return 1;
        }
        createSketch(argv[2]);
    } else if (command == "--distance") {
        if (argc < 4) {
            std::cerr << "Error: --distance requires at least two sketch files.\n";
            printUsage(argv[0]);
            return 1;
        }
        std::vector<std::string> files;
        for (int i = 2; i < argc; ++i) {
            files.push_back(argv[i]);
        }
        computeDistance(files);
    } else {
        std::cerr << "Error: Unknown command '" << command << "'.\n";
        printUsage(argv[0]);
        return 1;
    }
    return 0;
}