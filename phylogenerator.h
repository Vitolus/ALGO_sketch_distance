#ifndef ALGO_SKETCH_DISTANCE_PHYLOGENERATOR_H
#define ALGO_SKETCH_DISTANCE_PHYLOGENERATOR_H

#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::make_shared;
using std::unique_ptr;

/**
 * @brief a simple struct to represent a node in the phylogenetic tree
 */
struct Node{
    // using unique_ptr for automatic memory management of the tree. When a Node is destroyed, it automatically
    // destroys its children
    std::unique_ptr<Node> left = nullptr;
    std::unique_ptr<Node> right = nullptr;
    std::string name; // name of the sequence (for leaves)
    double branch_length = 0.0; // distance to parent
    int id = 0; // unique identifier for the node
    int cluster_size = 1; // for UPGMA: number of leaves in this cluster
    static inline int next_id = 0; // ensure each node gets a unique ID

    explicit Node(std::string n = "") : name(std::move(n)), id(next_id++){}
    [[nodiscard]] bool is_leaf() const{
        return left == nullptr && right == nullptr;
    }
};

// TODO: build also graphic view of the tree
/**
 * @brief a helper function to print the tree structure to the console
 */
inline void printTreeRecursive(const Node* node, const std::string& prefix, bool is_left){
    if(!node) return;

    cout << prefix;
    cout << (is_left ? "├──" : "└──" );

    // print the value of the node
    if (node->is_leaf()) {
        cout << node->name << " [branch_length:" << node->branch_length << "]" << endl;
    } else {
        cout << "Node " << node->id << " [branch_length:" << node->branch_length << "]" << endl;
    }

    // enter the next level of the tree
    // the right child is printed first to get a more natural tree layout in the console
    printTreeRecursive(node->right.get(), prefix + (is_left ? "│   " : "    "), false);
    printTreeRecursive(node->left.get(), prefix + (is_left ? "│   " : "    "), true);
}

/**
 * @brief visualize the tree in a text-based format
 */
inline void printTree(const Node* root){
    if(!root) return;
    cout << "\n--- Tree Visualization ---" << endl;
    if (root->is_leaf()) {
        cout << root->name << endl;
    } else {
        cout << "Root" << endl;
    }
    printTreeRecursive(root->right.get(), "", false);
    printTreeRecursive(root->left.get(), "", true);
    cout << "------------------------\n";
}

/**
 * @brief convert a tree to Newick format
 */
inline void toNewickRecursive(const Node* node, std::string& newick_string){
    if(!node) return;
    if(node->is_leaf()){
        newick_string += node->name + ":" + std::to_string(node->branch_length);
    }else{
        newick_string += "(";
        toNewickRecursive(node->left.get(), newick_string);
        newick_string += ",";
        toNewickRecursive(node->right.get(), newick_string);
        newick_string += ")";
        // the root node in many Newick formats doesn't have a branch length
        if(node->branch_length > 0.0){
            newick_string += ":" + std::to_string(node->branch_length);
        }
    }
}
inline std::string toNewick(const Node* root){
    std::string newick_string = "";
    toNewickRecursive(root, newick_string);
    newick_string += ";";
    return newick_string;
}

/**
 * @brief perform UPGMA clustering
 */
inline void buildUPGMATree(const std::vector<std::string>& names, const std::vector<std::vector<double>>& matrix){
    if(names.empty()) return;
    // initialize active clusters: one for each leaf node
    std::vector<std::unique_ptr<Node>> clusters;
    clusters.reserve(names.size());
    for(const auto& name : names){
        clusters.push_back(std::make_unique<Node>(name));
    }
    auto current_matrix = matrix;
    while (clusters.size() > 1) {
        // find the pair of clusters (i, j) with the minimum distance
        double min_dist = std::numeric_limits<double>::max();
        int min_i = -1;
        int min_j = -1;
        for(size_t i = 0; i < current_matrix.size(); i++){
            for(size_t j = i + 1; j < current_matrix.size(); j++){
                if(current_matrix[i][j] < min_dist){
                    min_dist = current_matrix[i][j];
                    min_i = static_cast<int>(i);
                    min_j = static_cast<int>(j);
                }
            }
        }
        // create a new parent node for clusters i and j
        auto new_node = std::make_unique<Node>();
        // branch lengths are half the distance between the clusters
        clusters[min_i]->branch_length = min_dist / 2.0;
        clusters[min_j]->branch_length = min_dist / 2.0;
        // transfer ownership of the child nodes to the new parent node
        new_node->left = std::move(clusters[min_i]);
        new_node->right = std::move(clusters[min_j]);
        new_node->cluster_size = new_node->left->cluster_size + new_node->right->cluster_size;
        // update the distance matrix
        std::vector<double> new_distances;
        for(size_t k = 0; k < clusters.size(); k++){
            if(k == min_i || k == min_j) continue;
            // weighted average distance for the new cluster
            double dist = (current_matrix[min_i][k] * new_node->left->cluster_size +
                           current_matrix[min_j][k] * new_node->right->cluster_size) /
                          (new_node->left->cluster_size + new_node->right->cluster_size);
            new_distances.push_back(dist);
        }
        // create the next smaller matrix
        std::vector<std::vector<double>> next_matrix(current_matrix.size() - 1,
            std::vector<double>(current_matrix.size() - 1));
        int r_new = 1;
        int c_new = 1;
        for(size_t r = 0; r < current_matrix.size(); r++){
            if(r == min_i || r == min_j) continue;
            c_new = r_new + 1;
            for(size_t c = r + 1; c < current_matrix.size(); c++){
                if(c == min_i || c == min_j) continue;
                next_matrix[r_new][c_new] = next_matrix[c_new][r_new] = current_matrix[r][c];
                c_new++;
            }
            r_new++;
        }
        for(size_t k = 0; k < new_distances.size(); k++){
             next_matrix[0][k+1] = next_matrix[k+1][0] = new_distances[k];
        }
        current_matrix = next_matrix;
        // update the list of active clusters. Replace the first merged cluster with the new node
        clusters[min_i] = std::move(new_node);
        // remove the second merged cluster
        clusters.erase(clusters.begin() + min_j);
    }
    cout << toNewick(clusters[0].get()) << endl;
    printTree(clusters[0].get());
}

/**
 * @brief perform Neighbor-Joining
 */
inline void buildNJTree(const std::vector<std::string>& names, const std::vector<std::vector<double>>& matrix){
    if(names.empty()) return;
    std::vector<std::unique_ptr<Node>> clusters;
    clusters.reserve(names.size());
    for(const auto& name : names){
        clusters.push_back(std::make_unique<Node>(name));
    }
    auto current_matrix = matrix;
    std::vector<int> active_indices(names.size());
    std::iota(active_indices.begin(), active_indices.end(), 0); // Fills with 0, 1, 2, ...
    while(clusters.size() > 2){
        const auto n = static_cast<int>(clusters.size());
        // calculate the 'u' values (net divergence) for each cluster
        std::vector<double> u(n, 0.0);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i != j) u[i] += current_matrix[i][j];
            }
        }
        // find the pair (i, j) that minimizes the Q-criterion
        double min_q = std::numeric_limits<double>::max();
        int min_i = -1;
        int min_j = -1;
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                if(const double q = (n - 2) * current_matrix[i][j] - u[i] - u[j]; q < min_q){
                    min_q = q;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        // create a new parent node
        auto new_node = std::make_unique<Node>();
        const double dist_ij = current_matrix[min_i][min_j];
        // 4. Calculate branch lengths from the new node to i and j
        clusters[min_i]->branch_length = (dist_ij / 2.0) + (u[min_i] - u[min_j]) / (2.0 * (n - 2));
        clusters[min_j]->branch_length = dist_ij - clusters[min_i]->branch_length;
        new_node->left = std::move(clusters[min_i]);
        new_node->right = std::move(clusters[min_j]);
        // update the distance matrix
        std::vector<double> new_distances;
        for(int k = 0; k < n; k++){
            if(k == min_i || k == min_j) continue;
            double dist = (current_matrix[min_i][k] + current_matrix[min_j][k] - dist_ij) / 2.0;
            new_distances.push_back(dist);
        }
        std::vector next_matrix(n - 1, std::vector<double>(n - 1));
        int r_new = 1;
        for(int r = 0; r < n; r++){
            if(r == min_i || r == min_j) continue;
            int c_new = r_new + 1;
            for(int c = r + 1; c < n; c++){
                if(c == min_i || c == min_j) continue;
                next_matrix[r_new][c_new] = next_matrix[c_new][r_new] = current_matrix[r][c];
                c_new++;
            }
            r_new++;
        }
        for(size_t k = 0; k < new_distances.size(); ++k){
             next_matrix[0][k+1] = next_matrix[k+1][0] = new_distances[k];
        }
        current_matrix = next_matrix;
        // update cluster list
        clusters[min_i] = std::move(new_node);
        clusters.erase(clusters.begin() + min_j);
    }
    // join the last two clusters
    if(clusters.size() == 2){
        clusters[0]->branch_length = clusters[1]->branch_length = current_matrix[0][1] / 2.0;
        auto root = std::make_unique<Node>();
        root->left = std::move(clusters[0]);
        root->right = std::move(clusters[1]);
        cout << toNewick(root.get()) << endl;
        printTree(root.get());
        return; // Exit after processing the final tree
    }
    // Handle cases with 0 or 1 initial items, or if the loop finishes with one cluster
    if(!clusters.empty()){
        cout << toNewick(clusters[0].get()) << endl;
        printTree(clusters[0].get());
    }
}

#endif //ALGO_SKETCH_DISTANCE_PHYLOGENERATOR_H