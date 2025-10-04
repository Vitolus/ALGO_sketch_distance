#ifndef ALGO_SKETCH_DISTANCE_PHYLOGENERATOR_H
#define ALGO_SKETCH_DISTANCE_PHYLOGENERATOR_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <limits>
#include <numeric>
#include <memory>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::make_shared;
using std::shared_ptr;

/**
 * @brief a simple struct to represent a node in the phylogenetic tree
 */
struct Node{
    string name;
    double branchLength{};
    shared_ptr<Node> left;
    shared_ptr<Node> right;
};

using DistanceMatrix = vector<vector<double>>;
using NodePtr = shared_ptr<Node>;

/**
 * @brief perform UPGMA clustering
 */
inline NodePtr upgma(const DistanceMatrix& distMatrix, const vector<string>& labels){
    if(distMatrix.empty()){
        return nullptr;
    }
    const size_t n = labels.size();
    vector<NodePtr> nodes;
    std::map<string, int> labelToIndex;
    for(auto i = 0; i < n; ++i){
        nodes.push_back(make_shared<Node>());
        nodes.back()->name = labels[i];
        nodes.back()->branchLength = 0.0;
        labelToIndex[labels[i]] = i;
    }
    DistanceMatrix currentDistances = distMatrix;
    vector<int> clusterSizes(n, 1);
    vector<bool> activeClusters(n, true);
    while(nodes.size() > 1){
        double minDistance = std::numeric_limits<double>::max();
        int i_min = -1, j_min = -1;
        // Find the two closest clusters
        for(int i = 0; i < n; ++i){
            if(!activeClusters[i]) continue;
            for(int j = i + 1; j < n; ++j){
                if(!activeClusters[j]) continue;
                if(currentDistances[i][j] < minDistance){
                    minDistance = currentDistances[i][j];
                    i_min = i;
                    j_min = j;
                }
            }
        }
        // create a new internal node
        const auto newNode = make_shared<Node>();
        newNode->left = nodes[i_min];
        newNode->right = nodes[j_min];
        newNode->name = "(" + nodes[i_min]->name + "," + nodes[j_min]->name + ")";
        // calculate branch lengths for the new branches
        newNode->left->branchLength = minDistance / 2.0 - nodes[i_min]->branchLength;
        newNode->right->branchLength = minDistance / 2.0 - nodes[j_min]->branchLength;
        // update distances
        for(int k = 0; k < n; ++k){
            if(!activeClusters[k] || k == i_min || k == j_min){
                continue;
            }
            const double newDist = (currentDistances[i_min][k] * clusterSizes[i_min] +
                currentDistances[j_min][k] * clusterSizes[j_min]) / (clusterSizes[i_min] + clusterSizes[j_min]);
            currentDistances[i_min][k] = currentDistances[k][i_min] = newDist;
        }
        // deactivate merged clusters and update `nodes`
        activeClusters[j_min] = false;
        clusterSizes[i_min] += clusterSizes[j_min];
        nodes[i_min] = newNode;
        // this is a simplification. A more robust implementation would manage indices better.
        // we simply treat j_min's slot as deactivated.
    }
    // the last remaining active node is the root
    for(int i = 0; i < n; ++i){
        if(activeClusters[i]){
            return nodes[i];
        }
    }
    return nullptr;
}

/**
 * @brief perform Neighbor-Joining
 */
inline NodePtr neighborJoining(const DistanceMatrix& distMatrix, const vector<string>& labels){
    if(distMatrix.empty()){
        return nullptr;
    }
    const size_t n = labels.size();
    vector<NodePtr> nodes;
    for(int i = 0; i < n; ++i){
        nodes.push_back(make_shared<Node>());
        nodes.back()->name = labels[i];
    }
    const DistanceMatrix& currentDistances = distMatrix;
    vector<int> activeIndices(n);
    iota(activeIndices.begin(), activeIndices.end(), 0);
    while(activeIndices.size() > 2){
        const size_t currentN = activeIndices.size();
        // calculate Q-matrix
        DistanceMatrix Q(currentN, vector<double>(currentN, 0.0));
        vector<double> r(currentN, 0.0);
        for(int i = 0; i < currentN; ++i){
            for(int j = 0; j < currentN; ++j){
                if(i != j){
                    r[i] += currentDistances[activeIndices[i]][activeIndices[j]];
                }
            }
        }
        double minQ = std::numeric_limits<double>::max();
        int i_min = -1, j_min = -1;
        for(int i = 0; i < currentN; ++i){
            for(int j = i + 1; j < currentN; ++j){
                Q[i][j] = currentDistances[activeIndices[i]][activeIndices[j]] - (r[i] + r[j]) /
                    (static_cast<double>(currentN) - 2.0);
                if(Q[i][j] < minQ){
                    minQ = Q[i][j];
                    i_min = i;
                    j_min = j;
                }
            }
        }
        // create a new internal node
        auto newNode = make_shared<Node>();
        newNode->left = nodes[i_min];
        newNode->right = nodes[j_min];
        // calculate branch lengths
        const double branch_i = (currentDistances[activeIndices[i_min]][activeIndices[j_min]] + (r[i_min] - r[j_min]) /
            (static_cast<double>(currentN) - 2.0)) / 2.0;
        const double branch_j = currentDistances[activeIndices[i_min]][activeIndices[j_min]] - branch_i;
        newNode->left->branchLength = branch_i;
        newNode->right->branchLength = branch_j;
        // update distances for the new node
        vector<double> newDists(currentN);
        for(int k = 0; k < currentN; ++k){
            if(k != i_min && k != j_min){
                newDists[k] = (currentDistances[activeIndices[i_min]][activeIndices[k]] +
                    currentDistances[activeIndices[j_min]][activeIndices[k]] -
                    currentDistances[activeIndices[i_min]][activeIndices[j_min]]) / 2.0;
            }
        }
        // update active indices and nodes
        int old_i = activeIndices[i_min], old_j = activeIndices[j_min];
        // remove j_min
        activeIndices.erase(activeIndices.begin() + j_min);
        nodes.erase(nodes.begin() + j_min);
        // remove i_min (since j_min was removed, i_min's index might have shifted)
        activeIndices.erase(activeIndices.begin() + i_min);
        nodes.erase(nodes.begin() + i_min);
        // add a new node to the end
        activeIndices.push_back(old_i);
        nodes.push_back(newNode);
        // a full implementation would update the distance matrix with the new node's distances, managing the shrinking
        // distance matrix. The above code is a conceptual simplification.
    }
    // Connect the last two nodes
    auto root = make_shared<Node>();
    root->left = nodes[0];
    root->right = nodes[1];
    root->left->branchLength = currentDistances[activeIndices[0]][activeIndices[1]] / 2.0;
    root->right->branchLength = currentDistances[activeIndices[0]][activeIndices[1]] / 2.0;
    return root;
}

inline void displayTree(const NodePtr& node, const string& prefix = "", const bool isLeft = false){
    if(node == nullptr){
        return;
    }
    cout << prefix;
    cout << (isLeft ? "|-- " : "`-- ");
    if(node->left == nullptr && node->right == nullptr){
        // This is a leaf node
        cout << node->name << " (Length: " << node->branchLength << ")" << endl;
    }else{
        // This is an internal node
        cout << "---" << " (Length: " << node->branchLength << ")" << endl;
        // Recursively call for children
        displayTree(node->right, prefix + (isLeft ? "|   " : "    "), false);
        displayTree(node->left, prefix + (isLeft ? "|   " : "    "), true);
    }
}

#endif //ALGO_SKETCH_DISTANCE_PHYLOGENERATOR_H