#ifndef ALGO_SKETCH_DISTANCE_FRACMINHASH_H
#define ALGO_SKETCH_DISTANCE_FRACMINHASH_H

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cstring>

class FracMinHash{
public:
    // scale in (0,1]; k must be <= 31 (2 bits per base fits in 62 bits)
    FracMinHash(double scale, unsigned k, uint64_t seed = 1469598103934665603ULL);

    // streaming: feed characters from alphabet {A,C,G,T} (case-insensitive). Non-ACGT resets window.
    void add_char(char c);

    // return number of retained hashes
    size_t sketch_size() const;

    // merge: sketches must use same scale and k
    void merge(const FracMinHash &other);

    // estimate Jaccard assuming same scale; returns value in [0,1]
    double jaccard(const FracMinHash &other) const;

    // distance = 1 - Jaccard
    double distance(const FracMinHash &other) const;

    // save/load binary sketch file (compact: first k, then scale as double, seed, then size and 8-byte hashes)
    void save(const std::string &filename) const;
    static FracMinHash load(const std::string &filename);

private:
    unsigned k_;
    double scale_;
    uint64_t seed_;
    uint64_t threshold_; // compare scrambled hash < threshold_

    // current rolling window
    uint64_t fw_hash_;   // forward 2-bit encoding in lower 2k bits
    uint64_t rc_hash_;   // reverse complement encoding
    unsigned filled_;    // number of consecutive valid bases in window

    // retained scrambled hashes (the sample)
    std::unordered_set<uint64_t> sketch_;

    // helpers
    /**
     * @brief Produce a pseudo-random-looking 64-bit hash from an integer representation of the k‑mer.
     * It reduces clustering and ensures uniformity across the 0..2^64−1 space so the fractional sampling
     * threshold behaves correctly.
     */
    static inline uint64_t splitmix64(uint64_t x);
    /**
     * @brief Used to encode bases compactly into a 2‑bit-per-base integer for the rolling k‑mer window 
     * and to detect invalid characters so the rolling window can be reset.
     */
    static inline int base_to_code(char c);
    /**
     * @brief Converts the canonical k‑mer integer into the final 64‑bit value used for fractional sampling
     * and storage in the sketch. Using a seed makes sketches reproducible and allows using the same 
     * randomized mapping across different runs/sketches.
     */
    static inline uint64_t scramble(uint64_t x, uint64_t seed);

    void add_kmer_from_window();
};

#endif //ALGO_SKETCH_DISTANCE_FRACMINHASH_H