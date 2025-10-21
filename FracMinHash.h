#ifndef ALGO_SKETCH_DISTANCE_FRACMINHASH_H
#define ALGO_SKETCH_DISTANCE_FRACMINHASH_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

class FracMinHash{
public:
    // scale in (0,1]; k must be <= 31 (2 bits per base fit in 62 bits)
    FracMinHash(const std::string &filename, double scale, unsigned k, uint64_t seed = 1469598103934665603ULL);

    /**
     * @brief add one base (ACGT); other characters reset the rolling window
     */
    void add_char(char c);

    /**
     * @brief number of hashes in the sketch
     */
    size_t sketch_size() const;

    /**
     * @brief merge: sketches must use the same scale and k
     */
    // void merge(const FracMinHash &other);

    /**
     * @brief estimate Jaccard assuming same scale; returns value in [0,1]
     */
    double jaccard(const FracMinHash &other) const;

    /**
     * @brief distance = 1 - Jaccard
     */
    double distance(const FracMinHash &other) const;

    /**
     * @brief save/load binary sketch file (compact: first k, then scale as double, seed, then size and 8-byte hashes)
     */
    void save(const std::string &filename) const;
    static FracMinHash load(const std::string &filename);

private:
    std::string filename_;
    unsigned k_;
    double scale_;
    uint64_t seed_;
    uint64_t threshold_; // compare scrambled hash < threshold_

    // current rolling window
    uint64_t fw_hash_;   // forward 2-bit encoding in lower 2k bits
    uint64_t rc_hash_;   // reverse complement encoding
    unsigned filled_;    // number of consecutive valid bases in a window

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