#ifndef ALGO_SKETCH_DISTANCE_FRACMINHASH_H
#define ALGO_SKETCH_DISTANCE_FRACMINHASH_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include "BloomFilter.h"

class FracMinHash{
public:
    // scale in (0,1]; k must be <= 31
    // bloom_size_bits: size of the bloom filter in bits. 400k bits = 50KB.
    // bloom_num_hashes: number of hash functions for the bloom filter.
    FracMinHash(std::string filename, double scale, uint8_t k, uint64_t seed = 1469598103934665603ULL,
        uint64_t bloom_size_bits = 400000, uint8_t bloom_num_hashes = 10);

    /**
     * @brief add one base (ACGT); other characters reset the rolling window
     */
    void add_char(char c);

    /**
     * @brief add a sequence of bases from a buffer
     */
    void add_sequence(const char* seq, size_t len);

    /**
     * @brief number of items added to the sketch
     */
    [[nodiscard]] size_t sketch_size() const;

    /**
     * @brief merge: sketches must use the same scale and k
     */
    void merge(const FracMinHash &other);

    /**
     * @brief estimate Jaccard assuming same scale; returns value in [0,1]
     */
    [[nodiscard]] double jaccard(const FracMinHash &other) const;

    /**
     * @brief distance = 1 - Jaccard
     */
    [[nodiscard]] double distance(const FracMinHash &other) const;

    /**
     * @brief save/load binary sketch file (compact: first k, then scale as double, seed, then bloom filter data)
     */
    void save(const std::string &filename) const;
    static FracMinHash load(const std::string &filename);

private:
    std::string filename_;
    uint8_t k_;
    double scale_;
    uint64_t seed_;
    uint64_t threshold_; // compare scrambled hash < threshold_

    // current rolling window
    uint64_t fw_hash_;   // forward 2-bit encoding in lower 2k bits
    uint64_t rc_hash_;   // reverse complement encoding
    unsigned filled_;    // number of consecutive valid bases in a window

    // Bloom filter for storing the sample of scrambled hashes
    BloomFilter sketch_;
    // Keep a separate count of items added, as the bloom filter doesn't store this.
    size_t sketch_item_count_;

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

    // This function is no longer needed as its body is directly inlined into add_sequence
    // inline void add_kmer_from_window();
};

#endif //ALGO_SKETCH_DISTANCE_FRACMINHASH_H