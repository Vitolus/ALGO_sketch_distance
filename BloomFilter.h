#ifndef ALGO_SKETCH_DISTANCE_BLOOMFILTER_H
#define ALGO_SKETCH_DISTANCE_BLOOMFILTER_H

#include <vector>
#include <cstdint>
#include <functional>

class BloomFilter {
public:
    BloomFilter(uint64_t size_in_bits, uint8_t num_hashes);

    // Adds an item to the filter by setting the corresponding bits.
    // This function is branchless and unconditionally sets bits for maximum performance.
    // Inlined for performance: reduces function call overhead when called in a tight loop.
    inline void add(const uint64_t& item);

    bool contains(const uint64_t& item) const;

    // Bitwise OR another filter into this one
    void merge(const BloomFilter& other);

    std::vector<bool> get_bits_for_saving() const;

    uint64_t size_in_bits() const;

    uint8_t num_hashes() const;

    // Allow FracMinHash to access private members for loading
    friend class FracMinHash;

private:
    std::vector<uint64_t> bits_;
    uint64_t size_in_bits_;
    uint8_t num_hashes_;

    // Generate the i-th hash for the item.
    // Here we use a simple double hashing approach.
    // Inlined for performance: reduces function call overhead.
    static inline uint64_t hash(const uint64_t& item, uint8_t i);
};

// Inlined hash function definition. Must be in the header for inlining to work across files.
inline uint64_t BloomFilter::hash(const uint64_t& item, uint8_t i) {
    // Use different seeds for different hash functions by XORing with a product of a large prime.
    return std::hash<uint64_t>{}(item ^ (i * 0x9e3779b97f4a7c15ULL));
}

// Inlined add function definition. Must be in the header.
inline void BloomFilter::add(const uint64_t& item) {
    // This branchless implementation unconditionally sets bits, which is faster than checking
    // before setting due to avoiding branch misprediction penalties, especially on dense filters.
    #pragma GCC unroll 8
    for (uint8_t i = 0; i < num_hashes_; ++i) {
        const uint64_t h = hash(item, i) % size_in_bits_;
        const uint64_t block_idx = h / 64;
        const uint64_t bit_mask = 1ULL << (h % 64);
        // Unconditionally set the bit.
        bits_[block_idx] |= bit_mask;
    }
    // This function is void because cheaply and accurately determining if a new bit was set
    // adds significant overhead. The calling code is designed to handle this by treating
    // the sketch item count as an over-estimation, which is an acceptable trade-off for the performance gain.
}

#endif //ALGO_SKETCH_DISTANCE_BLOOMFILTER_H