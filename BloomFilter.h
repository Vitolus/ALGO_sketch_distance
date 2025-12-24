#ifndef ALGO_SKETCH_DISTANCE_BLOOMFILTER_H
#define ALGO_SKETCH_DISTANCE_BLOOMFILTER_H

#include <cstdint>
#include <functional>
#include <vector>

class BloomFilter {
public:
    BloomFilter(uint64_t size_in_bits, uint8_t num_hashes);

    // Inlined for performance: reduces function call overhead
    /**
     * @brief adds an item to the filter by setting the corresponding bits
     */
    inline bool add(const uint64_t& item);

    [[nodiscard]] bool contains(const uint64_t& item) const;

    /**
     * @brief merges another BloomFilter into this one by performing a bitwise OR on their bit vectors
     */
    void merge(const BloomFilter& other);

    [[nodiscard]] uint64_t size_in_bits() const;

    [[nodiscard]] uint8_t num_hashes() const;

    // Allow FracMinHash to access private members for loading
    friend class FracMinHash;

private:
    std::vector<uint64_t> bits_;
    uint64_t size_in_bits_;
    uint8_t num_hashes_;

    // Inlined for performance: reduces function call overhead.
    /**
     * @brief Generate the i-th hash for the given item using double hashing
     */
    static inline uint64_t hash(const uint64_t& item, uint8_t i);
};

inline uint64_t BloomFilter::hash(const uint64_t& item, uint8_t i){
    // Use different seeds for different hash functions
    // This is a simple way to get multiple hashes from one item
    return std::hash<uint64_t>{}(item ^ (i * 0x9e3779b97f4a7c15ULL));
}

inline bool BloomFilter::add(const uint64_t& item){
    bool new_bit_set = false;
    for(uint8_t i = 0; i < num_hashes_; ++i){
        const uint64_t h = hash(item, i) % size_in_bits_;
        const uint64_t block_idx = h / 64;
        const uint64_t bit_mask = 1ULL << (h % 64);
        // Check if the bit is not already set
        if((bits_[block_idx] & bit_mask) == 0){
            new_bit_set = true;
            bits_[block_idx] |= bit_mask;
        }
    }
    return new_bit_set;
}

#endif //ALGO_SKETCH_DISTANCE_BLOOMFILTER_H