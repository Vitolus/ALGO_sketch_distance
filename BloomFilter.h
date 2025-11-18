#ifndef ALGO_SKETCH_DISTANCE_BLOOMFILTER_H
#define ALGO_SKETCH_DISTANCE_BLOOMFILTER_H

#include <vector>
#include <cstdint>
#include <functional>

class BloomFilter {
public:
    BloomFilter(uint64_t size_in_bits, uint8_t num_hashes);

    // Returns true if the item resulted in setting at least one new bit.
    bool add(const uint64_t& item);

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
    static uint64_t hash(const uint64_t& item, uint8_t i);
};

#endif //ALGO_SKETCH_DISTANCE_BLOOMFILTER_H