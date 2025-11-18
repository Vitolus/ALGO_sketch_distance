#include "BloomFilter.h"
#include <atomic>

BloomFilter::BloomFilter(uint64_t size_in_bits, uint8_t num_hashes)
    : size_in_bits_(size_in_bits), num_hashes_(num_hashes) {
    // Allocate enough 64-bit integers to hold size_in_bits
    bits_.resize((size_in_bits + 63) / 64, 0);
}

// The 'hash' and 'add' methods are defined in the header file (BloomFilter.h)
// to allow the compiler to inline them, which is a critical optimization for performance.

bool BloomFilter::contains(const uint64_t& item) const {
    for (uint8_t i = 0; i < num_hashes_; ++i) {
        const uint64_t h = hash(item, i) % size_in_bits_;
        if (!((bits_[h / 64] >> (h % 64)) & 1ULL)) {
            return false;
        }
    }
    return true;
}

void BloomFilter::merge(const BloomFilter& other) {
    if (size_in_bits_ != other.size_in_bits_ || num_hashes_ != other.num_hashes_) {
        return; // Or throw an exception
    }
    for (size_t i = 0; i < bits_.size(); ++i) {
        bits_[i] |= other.bits_[i];
    }
}

std::vector<bool> BloomFilter::get_bits_for_saving() const {
    std::vector<bool> bool_bits(size_in_bits_);
    for (uint64_t i = 0; i < size_in_bits_; ++i) {
        if ((bits_[i / 64] >> (i % 64)) & 1ULL) {
            bool_bits[i] = true;
        }
    }
    return bool_bits;
}

uint64_t BloomFilter::size_in_bits() const {
    return size_in_bits_;
}

uint8_t BloomFilter::num_hashes() const {
    return num_hashes_;
}
