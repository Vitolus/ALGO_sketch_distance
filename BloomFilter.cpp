#include "BloomFilter.h"
#include <atomic>

BloomFilter::BloomFilter(uint64_t size_in_bits, uint8_t num_hashes)
    : size_in_bits_(size_in_bits), num_hashes_(num_hashes) {
    // Allocate enough 64-bit integers to hold size_in_bits
    bits_.resize((size_in_bits + 63) / 64, 0);
}

uint64_t BloomFilter::hash(const uint64_t& item, uint8_t i) {
    // Use different seeds for different hash functions
    // This is a simple way to get multiple hashes from one item
    return std::hash<uint64_t>{}(item ^ (i * 0x9e3779b97f4a7c15ULL));
}

bool BloomFilter::add(const uint64_t& item) {
    bool new_bit_set = false;
    for (uint8_t i = 0; i < num_hashes_; ++i) {
        const uint64_t h = hash(item, i) % size_in_bits_;
        const uint64_t block_idx = h / 64;
        const uint64_t bit_mask = 1ULL << (h % 64);
        // Check if the bit is not already set
        if ((bits_[block_idx] & bit_mask) == 0) {
            new_bit_set = true;
            // Non-atomic OR is sufficient for single-threaded creation.
            // For multi-threading, this would require a lock or C++20's std::atomic_ref.
            bits_[block_idx] |= bit_mask;
        }
    }
    return new_bit_set;
}

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
