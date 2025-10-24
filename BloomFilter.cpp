#include "BloomFilter.h"

BloomFilter::BloomFilter(uint64_t size_in_bits, uint8_t num_hashes)
    : bits_(size_in_bits, false), num_hashes_(num_hashes) {}

uint64_t BloomFilter::hash(const uint64_t& item, uint8_t i) {
    // Use different seeds for different hash functions
    // This is a simple way to get multiple hashes from one item
    return std::hash<uint64_t>{}(item ^ (i * 0x9e3779b97f4a7c15ULL));
}

void BloomFilter::add(const uint64_t& item) {
    for (uint8_t i = 0; i < num_hashes_; ++i) {
        bits_[hash(item, i) % bits_.size()] = true;
    }
}

bool BloomFilter::contains(const uint64_t& item) const {
    for (uint8_t i = 0; i < num_hashes_; ++i) {
        if (!bits_[hash(item, i) % bits_.size()]) {
            return false;
        }
    }
    return true;
}

void BloomFilter::merge(const BloomFilter& other) {
    if (bits_.size() != other.bits_.size() || num_hashes_ != other.num_hashes_) {
        return; // Or throw an exception
    }
    for (size_t i = 0; i < bits_.size(); ++i) {
        if (other.bits_[i]) {
            bits_[i] = true;
        }
    }
}

const std::vector<bool>& BloomFilter::get_bits() const {
    return bits_;
}

uint64_t BloomFilter::size_in_bits() const {
    return bits_.size();
}

uint8_t BloomFilter::num_hashes() const {
    return num_hashes_;
}
