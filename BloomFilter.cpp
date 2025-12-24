#include "BloomFilter.h"

BloomFilter::BloomFilter(uint64_t size_in_bits, uint8_t num_hashes)
    : size_in_bits_(size_in_bits), num_hashes_(num_hashes){
    // Allocate enough 64-bit integers to hold size_in_bits
    bits_.resize((size_in_bits + 63) / 64, 0);
}

bool BloomFilter::contains(const uint64_t& item) const{
    for(uint8_t i = 0; i < num_hashes_; ++i){
        const uint64_t h = hash(item, i) % size_in_bits_;
        if(!((bits_[h / 64] >> (h % 64)) & 1ULL)){
            return false;
        }
    }
    return true;
}

void BloomFilter::merge(const BloomFilter& other){
    if(size_in_bits_ != other.size_in_bits_ || num_hashes_ != other.num_hashes_){
        return; // Or throw an exception
    }
    for(size_t i = 0; i < bits_.size(); ++i){
        bits_[i] |= other.bits_[i];
    }
}

uint64_t BloomFilter::size_in_bits() const{
    return size_in_bits_;
}

uint8_t BloomFilter::num_hashes() const{
    return num_hashes_;
}
