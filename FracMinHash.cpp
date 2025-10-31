#include "FracMinHash.h"
#include <iostream>
#include <utility>
#include <vector>

FracMinHash::FracMinHash(std::string filename, const double scale, const uint8_t k, const uint64_t seed,
    const uint64_t bloom_size_bits, const uint8_t bloom_num_hashes)
    : filename_(std::move(filename)), k_(k), scale_(scale), seed_(seed), fw_hash_(0), rc_hash_(0), filled_(0),
    sketch_(bloom_size_bits, bloom_num_hashes), sketch_item_count_(0){
    // require 2*k <= 62 to keep mask shifts safe with 1ULL << (2*k)
    if(k == 0 || k > 31) throw std::invalid_argument("k must be in 1..31");
    if(!(scale > 0.0 && scale <= 1.0)) throw std::invalid_argument("scale must be in (0,1]");
    // compute threshold = floor(scale * 2^64).
    // use long double and handle scale==1.0 to avoid rounding/exclusion issues
    if(scale_ == 1.0){
        threshold_ = std::numeric_limits<uint64_t>::max(); // include full range
    }else{
        // scale in (0,1). ld = floor(scale * 2^64)
        const long double prod = std::floor(std::ldexp(static_cast<long double>(scale_), 64));
        if (prod < 1.0L) threshold_ = 1;
        else if (prod > static_cast<long double>(std::numeric_limits<uint64_t>::max()))
            threshold_ = std::numeric_limits<uint64_t>::max();
        else
            threshold_ = static_cast<uint64_t>(prod);
    }
}

inline uint64_t FracMinHash::splitmix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL; // add golden ratio
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; // mix bits thoroughly
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; // more mixing
    x = x ^ (x >> 31); // final avalanche
    return x;
}

inline int FracMinHash::base_to_code(const char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1; // N or other -> invalid
    }
}

inline uint64_t FracMinHash::scramble(const uint64_t x, const uint64_t seed){
    // mix kmer integer with seed using splitmix64
    return splitmix64(x ^ seed);
}

void FracMinHash::add_char(const char c){
    const int code = base_to_code(c);
    if(code < 0){
        // reset rolling window on invalid base
        fw_hash_ = rc_hash_ = 0;
        filled_ = 0;
        return;
    }
    // update forward: shift left 2 bits, add code, keep only 2k bits
    fw_hash_ = ((fw_hash_ << 2) | static_cast<uint64_t>(code)) & ((1ULL << (2*k_)) - 1ULL);
    // update reverse complement: shift right 2 bits and add complement in highest positions
    // complement: A<->T, C<->G => code_comp = 3 - code
    const unsigned comp = 3 - static_cast<unsigned>(code);
    // rc_hash_ represented in lower 2k bits as a reverse complement of the current window
    rc_hash_ = (rc_hash_ >> 2) | (static_cast<uint64_t>(comp) << (2*(k_-1)));
    if(filled_ < k_){
        ++filled_;
        if(filled_ < k_) return; // need a full window
    }
    // add canonical k-mer
    add_kmer_from_window();
}

void FracMinHash::add_kmer_from_window(){
    const uint64_t canonical = fw_hash_ < rc_hash_ ? fw_hash_ : rc_hash_;
    if(const uint64_t s = scramble(canonical, seed_); s < threshold_){
        // add() returns true if the item was novel (at least one bit was flipped).
        if (sketch_.add(s)) {
            sketch_item_count_++;
        }
    }
}

size_t FracMinHash::sketch_size() const{
    // Returns the count of items added, not the bloom filter size in bits.
    // For a more accurate cardinality estimate from a full bloom filter,
    // a formula could be used, but this direct count is better when available.
    return sketch_item_count_;
}

void FracMinHash::merge(const FracMinHash &other){
    if(k_ != other.k_) throw std::invalid_argument("k mismatch in merge");
    if(scale_ != other.scale_) throw std::invalid_argument("scale mismatch in merge");
    if(seed_ != other.seed_) throw std::invalid_argument("seed mismatch in merge");
    if(sketch_.size_in_bits() != other.sketch_.size_in_bits() || sketch_.num_hashes() != other.sketch_.num_hashes()) {
        throw std::invalid_argument("Bloom filter parameters must match for merge");
    }
    sketch_.merge(other.sketch_);
    // Note: sketch_item_count_ becomes an estimate after merging.
    // A proper calculation would require cardinality estimation.
    // For simplicity, we can sum them, but it's an overestimation.
    sketch_item_count_ += other.sketch_item_count_;
}

double FracMinHash::jaccard(const FracMinHash &other) const{
    if(k_ != other.k_) throw std::invalid_argument("k mismatch in jaccard");
    if(scale_ != other.scale_) throw std::invalid_argument("scale mismatch in jaccard");
    if(seed_ != other.seed_) throw std::invalid_argument("seed mismatch in jaccard");
    if(sketch_.size_in_bits() != other.sketch_.size_in_bits() || sketch_.num_hashes() != other.sketch_.num_hashes()) {
        throw std::invalid_argument("Bloom filter parameters must match for Jaccard estimation");
    }
    const double m = static_cast<double>(sketch_.size_in_bits());
    const double k = static_cast<double>(sketch_.num_hashes());
    if (m == 0) return 1.0; // Avoid division by zero

    // Directly use the optimized bit vectors
    const auto& bits1 = sketch_.bits_;
    const auto& bits2 = other.sketch_.bits_;
    uint64_t set_bits1 = 0;
    uint64_t set_bits2 = 0;
    uint64_t union_set_bits = 0;

    for(size_t i = 0; i < bits1.size(); ++i) {
        // Use popcount to count set bits in parallel
        set_bits1 += __builtin_popcountll(bits1[i]);
        set_bits2 += __builtin_popcountll(bits2[i]);
        union_set_bits += __builtin_popcountll(bits1[i] | bits2[i]);
    }
    if (union_set_bits == 0) return 1.0; // Both empty -> identical
    // Cardinality estimation function
    auto estimate_cardinality = [m, k](uint64_t set_bits) -> double {
        if (set_bits == 0) return 0.0;
        if (set_bits >= m) return -1.0; // Saturated, estimate is unreliable
        return -m / k * std::log(1.0 - static_cast<double>(set_bits) / m);
    };
    double card1 = estimate_cardinality(set_bits1);
    double card2 = estimate_cardinality(set_bits2);
    double card_union = estimate_cardinality(union_set_bits);
    // If any estimate is unreliable (filter is saturated), fall back to the simpler bit-based Jaccard.
    if (card1 < 0 || card2 < 0 || card_union < 0) {
        uint64_t intersection_bits = set_bits1 + set_bits2 - union_set_bits;
        return static_cast<double>(intersection_bits) / static_cast<double>(union_set_bits);
    }
    // Estimate intersection size using the inclusion-exclusion principle
    double card_intersection = card1 + card2 - card_union;
    if (card_intersection < 0) card_intersection = 0;
    if (card_union < 1.0) return 1.0; // Avoid division by zero for very small cardinalities
    double jac = card_intersection / card_union;
    std::cout << filename_<< "," << other.filename_ << "    Estimated Intersection: " << card_intersection << ", Estimated Union: " << card_union << ", Jaccard Estimate: " << jac 
    << std::endl;
    return jac;
}

double FracMinHash::distance(const FracMinHash &other) const{
    double jac = jaccard(other);
    jac = std::max(0.0, std::min(1.0, jac)); // clamp to [0,1]
    return 1.0 - jac;
}

void FracMinHash::save(const std::string &filename) const{
    std::ofstream out(filename, std::ios::binary);
    if(!out) throw std::runtime_error("cannot open file for writing: " + filename);
    // header: magic + version
    constexpr char magic[4] = {'F', 'B', 'F', 1}; // FracMinHash Bloom Filter v1
    out.write(magic, 4);
    // params: k, scale, seed
    out.write(reinterpret_cast<const char*>(&k_), sizeof(k_));
    out.write(reinterpret_cast<const char*>(&scale_), sizeof(scale_));
    out.write(reinterpret_cast<const char*>(&seed_), sizeof(seed_));
    // bloom filter params: size in bits, num hashes
    const uint64_t bf_size_bits = sketch_.size_in_bits();
    const uint8_t bf_num_hashes = sketch_.num_hashes();
    out.write(reinterpret_cast<const char*>(&bf_size_bits), sizeof(bf_size_bits));
    out.write(reinterpret_cast<const char*>(&bf_num_hashes), sizeof(bf_num_hashes));
    // Serialize the bit vector by packing bits into bytes
    const auto bits = sketch_.get_bits_for_saving();
    const size_t num_bytes = (bits.size() + 7) / 8;
    std::vector<char> packed(num_bytes, 0);
    for(size_t i = 0; i < bits.size(); ++i) {
        if (bits[i]) {
            packed[i / 8] |= (1 << (i % 8));
        }
    }
    out.write(packed.data(), packed.size());
    out.close();
}

FracMinHash FracMinHash::load(const std::string &filename){
    std::ifstream in(filename, std::ios::binary);
    if(!in) throw std::runtime_error("cannot open file for reading: " + filename);
    char magic[4];
    in.read(magic, 4);
    if(in.gcount() != 4 || magic[0] != 'F' || magic[1] != 'B' || magic[2] != 'F'){
        throw std::runtime_error("invalid sketch file (magic mismatch for bloom filter format)");
    }
    uint8_t k;
    double scale;
    uint64_t seed;
    in.read(reinterpret_cast<char*>(&k), sizeof(k));
    in.read(reinterpret_cast<char*>(&scale), sizeof(scale));
    in.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    uint64_t bf_size_bits;
    uint8_t bf_num_hashes;
    in.read(reinterpret_cast<char*>(&bf_size_bits), sizeof(bf_size_bits));
    in.read(reinterpret_cast<char*>(&bf_num_hashes), sizeof(bf_num_hashes));
    FracMinHash fm(filename, scale, k, seed, bf_size_bits, bf_num_hashes);
    // Deserialize the bit vector
    const size_t num_bytes = (bf_size_bits + 7) / 8;
    std::vector<char> packed(num_bytes);
    in.read(packed.data(), num_bytes);
    if (static_cast<size_t>(in.gcount()) != num_bytes) {
        throw std::runtime_error("unexpected EOF while reading sketch bitfield");
    }
    // Manually set the bits in the new filter's internal uint64_t vector
    auto& internal_bits = fm.sketch_.bits_;
    size_t bit_count = 0;
    for(size_t i = 0; i < packed.size() && bit_count < bf_size_bits; ++i) {
        for(int j = 0; j < 8 && bit_count < bf_size_bits; ++j) {
            if ((packed[i] >> j) & 1) {
                internal_bits[bit_count / 64] |= (1ULL << (bit_count % 64));
            }
            bit_count++;
        }
    }
    // Note: sketch_item_count_ is not stored, so it will be 0 on load.
    // This is a limitation; for accurate cardinality, it would need to be stored
    // or estimated from the loaded filter.
    return fm;
}