#include "FracMinHash.h"
#include <cmath>
#include <iostream>

FracMinHash::FracMinHash(const double scale, const unsigned k, const uint64_t seed)
    : k_(k), scale_(scale), seed_(seed), fw_hash_(0), rc_hash_(0), filled_(0){
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
        // store the scrambled hash value (store s to simplify comparisons across sketches with the same seed / scale)
        sketch_.insert(s);
    }
}

size_t FracMinHash::sketch_size() const{
    return sketch_.size();
}

/* void FracMinHash::merge(const FracMinHash &other){
    if(k_ != other.k_) throw std::invalid_argument("k mismatch in merge");
    if(scale_ != other.scale_) throw std::invalid_argument("scale mismatch in merge");
    if(seed_ != other.seed_) throw std::invalid_argument("seed mismatch in merge");
    // union of sets
    for(auto h : other.sketch_) sketch_.insert(h);
} */

double FracMinHash::jaccard(const FracMinHash &other) const{
    if(k_ != other.k_) throw std::invalid_argument("k mismatch in jaccard");
    if(scale_ != other.scale_) throw std::invalid_argument("scale mismatch in jaccard");
    if(seed_ != other.seed_) throw std::invalid_argument("seed mismatch in jaccard");
    if(sketch_.empty() && other.sketch_.empty()) return 1.0; // both empty -> identical
    if(sketch_.empty() || other.sketch_.empty()) return 0.0; // one empty -> disjoint
    // compute intersection size
    const std::unordered_set<uint64_t> *small = &sketch_;
    const std::unordered_set<uint64_t> *big = &other.sketch_;
    if(small->size() > big->size()) std::swap(small, big);
    size_t inter = 0;
    for(auto h : *small) if(big->find(h) != big->end()) ++inter; // count intersection
    size_t uni = sketch_.size() + other.sketch_.size() - inter; // union size
    uni = static_cast<double>(uni);
    double jac = (uni == 0.0) ? 0.0 : static_cast<double>(inter) / uni; // jaccard index
    std::cout << "Intersection: " << inter << ", Union: " << uni << ", Jaccard: " << jac << std::endl;
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
    constexpr char magic[4] = {'F', 'M', 'H', 1};
    out.write(magic, 4);
    // k (1 bytes), scale (8 bytes double), seed (8 bytes)
    const uint8_t k8 = k_;
    out.write(reinterpret_cast<const char*>(&k8), sizeof(k8));
    const double scale_d = scale_;
    out.write(reinterpret_cast<const char*>(&scale_d), sizeof(scale_d));
    out.write(reinterpret_cast<const char*>(&seed_), sizeof(seed_));
    // number of hashes (8 bytes)
    const uint64_t n = sketch_.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    // write hashes (8 bytes each), unsorted order OK
    for(auto h : sketch_) out.write(reinterpret_cast<const char*>(&h), sizeof(h));
    out.close();
}

FracMinHash FracMinHash::load(const std::string &filename){
    std::ifstream in(filename, std::ios::binary);
    if(!in) throw std::runtime_error("cannot open file for reading: " + filename);
    char magic[4];
    in.read(magic, 4);
    if(in.gcount() != 4 || magic[0] != 'F' || magic[1] != 'M' || magic[2] != 'H'){
        throw std::runtime_error("invalid sketch file (magic mismatch)");
    }
    uint8_t k8;
    in.read(reinterpret_cast<char*>(&k8), sizeof(k8));
    double scale_d;
    in.read(reinterpret_cast<char*>(&scale_d), sizeof(scale_d));
    uint64_t seed;
    in.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    FracMinHash fm(scale_d, k8, seed);
    for(uint64_t i = 0; i < n; i++){
        uint64_t h;
        in.read(reinterpret_cast<char*>(&h), sizeof(h));
        if (!in) throw std::runtime_error("unexpected EOF while reading sketch");
        fm.sketch_.insert(h);
    }
    return fm;
}