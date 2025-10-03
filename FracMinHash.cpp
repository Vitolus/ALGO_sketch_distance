#include "FracMinHash.h"

FracMinHash::FracMinHash(const double scale, const unsigned k, const uint64_t seed)
    : k_(k), scale_(scale), seed_(seed), fw_hash_(0), rc_hash_(0), filled_(0){
    if(k == 0 || k > 31) throw std::invalid_argument("k must be in 1..31");
    if(!(scale > 0.0 && scale <= 1.0)) throw std::invalid_argument("scale must be in (0,1]");
    // set threshold as floor(scale * 2^64)
    const long double s = scale_;
    constexpr uint64_t max64 = std::numeric_limits<uint64_t>::max();
    const long double thr_ld = s * static_cast<long double>(max64);
    threshold_ = static_cast<uint64_t>(thr_ld);
    if(threshold_ == 0) threshold_ = 1; // avoid zero thresholds
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

void FracMinHash::merge(const FracMinHash &other){
    if(k_ != other.k_) throw std::invalid_argument("k mismatch in merge");
    if(scale_ != other.scale_) throw std::invalid_argument("scale mismatch in merge");
    if(seed_ != other.seed_) throw std::invalid_argument("seed mismatch in merge");
    // union of sets
    for(auto h : other.sketch_) sketch_.insert(h);
}

double FracMinHash::jaccard(const FracMinHash &other) const{
    if(k_ != other.k_) throw std::invalid_argument("k mismatch in jaccard");
    if(scale_ != other.scale_) throw std::invalid_argument("scale mismatch in jaccard");
    if(seed_ != other.seed_) throw std::invalid_argument("seed mismatch in jaccard");
    if(sketch_.empty() && other.sketch_.empty()) return 1.0;
    if(sketch_.empty() || other.sketch_.empty()) return 0.0;
    // compute intersection size
    const std::unordered_set<uint64_t> *small = &sketch_;
    const std::unordered_set<uint64_t> *big = &other.sketch_;
    if(small->size() > big->size()) std::swap(small, big);
    size_t inter = 0;
    for(auto h : *small) if(big->find(h) != big->end()) ++inter; // count intersection
    const size_t uni = sketch_.size() + other.sketch_.size() - inter; // union size
    return (uni == 0) ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni); // Jaccard index
}

double FracMinHash::distance(const FracMinHash &other) const{
    double j = jaccard(other);
    if(j < 0.0) j = 0.0;
    if(j > 1.0) j = 1.0;
    return 1.0 - j;
}

void FracMinHash::save(const std::string &filename) const{
    std::ofstream out(filename, std::ios::binary);
    if(!out) throw std::runtime_error("cannot open file for writing: " + filename);
    // header: magic + version
    constexpr char magic[4] = {'F', 'M', 'H', 0};
    out.write(magic, 4);
    // k (4 bytes), scale (8 bytes double), seed (8 bytes), threshold (8 bytes)
    const uint32_t k32 = k_;
    out.write(reinterpret_cast<const char*>(&k32), sizeof(k32));
    const double scale_d = scale_;
    out.write(reinterpret_cast<const char*>(&scale_d), sizeof(scale_d));
    out.write(reinterpret_cast<const char*>(&seed_), sizeof(seed_));
    out.write(reinterpret_cast<const char*>(&threshold_), sizeof(threshold_));
    // number of hashes (8 bytes)
    const uint64_t n = sketch_.size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    // write hashes (8 bytes each) - unsorted order OK
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
    uint32_t k32;
    in.read(reinterpret_cast<char*>(&k32), sizeof(k32));
    double scale_d;
    in.read(reinterpret_cast<char*>(&scale_d), sizeof(scale_d));
    uint64_t seed;
    in.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    uint64_t threshold;
    in.read(reinterpret_cast<char*>(&threshold), sizeof(threshold));
    uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    FracMinHash fm(scale_d, k32, seed);
    fm.threshold_ = threshold; // ensure an exact threshold
    for (uint64_t i = 0; i < n; ++i) {
        uint64_t h;
        in.read(reinterpret_cast<char*>(&h), sizeof(h));
        if (!in) throw std::runtime_error("unexpected EOF while reading sketch");
        fm.sketch_.insert(h);
    }
    return fm;
}