// Protein bottom-k sketches and same-marker distances for the rescue guide tree.
// Standalone C++17; no alignment, N*M subprocesses or Python pairwise loops.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Sketch = std::vector<uint64_t>;
using Species = std::vector<Sketch>;

static uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static void put(std::ostream& out, uint64_t x) {
    for (int i = 0; i < 8; ++i) out.put(static_cast<char>((x >> (8*i)) & 255));
}
static uint64_t get(std::istream& in) {
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) {
        int c = in.get();
        if (c == EOF) throw std::runtime_error("Truncated sketch");
        x |= static_cast<uint64_t>(c) << (8*i);
    }
    return x;
}
static Sketch sketch(const std::string& sequence, size_t k, size_t size) {
    const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
    uint64_t modulus = 1, code = 0;
    for (size_t i = 1; i < k; ++i) modulus *= 20;
    size_t valid = 0;
    Sketch hashes;
    hashes.reserve(sequence.size());
    for (char c : sequence) {
        size_t value = alphabet.find(c);
        if (value == std::string::npos) { valid = 0; code = 0; continue; }
        code = (code % modulus) * 20 + value;
        if (++valid >= k) hashes.push_back(mix(code));
    }
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    if (hashes.size() > size) hashes.resize(size);
    return hashes;
}
// Sample the bottom-k of the UNION, rather than intersecting truncated sets
// with a biased denominator. Short proteins use their entire k-mer set.
static double distance(const Sketch& a, const Sketch& b, size_t k, size_t size) {
    size_t i = 0, j = 0, total = 0, shared = 0;
    while (total < size && (i < a.size() || j < b.size())) {
        if (i < a.size() && j < b.size() && a[i] == b[j]) { ++shared; ++i; ++j; }
        else if (j == b.size() || (i < a.size() && a[i] < b[j])) ++i;
        else ++j;
        ++total;
    }
    // Zero shared hashes is saturation, not a missing marker. Preserve it.
    if (!shared) return 1.0;
    double similarity = static_cast<double>(shared) / total;
    return std::min(1.0, -std::log(2*similarity/(1+similarity))/k);
}
static void make_sketch(const std::string& input, const std::string& output, size_t k, size_t size) {
    std::ifstream in(input);
    std::ofstream out(output, std::ios::binary);
    if (!in || !out) throw std::runtime_error("Cannot open sketch input/output");
    size_t markers;
    if (!(in >> markers) || markers < 1 || markers > 100000) throw std::runtime_error("Invalid marker count");
    out.write("GGKMER01", 8); put(out, k); put(out, size); put(out, markers);
    std::string line;
    std::getline(in, line);
    for (size_t g = 0; g < markers; ++g) {
        if (!std::getline(in, line)) throw std::runtime_error("Missing marker input");
        Sketch values = sketch(line, k, size);
        put(out, values.size());
        for (uint64_t x : values) put(out, x);
    }
    if (std::getline(in, line)) throw std::runtime_error("Extra marker input");
    if (!out) throw std::runtime_error("Sketch write failed");
}
static void compare(const std::string& list, const std::string& output, size_t cpus) {
    std::ifstream paths(list);
    std::string path;
    std::vector<Species> species;
    size_t k = 0, size = 0, markers = 0;
    while (std::getline(paths, path)) {
        std::ifstream in(path, std::ios::binary);
        char magic[8]; in.read(magic, 8);
        if (!in || std::string(magic, 8) != "GGKMER01") throw std::runtime_error("Invalid sketch: " + path);
        size_t sk = get(in), ss = get(in), sm = get(in);
        if (sk < 3 || sk > 10 || ss < 1 || ss > 65536 || sm < 1 || sm > 100000)
            throw std::runtime_error("Invalid sketch dimensions");
        if (!species.empty() && (k != sk || size != ss || markers != sm)) throw std::runtime_error("Incompatible sketches");
        k = sk; size = ss; markers = sm;
        Species values(markers);
        for (Sketch& hashes : values) {
            size_t count = get(in);
            if (count > size) throw std::runtime_error("Invalid sketch size");
            hashes.resize(count);
            for (uint64_t& x : hashes) x = get(in);
            if (!std::is_sorted(hashes.begin(), hashes.end()) ||
                std::adjacent_find(hashes.begin(), hashes.end()) != hashes.end()) throw std::runtime_error("Unsorted sketch");
        }
        if (in.peek() != EOF) throw std::runtime_error("Extra sketch bytes");
        species.push_back(std::move(values));
    }
    size_t n = species.size();
    if (n < 3 || n > 10000) throw std::runtime_error("Need 3..10000 species");
    struct Result { double sum[2] = {0,0}; size_t count[2] = {0,0}, saturated = 0; };
    std::vector<Result> results(n*n);
    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    for (size_t t = 0; t < std::min(cpus, n); ++t) workers.emplace_back([&]() {
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= n) break;
            for (size_t j = i+1; j < n; ++j) {
                Result& r = results[i*n+j];
                for (size_t g = 0; g < markers; ++g) {
                    const Sketch& a = species[i][g]; const Sketch& b = species[j][g];
                    if (a.empty() || b.empty()) continue;
                    double d = distance(a,b,k,size);
                    r.sum[g%2] += d; ++r.count[g%2];
                    if (d == 1.0) ++r.saturated;
                }
            }
        }
    });
    for (auto& worker : workers) worker.join();
    std::ofstream out(output);
    if (!out) throw std::runtime_error("Cannot open distances output");
    out << "i\tj\tdistance\tshared\tsaturated\tpanel_0\tcount_0\tpanel_1\tcount_1\n" << std::setprecision(17);
    for (size_t i = 0; i < n; ++i) for (size_t j = i+1; j < n; ++j) {
        const Result& r = results[i*n+j]; size_t count = r.count[0] + r.count[1];
        out << i << '\t' << j << '\t' << (count ? (r.sum[0]+r.sum[1])/count : -1) << '\t' << count << '\t' << r.saturated;
        for (int panel = 0; panel < 2; ++panel) out << '\t' << (r.count[panel] ? r.sum[panel]/r.count[panel] : -1) << '\t' << r.count[panel];
        out << '\n';
    }
    if (!out) throw std::runtime_error("Distance write failed");
}
int main(int argc, char** argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "sketch") {
            size_t k = std::stoul(argv[4]), size = std::stoul(argv[5]);
            if (k < 3 || k > 10 || size < 1 || size > 65536) throw std::runtime_error("Invalid k/sketch size");
            make_sketch(argv[2],argv[3],k,size);
        } else if (argc == 5 && std::string(argv[1]) == "compare") {
            size_t cpus = std::stoul(argv[4]);
            if (!cpus || cpus > 1024) throw std::runtime_error("Invalid CPU count");
            compare(argv[2],argv[3],cpus);
        } else throw std::runtime_error("Usage: gg-kmer-distance sketch INPUT OUTPUT K SIZE | compare LIST OUTPUT CPUS");
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
