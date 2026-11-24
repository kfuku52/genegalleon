// Exact global protein-only affine alignment for OrionGeno rescue.
// The three states are diagonal (0), gap in the donor (1) and gap in the
// native peptide (2). A tie always chooses the lower state, including at
// the end, as in Biopython's first Gotoh global traceback.
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

namespace {

constexpr int gap_open = 10;
constexpr int gap_extend = 1;
constexpr int impossible = -100'000'000;

struct Workspace {
    std::vector<int> waves;
    std::vector<unsigned char> reversed_donor;
    std::vector<std::size_t> offsets;
    std::unique_ptr<unsigned char[]> trace;
    std::size_t capacity = 0;
};

inline int predecessor(const int m, const int x, const int y) {
    return m >= x && m >= y ? 0 : x >= y ? 1 : 2;
}

// ABI-neutral inline helpers let the same recurrence compile inside both
// baseline and ISA-specific entry points.
template<class Vector, class Scalar>
__attribute__((always_inline)) inline void load(Vector& value, const Scalar* pointer) {
    std::memcpy(&value, pointer, sizeof value);
}

template<class Scores>
struct Best {
    Scores score;
    Scores state;
};

template<class Scores>
__attribute__((always_inline)) inline void best(
    const Scores& m, const Scores& x, const Scores& y, Best<Scores>& result) {
    const Scores prefer_x = x >= y;
    const Scores gap = (prefer_x & x) | (~prefer_x & y);
    const Scores prefer_m = m >= gap;
    result.score = (prefer_m & m) | (~prefer_m & gap);
    result.state = ~prefer_m & (2 + prefer_x);
}

inline void boundary(int* m, int* x, int* y, int i, int up, int left) {
    m[i] = impossible;
    x[i] = up;
    y[i] = left;
}

template<class Scores, class Residues>
__attribute__((always_inline)) inline int align(
    const unsigned char* native, int n, const unsigned char* donor, int m, int* stats) {
    constexpr int lanes = sizeof(Scores) / sizeof(int);
    if (n <= 0 || m <= 0 || static_cast<std::int64_t>(n) * m > 25'000'000) {
        return 1;
    }
    try {
        thread_local Workspace workspace;
        const std::size_t stride = static_cast<std::size_t>(n) + 1;
        workspace.waves.resize(9 * stride);
        workspace.reversed_donor.resize(m);
        std::reverse_copy(donor, donor + m, workspace.reversed_donor.begin());
        workspace.offsets.resize(static_cast<std::size_t>(n) + m + 1);
        const auto cells = static_cast<std::size_t>(n) * m;
        if (workspace.capacity < cells) {
            workspace.trace.reset(new unsigned char[cells]);
            workspace.capacity = cells;
        }

        // DP on antidiagonals. Adjacent cells in this loop have no dependency
        // on one another and the compiler can issue SIMD loads/stores.
        int* older_m = workspace.waves.data();
        int* older_x = older_m + stride;
        int* older_y = older_x + stride;
        int* prev_m = older_y + stride;
        int* prev_x = prev_m + stride;
        int* prev_y = prev_x + stride;
        int* next_m = prev_y + stride;
        int* next_x = next_m + stride;
        int* next_y = next_x + stride;
        older_m[0] = 0;
        older_x[0] = older_y[0] = impossible;
        boundary(prev_m, prev_x, prev_y, 0, impossible, -gap_open);
        boundary(prev_m, prev_x, prev_y, 1, -gap_open, impossible);
        std::size_t offset = 0;
        for (int d = 2; d <= n + m; ++d) {
            const int first = std::max(1, d - m);
            const int last = std::min(n, d - 1);
            if (d <= m) {
                boundary(next_m, next_x, next_y, 0, impossible, -gap_open - (d - 1) * gap_extend);
            }
            if (d <= n) {
                boundary(next_m, next_x, next_y, d, -gap_open - (d - 1) * gap_extend, impossible);
            }
            workspace.offsets[d] = offset;
            auto* trace = workspace.trace.get() + offset;
            const auto* a = native + first - 1;
            const auto* b = workspace.reversed_donor.data() + (m - d + first);
            // Score independent cells in this antidiagonal together.
            int k = 0;
            for (; k + lanes - 1 <= last - first; k += lanes) {
                const int i = first + k;
                Scores from_m, from_x, from_y;
                Best<Scores> diag, up, left;
                load(from_m, older_m + i - 1);
                load(from_x, older_x + i - 1);
                load(from_y, older_y + i - 1);
                best(from_m, from_x, from_y, diag);
                load(from_m, prev_m + i - 1);
                load(from_x, prev_x + i - 1);
                load(from_y, prev_y + i - 1);
                best(from_m - gap_open, from_x - gap_extend, from_y - gap_open, up);
                load(from_m, prev_m + i);
                load(from_x, prev_x + i);
                load(from_y, prev_y + i);
                best(from_m - gap_open, from_x - gap_open, from_y - gap_extend, left);
                Residues raw_a, raw_b;
                load(raw_a, a + k);
                load(raw_b, b + k);
                const Scores residues_a = __builtin_convertvector(raw_a, Scores);
                const Scores residues_b = __builtin_convertvector(raw_b, Scores);
                const Scores score = diag.score - 2 + ((residues_a == residues_b) & 6);
                const Scores states = diag.state | (up.state << 2) | (left.state << 4);
                const Residues packed = __builtin_convertvector(states, Residues);
                std::memcpy(next_m + i, &score, sizeof score);
                std::memcpy(next_x + i, &up.score, sizeof up.score);
                std::memcpy(next_y + i, &left.score, sizeof left.score);
                std::memcpy(trace + k, &packed, sizeof packed);
            }
            for (; k <= last - first; ++k) {
                const int i = first + k;
                const int diagonal_m = older_m[i - 1];
                const int diagonal_x = older_x[i - 1];
                const int diagonal_y = older_y[i - 1];
                const int up_m = prev_m[i - 1] - gap_open;
                const int up_x = prev_x[i - 1] - gap_extend;
                const int up_y = prev_y[i - 1] - gap_open;
                const int left_m = prev_m[i] - gap_open;
                const int left_x = prev_x[i] - gap_open;
                const int left_y = prev_y[i] - gap_extend;
                next_m[i] = std::max({diagonal_m, diagonal_x, diagonal_y}) + (a[k] == b[k] ? 4 : -2);
                next_x[i] = std::max({up_m, up_x, up_y});
                next_y[i] = std::max({left_m, left_x, left_y});
                const int states = predecessor(diagonal_m, diagonal_x, diagonal_y)
                                 | (predecessor(up_m, up_x, up_y) << 2)
                                 | (predecessor(left_m, left_x, left_y) << 4);
                trace[k] = static_cast<unsigned char>(states);
            }
            offset += last - first + 1;
            std::swap(older_m, prev_m);
            std::swap(older_x, prev_x);
            std::swap(older_y, prev_y);
            std::swap(prev_m, next_m);
            std::swap(prev_x, next_x);
            std::swap(prev_y, next_y);
        }

        int i = n;
        int j = m;
        int state = predecessor(prev_m[n], prev_x[n], prev_y[n]);
        int matches = 0;
        int aligned = 0;
        int loss_a = 0;
        int loss_b = 0;
        int last_state = -1;
        int run_length = 0;
        while (i > 0 && j > 0) {
            const int d = i + j;
            const int index = i - std::max(1, d - m);
            const int trace = workspace.trace[workspace.offsets[d] + index];
            run_length = state == last_state ? run_length + 1 : 1;
            last_state = state;
            if (state == 0) {
                const auto residue = native[--i];
                const auto other = donor[--j];
                ++aligned;
                matches += residue == other && residue != 'X' && residue != 'B'
                                             && residue != 'Z' && residue != 'J';
            } else if (state == 1) {
                --i;
                loss_a |= run_length >= 9 && j < m;
            } else {
                --j;
                loss_b |= run_length >= 9 && i < n;
            }
            state = (trace >> (2 * state)) & 3;
        }
        stats[0] = matches;
        stats[1] = aligned;
        stats[2] = loss_a;
        stats[3] = loss_b;
        return 0;
    } catch (const std::bad_alloc&) {
        return 2;
    }
}


using Scorer = int (*)(const unsigned char*, int, const unsigned char*, int, int*);

// GCC vector types map to SSE2 and NEON at the supported baseline ISAs.
int baseline(const unsigned char* native, int n, const unsigned char* donor, int m, int* stats) {
    using Scores = int __attribute__((vector_size(16)));
    using Residues = unsigned char __attribute__((vector_size(4)));
    return align<Scores, Residues>(native, n, donor, m, stats);
}

#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("avx2")))
int avx2(const unsigned char* native, int n, const unsigned char* donor, int m, int* stats) {
    using Scores = int __attribute__((vector_size(32)));
    using Residues = unsigned char __attribute__((vector_size(8)));
    return align<Scores, Residues>(native, n, donor, m, stats);
}
#endif

struct Backend {
    Scorer scorer;
    int lanes;
};

const Backend& selected_backend() {
    static const Backend backend = [] {
#if defined(__x86_64__) || defined(__i386__)
        if (__builtin_cpu_supports("avx2")) {
            return Backend{&avx2, 8};
        }
#endif
        return Backend{&baseline, 4};
    }();
    return backend;
}

}  // namespace

#define STRINGIFY_IMPL(x) #x
#define STRINGIFY(x) STRINGIFY_IMPL(x)

extern "C" const char* gg_rescue_protein_alignment_source() {
    return STRINGIFY(GG_RESCUE_ALIGNMENT_SOURCE_SHA);
}

extern "C" int gg_rescue_protein_alignment_lanes() {
    return selected_backend().lanes;
}

extern "C" int gg_rescue_protein_alignment_baseline(
    const unsigned char* native, int n, const unsigned char* donor, int m, int* stats) {
    return baseline(native, n, donor, m, stats);
}

extern "C" int gg_rescue_protein_alignment_v1(
    const unsigned char* native, int n, const unsigned char* donor, int m, int* stats) {
    return selected_backend().scorer(native, n, donor, m, stats);
}
