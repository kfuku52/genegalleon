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

// GCC vector types map to native SSE2/NEON at the container's baseline ISA.
using Scores = int __attribute__((vector_size(16)));
using Residues = unsigned char __attribute__((vector_size(4)));

template<class Vector, class Scalar>
inline Vector load(const Scalar* pointer) {
    Vector value;
    std::memcpy(&value, pointer, sizeof value);
    return value;
}

inline Scores choose(Scores mask, Scores yes, Scores no) {
    return (mask & yes) | (~mask & no);
}

struct Best {
    Scores score;
    Scores state;
};

inline Best best(Scores m, Scores x, Scores y) {
    const Scores prefer_x = x >= y;
    const Scores gap = choose(prefer_x, x, y);
    const Scores prefer_m = m >= gap;
    return {choose(prefer_m, m, gap), ~prefer_m & (2 + prefer_x)};
}

inline void boundary(int* m, int* x, int* y, int i, int up, int left) {
    m[i] = impossible;
    x[i] = up;
    y[i] = left;
}

}  // namespace

#define STRINGIFY_IMPL(x) #x
#define STRINGIFY(x) STRINGIFY_IMPL(x)

extern "C" const char* gg_rescue_protein_alignment_source() {
    return STRINGIFY(GG_RESCUE_ALIGNMENT_SOURCE_SHA);
}

extern "C" int gg_rescue_protein_alignment_v1(const unsigned char* native, int n,
                                               const unsigned char* donor, int m, int* stats) {
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
            for (; k + 3 <= last - first; k += 4) {
                const int i = first + k;
                const Best diag = best(load<Scores>(older_m + i - 1), load<Scores>(older_x + i - 1),
                                       load<Scores>(older_y + i - 1));
                const Best up = best(load<Scores>(prev_m + i - 1) - gap_open,
                                     load<Scores>(prev_x + i - 1) - gap_extend,
                                     load<Scores>(prev_y + i - 1) - gap_open);
                const Best left = best(load<Scores>(prev_m + i) - gap_open,
                                       load<Scores>(prev_x + i) - gap_open,
                                       load<Scores>(prev_y + i) - gap_extend);
                const Scores residues_a = __builtin_convertvector(load<Residues>(a + k), Scores);
                const Scores residues_b = __builtin_convertvector(load<Residues>(b + k), Scores);
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
