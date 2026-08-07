// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

// Benchmark-only use of the complete implementation from the pinned VROOM
// tree. This translation unit intentionally does not include CatBLST headers.

#include "Matrixnok_fft_preview.h"

#include <benchmark/benchmark.h>

#include "bls12_381_ring_params.hpp"
#include "bench_bls12_381_common.hpp"
#include "ec.hpp"
#include "fp12.hpp"
#include "scalar_mult.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    using Ring = bls12_381_params::BlsRingMatrixNoK;
    using G1Ring = bls12_381_params::BlsRingMatrixNoKEC;
    using Fp2 = FP2Ring<Ring>;
    using Curve = G1<G1Ring>;
    using Point = Curve::ProjPoint;
    using G2Fp2 = FP2Ring<G1Ring>;
    using G2Curve = G2<G1Ring>;
    using G2Point = G2Curve::ProjPoint;
    using G2Gls = GLSMult<G1Ring>;
    using Fp12 = FP12Ring<Ring>;
    using Gt = Fp12::StandardElement;
    using Scalar = std::array<uint8_t, 32>;

    constexpr const char* field_modulus_decimal =
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";
    constexpr const char* scalar_modulus_decimal =
        "52435875175126190479447740508185965837690552500527637822603658699938581184513";
    constexpr const char* left_decimal =
        "10688470598384192112760140533683975028842723510450044479374037626765189060660250327871479868795648192029502413";
    constexpr const char* right_decimal =
        "6303071894749576408239319750618994437764284500129504072810508811036909324432781187590536039296235607926728003";
    constexpr const char* third_decimal =
        "21376941196768384225520281067367950057685447020900088958748075253530378121320500655742959737591296384059004826";
    constexpr const char* fourth_decimal =
        "18909215684248729224717959251856983313292853500388512218431526433110727973298343562771608117888706823780184009";
    constexpr Scalar dense_scalar = {
        0xa5, 0x5a, 0xc3, 0x3c, 0x96, 0x69, 0xf0, 0x0f,
        0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0,
        0x0f, 0xed, 0xcb, 0xa9, 0x87, 0x65, 0x43, 0x21,
        0x10, 0x32, 0x54, 0x76, 0x08, 0x07, 0x06, 0x3f,
    };

    void BM_VROOM_MatrixNoK_ModMul(benchmark::State& state)
    {
        const Ring ring(BigInt(field_modulus_decimal, 10));
        auto left = ring.from_bigint(BigInt(left_decimal, 10));
        auto right = ring.from_bigint(BigInt(right_decimal, 10));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.modmul(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_VROOM_MatrixNoK_FP2_Mul(benchmark::State& state)
    {
        const Ring ring(BigInt(field_modulus_decimal, 10));
        const Fp2 fp2(ring);
        auto left = Fp2::StandardElement(
            ring.from_bigint(BigInt(left_decimal, 10)),
            ring.from_bigint(BigInt(right_decimal, 10)));
        auto right = Fp2::StandardElement(
            ring.from_bigint(BigInt(third_decimal, 10)),
            ring.from_bigint(BigInt(fourth_decimal, 10)));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = fp2.mul(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_VROOM_MatrixNoK_FP2_Square(benchmark::State& state)
    {
        const Ring ring(BigInt(field_modulus_decimal, 10));
        const Fp2 fp2(ring);
        auto value = Fp2::StandardElement(
            ring.from_bigint(BigInt(left_decimal, 10)),
            ring.from_bigint(BigInt(right_decimal, 10)));

        benchmark::DoNotOptimize(value);
        for (auto _ : state) {
            auto result = fp2.square(value);
            benchmark::DoNotOptimize(result);
        }
    }

    template<std::size_t batch>
    void BM_VROOM_MatrixNoK_BatchModMul(benchmark::State& state)
    {
        const Ring ring(BigInt(field_modulus_decimal, 10));
        const auto left_value = ring.from_bigint(BigInt(left_decimal, 10));
        const auto right_value = ring.from_bigint(BigInt(right_decimal, 10));
        std::array<Ring::StandardElement, batch> left;
        std::array<Ring::StandardElement, batch> right;
        left.fill(left_value);
        right.fill(right_value);

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.template batch_modmul<batch>(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    struct RootData
    {
        std::vector<Scalar> scalars;
        std::vector<MatrixNoKFftGlv4Schedule> schedules;
        Scalar scale = {};
        MatrixNoKFftGlv4Schedule scale_schedule = {};
    };

    RootData make_roots(const std::size_t n, const bool inverse)
    {
        RootData roots;
        roots.scalars.resize(n);
        roots.schedules.resize(n);
        if (!catblst_matrixnok_fft_roots(
                reinterpret_cast<unsigned char*>(roots.scalars.data()),
                roots.schedules.data(), n, inverse, roots.scale.data(),
                &roots.scale_schedule)) {
            throw std::runtime_error("failed to prepare MatrixNoK FFT roots");
        }
        return roots;
    }

    template<class T>
    void bit_reverse_preview(std::vector<T>& values)
    {
        const int n = static_cast<int>(values.size());
        for (int i = 1, j = 0; i < n; ++i) {
            int bit = n >> 1;
            for (; (j & bit) != 0; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            if (i < j) {
                std::swap(values[i], values[j]);
            }
        }
    }

    class PublicG1Multiplier
    {
    public:
        PublicG1Multiplier(const Curve& curve, const G1Ring& ring)
            : curve_(curve), beta_(ring.from_hex(
                "0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe"))
        {}

        Point multiply(const Point& point, const Scalar& scalar, const G1Ring& ring) const
        {
            constexpr int table_log_size = 4;
            constexpr int table_size = 1 << table_log_size;
            constexpr int scalar_bits = 128;
            constexpr int window_bits = table_log_size + 1;
            constexpr int windows = (scalar_bits + table_log_size) / window_bits;

            const auto [low_bytes, high_bytes] = div_by_zz(scalar);
            const std::array<ScalarForExponent<scalar_bits>, 2> scalars = {
                ScalarForExponent<scalar_bits>(low_bytes),
                ScalarForExponent<scalar_bits>(high_bytes),
            };

            const PrecomputedTable<Curve, G1Ring, table_log_size> point_table(
                point, curve_, ring);
            std::array<Point, table_size> phi_data;
            for (int i = 0; i < table_size; ++i) {
                phi_data[i] = point_table.data[i];
                phi_data[i].x = ring.modmul(phi_data[i].x, beta_);
                phi_data[i].y = ring.standard_negate(phi_data[i].y);
            }
            const PrecomputedTable<Curve, G1Ring, table_log_size> phi_table(phi_data);
            const std::array<const PrecomputedTable<Curve, G1Ring, table_log_size>*, 2> tables = {
                &point_table, &phi_table,
            };

            std::array<std::array<std::pair<uint8_t, bool>, windows>, 2> digits = {};
            for (int table = 0; table < 2; ++table) {
                bool carry = false;
                for (int window = 0; window < windows; ++window) {
                    uint64_t value = scalars[table].template get_bits<window_bits>(
                        window * window_bits);
                    value += static_cast<uint64_t>(carry);
                    digits[table][window] = tables[table]->booth_encode(value);
                    carry = digits[table][window].second;
                }
                if (carry) {
                    throw std::logic_error("unexpected G1 GLV recoding carry");
                }
            }

            Point result = curve_.zero(ring);
            for (int window = windows - 1; window >= 0; --window) {
                for (int table = 0; table < 2; ++table) {
                    const auto [index, negative] = digits[table][window];
                    if (index != 0) {
                        Point selected = tables[table]->data[index - 1];
                        if (negative) {
                            selected.y = ring.standard_negate(selected.y);
                        }
                        result = curve_.add_point(result, selected, ring);
                    }
                }
                if (window != 0) {
                    for (int bit = 0; bit < window_bits; ++bit) {
                        result = curve_.double_point(result, ring);
                    }
                }
            }
            return result;
        }

    private:
        Curve curve_;
        G1Ring::StandardElement beta_;
    };

    Point g1_add(const Point& left, const Point& right, const Curve& curve, const G1Ring& ring)
    {
        return curve.add_point(left, right, ring);
    }

    Point g1_sub(const Point& left, Point right, const Curve& curve, const G1Ring& ring)
    {
        right.y = ring.standard_negate(right.y);
        return curve.add_point(left, right, ring);
    }

    void g1_fft(
        std::vector<Point>& values,
        const RootData& roots,
        const PublicG1Multiplier& multiplier,
        const Curve& curve,
        const G1Ring& ring,
        const bool inverse)
    {
        bit_reverse_preview(values);
        const int n = static_cast<int>(values.size());
        for (int length = 2; length <= n;) {
            const int half = length / 2;
            for (int start = 0; start < n; start += length) {
                const Point u0 = values[start];
                const Point t0 = values[start + half];
                values[start] = g1_add(u0, t0, curve, ring);
                values[start + half] = g1_sub(u0, t0, curve, ring);
                for (int j = 1; j < half; ++j) {
                    const std::size_t root_index =
                        static_cast<std::size_t>(j * (n / length));
                    const Point u = values[start + j];
                    const Point t = multiplier.multiply(
                        values[start + j + half], roots.scalars[root_index], ring);
                    values[start + j] = g1_add(u, t, curve, ring);
                    values[start + j + half] = g1_sub(u, t, curve, ring);
                }
            }
            if (length == n) {
                break;
            }
            length <<= 1;
        }
        if (inverse) {
            for (Point& value : values) {
                value = multiplier.multiply(value, roots.scale, ring);
            }
        }
    }

    Gt gt_exp(
        const Gt& base,
        const MatrixNoKFftGlv4Schedule& schedule,
        const Fp12& fp12,
        const Ring& ring)
    {
        if (schedule.sched_len == 0) {
            return fp12.one(ring);
        }
        std::array<Gt, 4> basis = {
            base,
            fp12.fp12_flat_frobenius_map1(base, ring),
            fp12.fp12_flat_frobenius_map2(base, ring),
            fp12.fp12_flat_frobenius_map3(base, ring),
        };
        for (int i = 0; i < 4; ++i) {
            if (((schedule.sign_mask >> i) & 1u) != 0) {
                basis[i] = fp12.fp12_flat_conjugate(basis[i], ring);
            }
        }

        std::array<Gt, 16> table = {};
        table[0] = fp12.one(ring);
        table[1] = basis[0];
        table[2] = basis[1];
        table[4] = basis[2];
        table[8] = basis[3];
        table[3] = fp12.fp12_flat_mul(table[1], table[2], ring);
        table[5] = fp12.fp12_flat_mul(table[1], table[4], ring);
        table[6] = fp12.fp12_flat_mul(table[2], table[4], ring);
        table[7] = fp12.fp12_flat_mul(table[3], table[4], ring);
        table[9] = fp12.fp12_flat_mul(table[1], table[8], ring);
        table[10] = fp12.fp12_flat_mul(table[2], table[8], ring);
        table[11] = fp12.fp12_flat_mul(table[3], table[8], ring);
        table[12] = fp12.fp12_flat_mul(table[4], table[8], ring);
        table[13] = fp12.fp12_flat_mul(table[5], table[8], ring);
        table[14] = fp12.fp12_flat_mul(table[6], table[8], ring);
        table[15] = fp12.fp12_flat_mul(table[7], table[8], ring);

        auto out = fp12.recast_common(table[schedule.schedule[0] & 0x0f]);
        for (std::size_t op_index = 1; op_index < schedule.sched_len; ++op_index) {
            const uint8_t op = schedule.schedule[op_index];
            for (int square = 0; square < (op >> 4); ++square) {
                out = fp12.recast_common(fp12.fp12_flat_cyclotomic_sqr(out, ring));
            }
            const unsigned table_index = op & 0x0f;
            if (table_index != 0) {
                out = fp12.recast_common(
                    fp12.fp12_flat_mul(out, table[table_index], ring));
            }
        }
        return fp12.fp12_flat_mul(fp12.one(ring), out, ring);
    }

    Gt gt_add(const Gt& left, const Gt& right, const Fp12& fp12, const Ring& ring)
    {
        return fp12.fp12_flat_mul(left, right, ring);
    }

    Gt gt_sub(const Gt& left, const Gt& right, const Fp12& fp12, const Ring& ring)
    {
        return fp12.fp12_flat_mul(
            left, fp12.fp12_flat_conjugate(right, ring), ring);
    }

    void gt_fft(
        std::vector<Gt>& values,
        const RootData& roots,
        const Fp12& fp12,
        const Ring& ring,
        const bool inverse)
    {
        bit_reverse_preview(values);
        const int n = static_cast<int>(values.size());
        for (int length = 2; length <= n;) {
            const int half = length / 2;
            for (int start = 0; start < n; start += length) {
                const Gt u0 = values[start];
                const Gt t0 = values[start + half];
                values[start] = gt_add(u0, t0, fp12, ring);
                values[start + half] = gt_sub(u0, t0, fp12, ring);
                for (int j = 1; j < half; ++j) {
                    const std::size_t root_index =
                        static_cast<std::size_t>(j * (n / length));
                    const Gt u = values[start + j];
                    const Gt t = gt_exp(
                        values[start + j + half], roots.schedules[root_index],
                        fp12, ring);
                    values[start + j] = gt_add(u, t, fp12, ring);
                    values[start + j + half] = gt_sub(u, t, fp12, ring);
                }
            }
            if (length == n) {
                break;
            }
            length <<= 1;
        }
        if (inverse) {
            for (Gt& value : values) {
                value = gt_exp(value, roots.scale_schedule, fp12, ring);
            }
        }
    }

    std::string coefficient_hex(const unsigned char* bytes)
    {
        constexpr char digits[] = "0123456789abcdef";
        std::string out;
        out.reserve(96);
        for (std::size_t i = 0; i < 48; ++i) {
            out.push_back(digits[bytes[i] >> 4]);
            out.push_back(digits[bytes[i] & 0x0f]);
        }
        return out;
    }

    Point g1_generator(const G1Ring& ring)
    {
        std::array<unsigned char, 2 * 48> coordinates = {};
        catblst_matrixnok_g1_generator_coordinates(coordinates.data());
        const std::string x = coefficient_hex(coordinates.data());
        const std::string y = coefficient_hex(coordinates.data() + 48);
        return Point(
            ring.from_bigint_exact(BigInt(x.c_str(), 16)),
            ring.from_bigint_exact(BigInt(y.c_str(), 16)),
            ring.one());
    }

    std::vector<Point> make_g1_data(
        const std::size_t n, const Curve& curve, const G1Ring& ring)
    {
        std::vector<Point> out(n);
        const Point generator = g1_generator(ring);
        out[0] = generator;
        Point current = curve.double_point(generator, ring);
        for (std::size_t i = 1; i < n; ++i) {
            out[i] = current;
            current = curve.add_point(current, generator, ring);
        }
        return out;
    }

    void BM_VROOM_MatrixNoK_G1Add(benchmark::State& state)
    {
        const G1Ring ring(BigInt(field_modulus_decimal, 10));
        const Curve curve;
        const Point generator = g1_generator(ring);
        Point current = generator;
        for (auto _ : state) {
            current = curve.add_point(current, generator, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_VROOM_MatrixNoK_G1Double(benchmark::State& state)
    {
        const G1Ring ring(BigInt(field_modulus_decimal, 10));
        const Curve curve;
        Point current = g1_generator(ring);
        for (auto _ : state) {
            current = curve.double_point(current, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    G2Point g2_generator(const G1Ring& ring, const G2Fp2& fp2ring)
    {
        AffinePoint<G1Ring::StandardElement> g1;
        AffinePoint<G2Fp2::StandardElement> g2;
        bench_bls12_381::setup_g1_g2_generators(ring, g1, g2);
        return {g2.x, g2.y, fp2ring.one()};
    }

    void BM_VROOM_MatrixNoK_G2Add(benchmark::State& state)
    {
        const G1Ring ring(BigInt(field_modulus_decimal, 10));
        const G2Fp2 fp2ring(ring);
        const G2Curve curve;
        const auto generator = g2_generator(ring, fp2ring);
        auto current = generator;
        for (auto _ : state) {
            current = curve.add_point(current, generator, fp2ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_VROOM_MatrixNoK_G2Double(benchmark::State& state)
    {
        const G1Ring ring(BigInt(field_modulus_decimal, 10));
        const G2Fp2 fp2ring(ring);
        const G2Curve curve;
        auto current = g2_generator(ring, fp2ring);
        for (auto _ : state) {
            current = curve.double_point(current, fp2ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_VROOM_MatrixNoK_G2GlsCt(benchmark::State& state)
    {
        const G1Ring ring(BigInt(field_modulus_decimal, 10));
        const G2Fp2 fp2ring(ring);
        const G2Curve curve;
        const G2Gls gls(ring);
        auto current = g2_generator(ring, fp2ring);
        for (auto _ : state) {
            current = gls.multiply(
                current, dense_scalar, fp2ring, curve, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    Gt gt_generator(const Ring& ring)
    {
        std::array<unsigned char, 12 * 48> bytes = {};
        catblst_matrixnok_gt_generator_coefficients(bytes.data());
        Gt out = {};
        for (std::size_t i = 0; i < out.size(); ++i) {
            const std::string hex = coefficient_hex(bytes.data() + i * 48);
            out[i] = ring.from_bigint_exact(BigInt(hex.c_str(), 16));
        }
        return out;
    }

    std::vector<Gt> make_gt_data(const std::size_t n, const Fp12& fp12, const Ring& ring)
    {
        std::vector<Gt> out(n);
        const Gt generator = gt_generator(ring);
        out[0] = fp12.one(ring);
        Gt current = generator;
        for (std::size_t i = 1; i < n; ++i) {
            out[i] = current;
            current = fp12.fp12_flat_mul(current, generator, ring);
        }
        return out;
    }

    bool equal_g1(const Point& left, const Point& right, const G1Ring& ring)
    {
        const auto left_x = ring.modmul(left.x, right.z);
        const auto right_x = ring.modmul(right.x, left.z);
        const auto left_y = ring.modmul(left.y, right.z);
        const auto right_y = ring.modmul(right.y, left.z);
        return ring.to_bigint_exact(left_x) == ring.to_bigint_exact(right_x) &&
            ring.to_bigint_exact(left_y) == ring.to_bigint_exact(right_y);
    }

    bool is_on_curve_g1(const Point& point, const G1Ring& ring)
    {
        const BigInt modulus(field_modulus_decimal, 10);
        const BigInt x = ring.to_bigint_exact(point.x);
        const BigInt y = ring.to_bigint_exact(point.y);
        const BigInt z = ring.to_bigint_exact(point.z);
        BigInt residual = (y * y * z - x * x * x - BigInt(4) * z * z * z) % modulus;
        if (residual < 0) {
            residual += modulus;
        }
        return residual == 0;
    }

    bool equal_gt(const Gt& left, const Gt& right, const Ring& ring)
    {
        for (std::size_t i = 0; i < left.size(); ++i) {
            if (ring.to_bigint_exact(left[i]) != ring.to_bigint_exact(right[i])) {
                return false;
            }
        }
        return true;
    }

    BigInt scalar_bigint(const Scalar& scalar)
    {
        BigInt out(0);
        for (std::size_t i = scalar.size(); i-- > 0;) {
            out = (out << 8) + BigInt(static_cast<unsigned int>(scalar[i]));
        }
        return out;
    }

    bool validate_g1_round_trip(
        const Curve& curve, const PublicG1Multiplier& multiplier, const G1Ring& ring)
    {
        constexpr std::size_t n = 256;
        const RootData forward_roots = make_roots(n, false);
        const RootData inverse_roots = make_roots(n, true);
        const auto original = make_g1_data(n, curve, ring);
        auto actual = original;
        g1_fft(actual, forward_roots, multiplier, curve, ring, false);
        g1_fft(actual, inverse_roots, multiplier, curve, ring, true);
        for (std::size_t i = 0; i < n; ++i) {
            if (!equal_g1(actual[i], original[i], ring)) {
                return false;
            }
        }
        return true;
    }

    int validate_g1_basecase(
        const Curve& curve, const PublicG1Multiplier& multiplier, const G1Ring& ring)
    {
        constexpr std::size_t n = 2;
        const RootData forward_roots = make_roots(n, false);
        const RootData inverse_roots = make_roots(n, true);
        const auto original = make_g1_data(n, curve, ring);
        const GLVMult<G1Ring> reference(curve, ring);
        if (!is_on_curve_g1(original[0], ring) || !is_on_curve_g1(original[1], ring)) {
            return 1;
        }
        if (!equal_g1(original[1],
                reference.multiply(original[0], BigInt(2), ring), ring)) {
            return 2;
        }
        std::vector<Point> actual = original;
        g1_fft(actual, forward_roots, multiplier, curve, ring, false);
        if (!equal_g1(actual[0], reference.multiply(original[0], BigInt(3), ring), ring)) {
            return 3;
        }
        if (!equal_g1(actual[1], curve.negate(original[0], ring), ring)) {
            return 4;
        }
        g1_fft(actual, inverse_roots, multiplier, curve, ring, true);
        if (!equal_g1(actual[0], original[0], ring)) {
            return 5;
        }
        if (!equal_g1(actual[1], original[1], ring)) {
            return 6;
        }
        return 0;
    }

    bool validate_g1_multiplier(
        const Curve& curve, const PublicG1Multiplier& multiplier, const G1Ring& ring)
    {
        const Point generator = g1_generator(ring);
        const GLVMult<G1Ring> reference(curve, ring);
        RootData roots = make_roots(256, false);
        roots.scalars[0][0] = 0;
        roots.scalars[1].fill(0);
        roots.scalars[1][0] = 1;
        roots.scalars[2].fill(0);
        roots.scalars[2][0] = 2;
        for (const std::size_t index : {std::size_t(0), std::size_t(1),
                 std::size_t(2), std::size_t(17), std::size_t(127)}) {
            const Point expected = reference.multiply(
                generator, scalar_bigint(roots.scalars[index]), ring);
            const Point actual = multiplier.multiply(generator, roots.scalars[index], ring);
            if (!equal_g1(actual, expected, ring)) {
                return false;
            }
        }
        const RootData inverse_two = make_roots(2, true);
        const BigInt expected_inverse_two =
            (BigInt(scalar_modulus_decimal, 10) + BigInt(1)) / BigInt(2);
        if (!equal_g1(
                multiplier.multiply(generator, inverse_two.scale, ring),
                reference.multiply(generator, expected_inverse_two, ring), ring)) {
            return false;
        }
        return true;
    }

    bool validate_gt_round_trip(const Fp12& fp12, const Ring& ring)
    {
        constexpr std::size_t n = 256;
        const RootData forward_roots = make_roots(n, false);
        const RootData inverse_roots = make_roots(n, true);
        const auto original = make_gt_data(n, fp12, ring);
        auto actual = original;
        gt_fft(actual, forward_roots, fp12, ring, false);
        gt_fft(actual, inverse_roots, fp12, ring, true);
        for (std::size_t i = 0; i < n; ++i) {
            if (!equal_gt(actual[i], original[i], ring)) {
                return false;
            }
        }
        return true;
    }

    void BM_VROOM_MatrixNoK_G1Forward(benchmark::State& state)
    {
        static const int basecase_failure = [] {
            const G1Ring validation_ring(BigInt(field_modulus_decimal, 10));
            const Curve validation_curve;
            const PublicG1Multiplier validation_multiplier(validation_curve, validation_ring);
            return validate_g1_basecase(
                validation_curve, validation_multiplier, validation_ring);
        }();
        if (basecase_failure != 0) {
            static const std::string message = "VROOM MatrixNoK G1 FFT basecase failed at stage " +
                std::to_string(basecase_failure);
            state.SkipWithError(message.c_str());
            return;
        }
        static const bool multiplier_valid = [] {
            const G1Ring validation_ring(BigInt(field_modulus_decimal, 10));
            const Curve validation_curve;
            const PublicG1Multiplier validation_multiplier(validation_curve, validation_ring);
            return validate_g1_multiplier(
                validation_curve, validation_multiplier, validation_ring);
        }();
        if (!multiplier_valid) {
            state.SkipWithError("VROOM MatrixNoK public G1 multiplier check failed");
            return;
        }
        static const bool round_trip_valid = [] {
            const G1Ring validation_ring(BigInt(field_modulus_decimal, 10));
            const Curve validation_curve;
            const PublicG1Multiplier validation_multiplier(validation_curve, validation_ring);
            return validate_g1_round_trip(
                validation_curve, validation_multiplier, validation_ring);
        }();
        if (!round_trip_valid) {
            state.SkipWithError("VROOM MatrixNoK G1 FFT round trip failed");
            return;
        }

        const G1Ring ring(BigInt(field_modulus_decimal, 10));
        const Curve curve;
        const PublicG1Multiplier multiplier(curve, ring);
        const std::size_t n = static_cast<std::size_t>(state.range(0));
        const RootData roots = make_roots(n, false);
        const auto original = make_g1_data(n, curve, ring);
        std::vector<Point> data(n);

        for (auto _ : state) {
            state.PauseTiming();
            data = original;
            state.ResumeTiming();
            const auto start = std::chrono::steady_clock::now();
            g1_fft(data, roots, multiplier, curve, ring, false);
            const auto finish = std::chrono::steady_clock::now();
            state.SetIterationTime(std::chrono::duration<double>(finish - start).count());
            benchmark::DoNotOptimize(data.data());
        }
    }

    void BM_VROOM_MatrixNoK_GTInverse(benchmark::State& state)
    {
        static const bool round_trip_valid = [] {
            const Ring validation_ring(BigInt(field_modulus_decimal, 10));
            const Fp12 validation_fp12(validation_ring);
            return validate_gt_round_trip(validation_fp12, validation_ring);
        }();
        if (!round_trip_valid) {
            state.SkipWithError("VROOM MatrixNoK GT FFT round trip failed");
            return;
        }

        const Ring ring(BigInt(field_modulus_decimal, 10));
        const Fp12 fp12(ring);
        const std::size_t n = static_cast<std::size_t>(state.range(0));
        const RootData forward_roots = make_roots(n, false);
        const RootData inverse_roots = make_roots(n, true);
        auto frequency = make_gt_data(n, fp12, ring);
        gt_fft(frequency, forward_roots, fp12, ring, false);
        std::vector<Gt> data(n);

        for (auto _ : state) {
            state.PauseTiming();
            data = frequency;
            state.ResumeTiming();
            const auto start = std::chrono::steady_clock::now();
            gt_fft(data, inverse_roots, fp12, ring, true);
            const auto finish = std::chrono::steady_clock::now();
            state.SetIterationTime(std::chrono::duration<double>(finish - start).count());
            benchmark::DoNotOptimize(data.data());
        }
    }

    void BM_VROOM_MatrixNoK_FP12Mul(benchmark::State& state)
    {
        const Ring ring(BigInt(field_modulus_decimal, 10));
        const Fp12 fp12(ring);
        const auto values = make_gt_data(3, fp12, ring);
        for (auto _ : state) {
            auto result = fp12.fp12_flat_mul(values[1], values[2], ring);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_VROOM_MatrixNoK_CyclotomicSquare(benchmark::State& state)
    {
        const Ring ring(BigInt(field_modulus_decimal, 10));
        const Fp12 fp12(ring);
        const auto value = make_gt_data(2, fp12, ring)[1];
        for (auto _ : state) {
            auto result = fp12.fp12_flat_cyclotomic_sqr(value, ring);
            benchmark::DoNotOptimize(result);
        }
    }

    BENCHMARK(BM_VROOM_MatrixNoK_G1Forward)
        ->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)
        ->UseManualTime()->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_VROOM_MatrixNoK_G1Add);
    BENCHMARK(BM_VROOM_MatrixNoK_G1Double);
    BENCHMARK(BM_VROOM_MatrixNoK_G2Add);
    BENCHMARK(BM_VROOM_MatrixNoK_G2Double);
    BENCHMARK(BM_VROOM_MatrixNoK_G2GlsCt);
    BENCHMARK(BM_VROOM_MatrixNoK_FP12Mul);
    BENCHMARK(BM_VROOM_MatrixNoK_CyclotomicSquare);
    BENCHMARK(BM_VROOM_MatrixNoK_GTInverse)
        ->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)
        ->UseManualTime()->Unit(benchmark::kMillisecond);
    BENCHMARK(BM_VROOM_MatrixNoK_ModMul);
    BENCHMARK(BM_VROOM_MatrixNoK_FP2_Mul);
    BENCHMARK(BM_VROOM_MatrixNoK_FP2_Square);
    BENCHMARK_TEMPLATE(BM_VROOM_MatrixNoK_BatchModMul, 1);
    BENCHMARK_TEMPLATE(BM_VROOM_MatrixNoK_BatchModMul, 2);
    BENCHMARK_TEMPLATE(BM_VROOM_MatrixNoK_BatchModMul, 4);
    BENCHMARK_TEMPLATE(BM_VROOM_MatrixNoK_BatchModMul, 8);
    BENCHMARK_TEMPLATE(BM_VROOM_MatrixNoK_BatchModMul, 12);
}
