#pragma once
// File ncTools.h

#ifndef NC_TOOLS_H
#define NC_TOOLS_H

/*
* ncTools -- Non Cryptographic tools
* 
* The goal of ncTools is to provide high quality, efficient non-cryptographic tools. This initial
* set will be expanded as time allows. The random number generators have been tested with PractRand,
* and the hash functions have been tested with SMHasher, rurban version.
* 
* Contents
*		wyrand ... ultra-fast 64-bit state RNG
*		RNG256 ... large-state, long-period RNG
*		CompactHash ... 10 GBPS hasher
*		CompactHash_streaming ... streaming version of CompactHash
*/

#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <intrin.h> // _umul128
#include <limits>
#include <random>   // std::random_device
#include <string_view>
#include <vector>

///////////////
// utilities //
///////////////
namespace ncTools {

	// wymix mixing functions
	// ref: https://github.com/Nicoshev/rapidhash
    #if _MSC_VER
        [[nodiscard]] inline uint64_t wymix_fast(uint64_t a, uint64_t b) noexcept {
            uint64_t hi;
            uint64_t lo = _umul128(a, b, &hi);
            return (lo ^ hi);
        }
        [[nodiscard]] inline uint64_t wymix_protected(uint64_t a, uint64_t b) noexcept {
            uint64_t hi;
            uint64_t lo = _umul128(a, b, &hi);
            return (lo ^ hi) ^ (a ^ b);
        }
        /**
         * @brief Modified mixing function to preserve full state involvement.
         *
         * Borrowed from rapidhash protected mode, but with one deliberate change:
         *
         * When this mixer is called with B = A ^ SECRET (as in the first mixing step),
         * the original protected-mode final fold
         *
         *     (lo ^ hi) ^ (A ^ B)
         *
         * collapses to:
         *
         *     (lo ^ hi) ^ (A ^ (A ^ SECRET))
         *   = (lo ^ hi) ^ SECRET
         *
         * This discards any contribution from A/state in the final xor, weakening diffusion.
         *
         * To restore mixing with the full state value, we replace the xor-fold with:
         *
         *     A ^ (lo ^ hi)
         *
         * This keeps A directly involved in the output while maintaining strong avalanche.
         */
        [[nodiscard]] inline uint64_t wymix_modified(uint64_t a, uint64_t b) noexcept {
		    uint64_t hi;
		    uint64_t lo = _umul128(a, b, &hi);
		    return (lo ^ hi) ^ (a);
	    }
    #else
        [[nodiscard]] inline uint64_t wymix_fast(uint64_t a, uint64_t b) noexcept {
            unsigned __int128 product = (unsigned __int128)a * b;
            uint64_t lo = (uint64_t)product;
            uint64_t hi = (uint64_t)(product >> 64);
            return (lo ^ hi);
        }
        [[nodiscard]] inline uint64_t wymix_protected(uint64_t a, uint64_t b) noexcept {
            unsigned __int128 product = (unsigned __int128)a * b;
            uint64_t lo = (uint64_t)product;
            uint64_t hi = (uint64_t)(product >> 64);
            return (lo ^ hi) ^ (a ^ b);
        }
        [[nodiscard]] inline uint64_t wymix_modified(uint64_t a, uint64_t b) noexcept {
		    unsigned __int128 product = (unsigned __int128)a * b;
		    uint64_t lo = (uint64_t)product;
		    uint64_t hi = (uint64_t)(product >> 64);
		    return (lo ^ hi) ^ (a);
	    }
    #endif

	// John Maiga's mixing function
	// ref: https://jonkagstrom.com/mx3/mx3_rev2.html
	uint64_t mx3(uint64_t x) {
		x ^= x >> 32;
		x *= 0xbea225f9eb34556d;
		x ^= x >> 29;
		x *= 0xbea225f9eb34556d;
		x ^= x >> 32;
		x *= 0xbea225f9eb34556d;
		x ^= x >> 29;
		return x;
	}
}



////////////
// wyrand //
////////////
namespace ncTools {
    /**
     * @file wyrand.h
     * @brief A fast, high-quality 64-bit non-cryptographic pseudo-random number generator.
     *
     * This implementation is heavily inspired by Wang Yi's `wyrand` (the foundation of wyhash and
     * rapidhash), but includes a deliberate modification to the mixing step.
     *
     * **Important disclaimer**
     * Do not confuse this with the `wyrand` used in current rapidhash. The modification made here
     * (in `rapid_mix`) is the author's sole responsibility. Any defect in this variant does not
     * imply a weakness in Wang Yi's original design.
     *
     * This version has passed PractRand testing to **4 TB** with no failures reported.
     * It is suitable for many non-cryptographic applications including:
     * - Monte Carlo simulations
     * - Procedural content generation
     * - Randomized algorithms and data structures
     * - Games and simulations
     *
     * **Key characteristics**
     * - Extremely low latency (~1–2 cycles per call with good inlining on modern x86/ARM)
     * - Full 2⁶⁴ period (guaranteed by the odd Weyl-sequence increment)
     * - Minimal state: only 8 bytes
     * - Built-in support for creating large numbers of statistically independent streams
     * - Unbiased bounded random integers via rejection sampling
     * - Statically constexpr where meaningful
     *
     * **Not suitable for cryptographic or security-sensitive applications.**
     *
     * @note The core mixing function (`rapid_mix`) is a slight variant of the protected mode
     *       from an earlier rapidhash version (v1), modified to preserve better state diffusion
     *       in the specific `state` ⊕ `SECRET` feed pattern.
     *
     * ### Public Interface Summary
     *
     * **Constructors**
     * - `wyrand()`                       → non-deterministic seed (std::random_device)
     * - `explicit wyrand(uint64_t seed)` → deterministic seed
     *
     * **Random Number Generation**
     * - `result_type operator()()`      → next 64-bit value in [0, UINT64_MAX]
     * - `uniform(uint64_t limit)`       → unbiased value in [0, limit-1]
     * - `uniform(uint64_t lo, uint64_t hi)` → unbiased value in [lo, hi] (inclusive)
     *
     * **Stream Jumping / Parallelism**
     * - `discard(uint64_t n)`           → advance state by n steps
     * - `big_jump()`                    → advance by 2³² steps (~4 billion streams possible)
     * - `huge_jump()`                   → advance by 2⁴⁰ steps (~16 million streams possible)
     *
     * **Bounds (for generic code / traits)**
     * - `static constexpr min()`        → 0
     * - `static constexpr max()`        → UINT64_MAX
     *
     * **State & Seeding**
     * - `get_state()`                   → current internal state
     * - `set_state(uint64_t s)`         → set internal state directly
     * - `reseed(uint64_t seed)`         → deterministic reseed
     * - `reseed()`                      → non-deterministic reseed
     *
     * **Convenience Factory**
     * - `static std::vector<wyrand> create_multiple(size_t count, uint64_t base_seed = 0)`
     *   → creates `count` independent instances, each separated by 2⁴⁰ steps
     *     (supports up to ≈16.7 million streams)
     *
     * @see https://github.com/wangyi-fudan/wyhash (original wyhash/wyrand inspiration)
     * @see PractRand test reports for quality validation context
     *
     * Performance Test
     * ----------------
     *        Generator: wyrand
     *       cpu_id: 0
     *       rate: 12.9368GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 1
     *       rate: 12.8661GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 2
     *       rate: 13.0708GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 3
     *       rate: 13.1105GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 4
     *       rate: 12.9796GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 5
     *       rate: 12.9781GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 6
     *       rate: 13.0728GB/s
     *
     *       Generator: wyrand
     *       cpu_id: 7
     *       rate: 13.1534GB/s
     */

     /**
      * @class wyrand
      * @brief Fast 64-bit non-cryptographic PRNG with stream-splitting support.
      *
      * Usage examples:
      * @code
      *   wyrand rng;                       // non-deterministic seed
      *   uint64_t r = rng();               // get next value
      *
      *   wyrand rng2(12345);               // deterministic seed
      *
      *   auto gens = wyrand::create_multiple(5000, 98765ULL);
      *   // gens[0](), gens[1](), ... are independent streams
      * @endcode
      */
    class wyrand {
    private:
        uint64_t state;

        static constexpr uint64_t INCREMENT = 0x2d358dccaa6c78a5ull;
        static constexpr uint64_t WYRAND_SECRET = 0x8bb84b93962eacc9ull;

    public:
        using result_type = uint64_t;

        /// Default constructor: seeds non-deterministically via std::random_device.
        explicit wyrand() noexcept { reseed(); }

        /// Constructor with explicit 64-bit seed.
        explicit wyrand(uint64_t seed) noexcept { reseed(seed); }

        /// Generate the next pseudo-random 64-bit value.
        inline result_type operator()() noexcept {
            state += INCREMENT;
            return wymix_modified(state, state ^ WYRAND_SECRET);
        }

        /// Advance the state by @p n steps (equivalent to calling operator() n times, but much faster).
        inline void discard(uint64_t n) noexcept {
            state += n * INCREMENT;
        }

        /// Advance by exactly 2³² steps.
        /// Useful for creating up to ~4 billion independent streams, each safe for ~4 billion values.
        inline void big_jump() noexcept {
            discard(1ull << 32);
        }

        /// Advance by exactly 2⁴⁰ steps.
        /// Useful for creating up to ~16 million independent streams, each safe for ~1 trillion values.
        inline void huge_jump() noexcept {
            discard(1ull << 40);
        }

        static constexpr result_type min() noexcept { return 0; }
        static constexpr result_type max() noexcept { return ~uint64_t(0); }

        /// Get current internal state (for saving/checkpointing).
        inline uint64_t get_state() const noexcept { return state; }

        // Set internal state directly to an arbitrary 64-bit value.
        // 
        // @warning This bypasses all avalanche/mixing logic!
        //          Setting a simple/poorly distributed value (e.g. small integers,
        //          sequential counters, low-entropy data) can result in weak or
        //          correlated initial output sequences.
        // 
        // Use this ONLY for:
        //   - Restoring a previously saved state (via get_state())
        //   - Implementing custom stream positioning or jump tables
        //   - Advanced testing / debugging scenarios
        // 
        // For normal seeding (including from user-provided seeds), ALWAYS prefer
        // reseed(uint64_t) or the seeded constructor instead — they apply strong
        // mixing (SplitMix64) to ensure high-quality starting states.
        inline void set_state(uint64_t s) noexcept { state = s; }

        // Deterministic reseed with a 64-bit value.
        // 
        // Applies a strong SplitMix64 mixer to avalanche the input seed,
        // producing a high-quality initial state even from poor or patterned
        // seeds (small numbers, sequential values, etc.).
        // 
        // @recommended This is the preferred way to initialize or re-initialize
        //              the generator from a user-provided or computed 64-bit seed.
        inline void reseed(uint64_t seed) noexcept {
            // Use SplitMix64 mixer
            uint64_t  z = seed + 0x9e3779b97f4a7c15;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            z = z ^ (z >> 31);

            state = z;
        }

        // Non-deterministic reseed using entropy from std::random_device.
        // 
        // Internally mixes the raw entropy bits using the same strong avalanche
        // as reseed(uint64_t) to ensure good starting quality.
        inline void reseed() noexcept {
            std::random_device rd;
            uint64_t raw = ((uint64_t)rd() << 32) | rd();
            reseed(raw);  // reuse the mixer
        }

        /**
         * @brief Uniform random integer in [0, limit-1].
         *
         * Produces a perfectly unbiased value using rejection sampling.
         * Avoids the bias present in a naive `rng() % limit`.
         *
         * Special cases:
         *  - limit == 1         → always returns 0
         *  - limit == UINT64_MAX → returns full 64-bit value directly (no bias, no rejection)
         *  - limit > UINT64_MAX  → undefined behavior (do not use)
         *
         * @param limit Must be > 0 and ≤ UINT64_MAX
         * @return r such that 0 ≤ r < limit, uniformly distributed.
         */
        inline uint64_t uniform(uint64_t limit) noexcept {
            assert(limit > 0);
            if (limit <= 1) return 0;
            if (limit == UINT64_MAX) return operator()();

            uint64_t threshold = (max() / limit) * limit;
            uint64_t draw;
            do {
                draw = operator()();
            } while (draw >= threshold);

            return draw % limit;
        }

        /**
         * @brief Uniform random integer in [lo, hi], inclusive.
         *
         * Produces a perfectly unbiased value in the closed interval [lo, hi].
         * Internally delegates to the single-argument uniform() for unbiased sampling.
         *
         * @param lo Lower bound (inclusive)
         * @param hi Upper bound (inclusive)
         * @return r such that lo ≤ r ≤ hi, uniformly distributed.
         * @note If lo > hi, the arguments are automatically swapped.
         * @note The full range [0, UINT64_MAX] is supported without overflow issues.
         */
        inline uint64_t uniform(uint64_t lo, uint64_t hi) noexcept {
            if (lo > hi) std::swap(lo, hi);

            // Special case: full 64-bit range → avoid overflow in hi - lo + 1
            if (lo == 0 && hi == UINT64_MAX) {
                return operator()();
            }

            // Normal case: range size fits safely in uint64_t
            uint64_t range = hi - lo + 1;  // safe now that full-range is handled
            return lo + uniform(range);
        }

        // Return a uniformly distributed double in the range [0.0, 1.0)
        inline double uni() noexcept {
            uint64_t j = (*this)();

            // Clear 12 bits (1 sign bit, 11 exponent bits).
            // Keep the lower 52 bits.
            j &= 0x000fffffffffffffull;

            // Inject the exponent 1023, proper for numbers between 1 and 2
            // The exponent = 1023 (1.0)
            j |= 0x3ff0000000000000ull;

            // Reinterpret bits as a double, and shift 
            double d;
            memcpy(&d, &j, 8);
            return d - 1.0;// Shift from [1, 2) to [0, 1)
        }

        /**
         * @brief Factory to create multiple independent wyrand instances.
         *
         * Each instance is separated by 2⁴⁰ steps in the sequence.
         * Supports up to 2²⁴ (≈16.7 million) streams.
         * Each stream can safely generate up to roughly 2⁴⁰ values with extremely
         * low risk of detectable overlap or correlation (based on 4 TB single-stream PractRand).
         *
         * @param count     Number of generators to create (1 ≤ count ≤ 2²⁴)
         * @param base_seed Optional base seed (defaults to 0)
         * @return Vector of independent wyrand objects
         *
         * Example:
         * @code
         *   auto gens = wyrand::create_multiple(10000, 0xdeadbeefULL);
         *   // Use gens[i]() in parallel threads/tasks
         * @endcode
         */
        static std::vector<wyrand> create_multiple(size_t count, uint64_t base_seed = 0) {
            assert(count > 0 && count <= (1ull << 24)); // enforce meaningful limit

            std::vector<wyrand> gens;
            gens.reserve(count);

            wyrand current(base_seed);
            gens.push_back(current);                    // first stream: base seed

            for (size_t i = 1; i < count; ++i) {
                current.huge_jump();                    // advance 2⁴⁰ steps for next independent stream
                gens.push_back(current);
            }

            return gens;
        }
    };
}

/////////////////
// CompactHash //
/////////////////
namespace ncTools {
	// Constants derived from the first ~50 8-byte chunks of the fractional part of φ (golden ratio).
	// These are not specially tuned — almost any well-distributed 64-bit constants would likely work.
	// The real quality comes from wymix + counter-based lanes rather than constant choice.
	// 
	// Source of phi digits:  https://www.numberworld.org/constants.html
	// --- Seed Mixing ---
	constexpr static uint64_t SEED_MIX_0 = 0xddcb5d18576d77b6ull;
	constexpr static uint64_t SEED_MIX_1 = 0x00faa684ecba752dull;

	// --- Counter Initialization ---
	constexpr static uint64_t INITIAL_COUNTER_0 = 0xf39cc0605cedc834ull;
	constexpr static uint64_t INITIAL_COUNTER_1 = 0x1082276bf3a27251ull;

	// --- Per-block Counter Increments (must be odd) ---
	constexpr static uint64_t COUNTER_INCREMENT_0 = 0x9e3779b97f4a7c15ull;
	constexpr static uint64_t COUNTER_INCREMENT_1 = 0xa2e834c5893a39ebull;

	// --- Finalization Mixing ---
	constexpr static uint64_t FINAL_MIX_0 = 0x0347045b5bf1827full;
	constexpr static uint64_t FINAL_MIX_1 = 0x01886f0928403002ull;

	/**
		 * Fast, non-cryptographic 128-bit hash function with good avalanche properties.
		 *
		 * - 100% C++23 compliant, ~30 LOC (not counting constants and wymix definition)
		 * - Borrows mixing function from wyhash/rapidhash, with PHI-derived constants.
		 * - Simple 2-lane architecture
		 * - Bulk throughput: ~9.8-10.0 GB/s @ 3 GHz on x64 (unrolled/aligned cases)
		 * - Small key latency: ~29-34 cycles (1–31 bytes)
		 * - Uses Microsoft _umul128 intrinsic or GCC/Clang __int128 for 64x64->128 bit multiplication.
		 * - Passes full rurban smhasher suite with zero collisions and worst-case avalanche bias < 0.85%
		 */
	[[nodiscard]] inline std::array<uint64_t, 2> CompactHash(const void* x, const size_t size, const uint64_t seed = 0) noexcept {
		/*
		Performance:
			Bulk: 9.8-10.0 GB/s
			Small Key: 29-34 cycles/hash
		*/
		assert(x != nullptr || size == 0);

		// Initial state: asymmetric seed mixing (one XOR, one ADD) to break symmetry early
		uint64_t state[2] = { SEED_MIX_0 ^ seed, SEED_MIX_1 + seed };
		uint64_t counter[2] = { INITIAL_COUNTER_0 ^ seed, INITIAL_COUNTER_1 + seed };

		const uint8_t* data = static_cast<const uint8_t*>(x);
		size_t remaining = size;

		// Main input loop: process 16 byte chunks
		while (remaining >= 16) {
			uint64_t word[2];
			memcpy(word, data, 16);
			// Mix current block + incrementing per-lane counter (odd increments prevent short cycles)
			state[0] = wymix_protected(state[0], word[0] ^ (counter[0] += COUNTER_INCREMENT_0));
			state[1] = wymix_protected(state[1], word[1] ^ (counter[1] += COUNTER_INCREMENT_1));
			remaining -= 16;
			data += 16;
		}

		// Tail handling: pad with 0x80 byte + inject folded length into the **other** lane
		// → ensures length never overwrites the padding bit
		uint64_t tail[2] = { 0, 0 };
		if (remaining > 0)
			memcpy(tail, data, remaining);
		((uint8_t*)tail)[remaining] = 0x80;		// marker bit

		uint64_t L = size * 0xbea225f9eb34556dull;
		L ^= L >> 29;
		tail[remaining < 8] ^= L;  // folds length into the proper word

		// Mix tail into state
		state[0] = wymix_protected(state[0], tail[0] ^ (counter[0] += COUNTER_INCREMENT_0));
		state[1] = wymix_protected(state[1], tail[1] ^ (counter[1] += COUNTER_INCREMENT_1));

		// Finalization: cross-mix the two lanes (bidirectional) + extra counter steps
		// Ensures full dependency between lanes before output
		state[0] = wymix_protected(state[0], state[1] ^ (counter[0] += COUNTER_INCREMENT_0));
		state[1] = wymix_protected(state[1], state[0] ^ (counter[1] += COUNTER_INCREMENT_1));

		return { state[0], state[1] };
	}
}


////////////
// RNG256 //
////////////
namespace ncTools {
	/*
	 * RNG256 - A high-speed, 256-bit non-cryptographic PRNG.
	 *
	 * Logic: 256-bit Weyl Sequence transition + Folded-Multiplication Mixer.
	 * Period: 2^256 (~1.15 x 10^77 values).
	 *
	 * Throughput: ~4.2 GB/s per core (measured on Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz (3.00 GHz))
	 * Parallelism: Supports up to 2^128 independent streams via huge_jump().
	 *
	 * Statistical Quality: Passes PractRand, 2 TB. (One "Unusual" finding at 32 GB, not carried through
	 * into larger tests.)
	 * 		./RNG256.exe | ./RNG_test.exe stdin64 -tf 2 -te 1 -tlmax 2TB -multithreaded
	 *
	 * 		RNG_test using PractRand version 0.94
	 * 		RNG = RNG_stdin64, seed = unknown
	 * 		test set = expanded, folding = extra
	 *
	 * 		length= 256 megabytes (2^28 bytes), time= 2.9 seconds	no anomalies in 1151 test result(s)
	 * 		length= 512 megabytes (2^29 bytes), time= 7.9 seconds	no anomalies in 1220 test result(s)
	 * 		length= 1 gigabyte (2^30 bytes), time= 15.7 seconds		no anomalies in 1294 test result(s)
	 * 		length= 2 gigabytes (2^31 bytes), time= 29.3 seconds	no anomalies in 1368 test result(s)
	 * 		length= 4 gigabytes (2^32 bytes), time= 52.9 seconds	no anomalies in 1451 test result(s)
	 * 		length= 8 gigabytes (2^33 bytes), time= 106 seconds		no anomalies in 1543 test result(s)
	 * 		length= 16 gigabytes (2^34 bytes), time= 209 seconds	no anomalies in 1636 test result(s)
	 * 		length= 32 gigabytes (2^35 bytes), time= 382 seconds
	 *			Test Name                         Raw       Processed     Evaluation
	 * 			[Low8/32]BCFN(0+8,13-1,T)         R= +11.6  p =  1.1e-5   unusual
	 * 			...and 1716 test result(s) without anomalies
	 * 		length= 64 gigabytes (2^36 bytes), time= 787 seconds	no anomalies in 1805 test result(s)
	 * 		length= 128 gigabytes (2^37 bytes), time= 1596 seconds	no anomalies in 1902 test result(s)
	 * 		length= 256 gigabytes (2^38 bytes), time= 2907 seconds	no anomalies in 1983 test result(s)
	 * 		length= 512 gigabytes (2^39 bytes), time= 6031 seconds	no anomalies in 2058 test result(s)
	 * 		length= 1 terabyte (2^40 bytes), time= 12006 seconds	no anomalies in 2132 test result(s)
	 * 		length= 2 terabytes (2^41 bytes), time= 23540 seconds	no anomalies in 2194 test result(s)
	 *
	 *		PractRand command completed successfully
	 *
	 * Usage:
	 *	RNG256 gen(42);               // Seed with uint64_t
	 *	uint64_t val = gen();         // Generate
	 *
	 *	RNG256 stream2 = gen.huge_jump(); // Create independent parallel stream
	 */

	class RNG256 {
		using u32 = uint32_t;
		using u64 = uint64_t;

		// 256 bit state
		u64 state[4];

#if 1
		// 256-bit Weyl increment: four 64-bit limbs, each odd and derived from golden-ratio constants
		// Source: https://www.numberworld.org/constants.html
		// A high-quality, full-rank, well-diffused increment is essential for the lightweight
		// rapid_mix fold to produce statistically excellent output in this wide additive design.
		// See negative control below for what happens with a trivial increment.
		static constexpr u64 INCR[] = {
			0x9e3779b97f4a7c15ull,
			0xf39cc0605cedc834ull,
			0x1082276bf3a27251ull,
			0xf86c6a11d0c18e95ull
		};
#else
		// Negative control: trivial increment {1,0,0,0}
		// Reduces generator to near-trivial 64-bit counter with frozen upper 192 bits
		// PractRand (single-stream) fails dramatically at only 64 MB:
		//   Extreme BCFN failures (p ≈ 0 to 1e-9947), mod3n, BRank, TMFn, DC6,
		//   FPF, Low1/8/16/32/64, Gap-16 — massive biases and linear dependencies
		// Demonstrates that multi-limb diffusion in the increment is critical
		static constexpr u64 INCR[] = { 1, 0, 0, 0 };
#endif


#if 1
		// SECRET additives: four 64-bit constants injected before cross-XOR
		// Zeroing them out produced clean 2 GB output (single-stream),
		// suggesting they are not strictly required at modest scales.
		// However, they are extremely cheap (~4 add instructions) and may help
		// at very large volumes (64 GB+) by breaking subtle linear/low-entropy patterns.
		// Retained for robustness until proven unnecessary in larger runs.
		static constexpr u64 SECRET[] = { 0x2767f0b153d27b7full, 0x0347045b5bf1827full, 0x01886f0928403002ull, 0xc1d64ba40f335e36ull };
#else
		static constexpr u64 SECRET[] = { 0,0,0,0 };
#endif


	public:
		using result_type = u64;
		static inline u64 min() noexcept {
			return 0;
		}
		static inline u64 max() noexcept {
			return std::numeric_limits<result_type>::max();
		}

		//
		// Constructors
		// 

		// Construct with 1 seed value
		RNG256(u64 seed) {
			reseed(seed);
		}

		// Construct from a seed array
		template <size_t N>
		RNG256(std::array<uint64_t, N> seeds) {
			reseed(seeds);
		}

		// Construct from a random seed
		RNG256() {
			reseed();
		}

		//
		// Random number generation
		//

		inline u64 operator()()noexcept {
			increment();

			return wymix_fast (
				(state[0] + SECRET[0]) ^ (state[3] + SECRET[3]),
				(state[1] + SECRET[1]) ^ (state[2] + SECRET[2])
			);
		}

		// convenience RNG draws -- with explicit size
		inline u64 draw64() { return this->operator()(); }
		inline u32 draw32() { return (uint32_t)(this->operator()()); }
		/**
	 * @brief Uniform random integer in [0, limit-1].
	 *
	 * Produces a perfectly unbiased value using rejection sampling.
	 * Avoids the bias present in a naive `rng() % limit`.
	 *
	 * Special cases:
	 *  - limit == 1         → always returns 0
	 *  - limit == UINT64_MAX → returns full 64-bit value directly (no bias, no rejection)
	 *  - limit > UINT64_MAX  → undefined behavior (do not use)
	 *
	 * @param limit Must be > 0 and ≤ UINT64_MAX
	 * @return r such that 0 ≤ r < limit, uniformly distributed.
	 */
		inline uint64_t uniform(uint64_t limit) noexcept {
			assert(limit > 0);
			if (limit <= 1) return 0;
			if (limit == UINT64_MAX) return operator()();

			uint64_t threshold = (max() / limit) * limit;
			uint64_t draw;
			do {
				draw = operator()();
			} while (draw >= threshold);

			return draw % limit;
		}

		/**
		 * @brief Uniform random integer in [lo, hi], inclusive.
		 *
		 * Produces a perfectly unbiased value in the closed interval [lo, hi].
		 * Internally delegates to the single-argument uniform() for unbiased sampling.
		 *
		 * @param lo Lower bound (inclusive)
		 * @param hi Upper bound (inclusive)
		 * @return r such that lo ≤ r ≤ hi, uniformly distributed.
		 * @note If lo > hi, the arguments are automatically swapped.
		 * @note The full range [0, UINT64_MAX] is supported without overflow issues.
		 */
		inline uint64_t uniform(uint64_t lo, uint64_t hi) noexcept {
			if (lo > hi) std::swap(lo, hi);

			// Special case: full 64-bit range → avoid overflow in hi - lo + 1
			if (lo == 0 && hi == UINT64_MAX) {
				return operator()();
			}

			// Normal case: range size fits safely in uint64_t
			uint64_t range = hi - lo + 1;  // safe now that full-range is handled
			return lo + uniform(range);
		}

		// Return a uniformly distributed double in the range [0.0, 1.0)
		inline double uni() noexcept {
			uint64_t j = (*this)();

			// Clear 12 bits (1 sign bit, 11 exponent bits).
			// Keep the lower 52 bits.
			j &= 0x000fffffffffffffull;

			// Inject the exponent 1023, proper for numbers between 1 and 2
			// The exponent = 1023 (1.0)
			j |= 0x3ff0000000000000ull;

			// Reinterpret bits as a double, and shift 
			double d;
			memcpy(&d, &j, 8);
			return d - 1.0;// Shift from [1, 2) to [0, 1)
		}

		//
		// Reseed functions.
		//

		template <size_t N>
		inline void reseed(std::array<uint64_t, N> seeds) noexcept {
			auto hash = CompactHash((const void*)seeds.data(), sizeof(uint64_t) * seeds.size(), 42);
			wyrand gen(hash[0]);
			for (int i = 0; i < 4; ++i) {
				state[i] = gen();
				gen.discard(hash[1]);
			}
		}

		inline void reseed(uint64_t seed0) noexcept {
			uint64_t s = mx3(seed0 ^ 0xf06ad7ae9717877eull);  // or + some other secret
			std::array<uint64_t, 4> seeds = { s, mx3(s), mx3(s ^ 1), mx3(s ^ 2) };
			reseed(seeds);
		}

		inline void reseed() noexcept {
			std::array<uint64_t, 4> seeds;
			std::random_device rd;
			seeds[0] = ((u64)rd() << 32) | rd();
			seeds[1] = ((u64)rd() << 32) | rd();
			seeds[2] = ((u64)rd() << 32) | rd();
			seeds[3] = ((u64)rd() << 32) | rd();
			reseed(seeds);
		}

		//
		// Jump Functions
		// 

		// discard(nsteps) is equivalent to calling operator() nsteps times, but is much faster
		RNG256& discard(uint64_t nsteps) noexcept {
			if (nsteps == 0) return *this;

			for (int j = 0; j < 4; ++j) {
				u64 lo, hi;
				lo = _umul128(INCR[j], nsteps, &hi);
				addcarry(state, j, lo);
				if (j + 1 < 4)
					addcarry(state, j + 1, hi);
			}

			return *this;
		}

		// Discard an arbitrary 256 bit number of steps.
		// Note: big_jump and huge_jump are faster ways to make very large jumps, so only use this
		// when an arbitrary large step is needed (and a 2^64 or 2^128 step is insufficient).
		RNG256& discard(const std::array<u64, 4>& nsteps) noexcept {
			// product = (nsteps * INCR) mod (2^256)
			std::array<u64, 4>product = { 0 };
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < 4; ++j) {
					u64 lo, hi;
					if (i + j < 4) { // don't need higher order terms - exceeds 2^256-1
						lo = _umul128(INCR[i], nsteps[j], &hi);
						// We know the low 64 bits is needed because we are inside the if() conditional.
						addcarry(product.data(), i + j, lo);
						if (i + j + 1 < 4)
							// If the high 64 bits exceed2 2^256-1 we don't need to add it
							addcarry(product.data(), i + j + 1, hi);
					}
				}
			}

			// state += product
			unsigned char c = 0;
			c = _addcarryx_u64(c, state[0], product[0], &state[0]);
			c = _addcarryx_u64(c, state[1], product[1], &state[1]);
			c = _addcarryx_u64(c, state[2], product[2], &state[2]);
			c = _addcarryx_u64(c, state[3], product[3], &state[3]);

			return *this;
		}

		// big_jump is equivalent to calling operator() 2^64 times
		RNG256& big_jump() {
			// Multiplying INCR by 2^64 is equivalent to shifting it left by 64 bits.
			// So, we can add it to the upper 192 bits of state as follows:
			addcarry(1, INCR[0]);
			addcarry(2, INCR[1]);
			addcarry(3, INCR[2]);
			return *this;
		}

		// big_jump is equivalent to calling operator() 2^128 times
		RNG256& huge_jump() {
			// Multiplying INCR by 2^64 is equivalent to shifting it left by 128 bits.
			// So, we can add it to the upper 128 bits of state as follows:
			addcarry(2, INCR[0]);
			addcarry(3, INCR[1]);
			return *this;
		}

		/**
		 * @brief Factory to create multiple independent RNG256 instances.
		 *
		 * - Each instance is separated by 2^128 steps in the sequence.
		 * - Possible supports for up to 2^128 streams with multiple calls, but
		 *   "size_t count" is usually a 64 bit integer.
		 * - Each stream can safely generate up to roughly 2^128 values.
		 *
		 * @param count     Number of generators to create (1 ≤ count < 2^32)
		 * @param s0 Optional first base seed (defaults to 0)
		 * @param s1 Optional second base seed (defaults to 0)
		 * @return Vector of independent RNG256 objects
		 *
		 * Example:
		 * @code
		 *   auto gens = RNG256::create_multiple(10000, 1776, 1492);
		 *   // Use gens[i]() in parallel threads/tasks
		 * @endcode
		 */
		 // Create multiple RNGs, explicit seed(s)
		static std::vector<RNG256> create_multiple(
			size_t count, // how many RNG256's we will create
			u64 seed      // seed
		) {
			// For now (2026), 2^32 is a safe upper bound for vector allocation and stream count
			assert(count > 0 && count <= 0xFFFFFFFFull);

			std::vector<RNG256> gens;
			gens.reserve(count);

			RNG256 current(seed);
			gens.push_back(current);
			for (size_t i = 1; i < count; ++i) {
				current.huge_jump();// Each generator's sequence starts 2^128 steps after the previous
				gens.push_back(current);
			}
			return gens;
		}

	private:
		// Move the state one increment forward
		unsigned char increment() {
			// state += INCR
			unsigned char c = 0;
			c = _addcarryx_u64(c, state[0], INCR[0], &state[0]);
			c = _addcarryx_u64(c, state[1], INCR[1], &state[1]);
			c = _addcarryx_u64(c, state[2], INCR[2], &state[2]);
			c = _addcarryx_u64(c, state[3], INCR[3], &state[3]);
			return c;
		}


		// helper function -- assists discard()
		inline void addcarry(int branch, u64 value) noexcept {
			addcarry(state, branch, value);
		}

		/**
		 * Add a value to an arbitrary branch, and propagate carry bits upward.
		 *
		 * Performs: x += (value << (64 * branch)) with 256-bit carry propagation.
		 *
		 * Logic: Entering at 'branch' adds the initial value to that limb.
		 * By then setting value = 0, the subsequent fallthrough cases effectively
		 * become: c_out = x[n] + 0 + c_in.
		 * This ripples the carry bit through the higher-order limbs until it is
		 * absorbed or the state ends.
		 *
		 * Arguments:
		 *		x      ... pointer, valid for range x[0,..., 3]
		 *		branch ... index into x
		 *		value  ... value to be added to x[branch]
		 */
		inline void addcarry(u64* x, int branch, u64 value) noexcept {
			unsigned char c = 0;

			switch (branch) {
			case 0: c = _addcarryx_u64(c, x[0], value, &x[0]); value = 0; [[fallthrough]];
			case 1: c = _addcarryx_u64(c, x[1], value, &x[1]); value = 0; [[fallthrough]];
			case 2: c = _addcarryx_u64(c, x[2], value, &x[2]); value = 0; [[fallthrough]];
			case 3: c = _addcarryx_u64(c, x[3], value, &x[3]);
				break;
			}
		}
	};

	// TestHarness to allow PractRand evaluation of interleaved generators created with
	// the create_multiple() factory. This is particularly challenging because the factory
	// uses huge_jump, which modified only state[2] and state[3], so a much stronger dependency
	// on the quality of the mixing function is needed to pass the test than with a single
	// instance of the RNG256.
	// 
	// In other words, if your mixer is only "good enough" for single-stream usage, it will very 
	// likely fail this interleaved huge-jump test — even if it passes 5–20 TB single-stream.
	class TestRNGHuge {
		/*
		PractRand 2TB result:

			RNG_test using PractRand version 0.94
			RNG = RNG_stdin64, seed = unknown
			test set = expanded, folding = extra

			length= 256 megabytes (2^28 bytes), time= 2.5 seconds	  no anomalies in 1151 test result(s)
			length= 512 megabytes (2^29 bytes), time= 7.2 seconds	  no anomalies in 1220 test result(s)
			length= 1 gigabyte (2^30 bytes), time= 14.3 seconds		no anomalies in 1292 test result(s)
			length= 2 gigabytes (2^31 bytes), time= 26.2 seconds	  no anomalies in 1368 test result(s)
			length= 4 gigabytes (2^32 bytes), time= 46.6 seconds	  no anomalies in 1449 test result(s)
			length= 8 gigabytes (2^33 bytes), time= 91.9 seconds	  no anomalies in 1539 test result(s)
			length= 16 gigabytes (2^34 bytes), time= 177 seconds	  no anomalies in 1639 test result(s)
			length= 32 gigabytes (2^35 bytes), time= 324 seconds	  no anomalies in 1715 test result(s)
			length= 64 gigabytes (2^36 bytes), time= 684 seconds	  no anomalies in 1809 test result(s)
			length= 128 gigabytes (2^37 bytes), time= 1347 seconds	  no anomalies in 1902 test result(s)
			length= 256 gigabytes (2^38 bytes), time= 2459 seconds	  no anomalies in 1983 test result(s)
			length= 512 gigabytes (2^39 bytes), time= 5187 seconds
			  Test Name                         Raw       Processed     Evaluation
			  [Low1/32]Gap-16:B                 R=  -5.6  p =1-1.0e-4   unusual
			  ...and 2057 test result(s) without anomalies
			length= 1 terabyte (2^40 bytes), time= 10365 seconds	  no anomalies in 2132 test result(s)
			length= 2 terabytes (2^41 bytes), time= 20944 seconds	  no anomalies in 2194 test result(s)

	  */
		constexpr static unsigned int NGENS = 8; // Must be power of 2
		constexpr static uint64_t SEED = 12345;

		// Interleaved, with "huge_jump" (2^128) separation between lanes
		std::vector<RNG256> generators;

		int index = 0;

	public:
		TestRNGHuge()
			: generators(RNG256::create_multiple((size_t)NGENS, SEED))
		{
		}

		uint64_t operator()() {
			index = (index + 1) & (NGENS - 1); // Faster than % operator. This is why NGENS needs to be a power of 2
			return generators[index]();
		}
	};
}

///////////////////////////
// CompactHash_streaming //
///////////////////////////
namespace ncTools {
	class CompactHash_streaming
	{
		/*
		Class-based version of CompactHash. More versatile, and very slightly slower, than the original.
		Returns identical hash values.

		Performance:
			Bulk: 9.2-9.7 GB/s
			Small Key: 37-52 cycles/hash
		*/
		// Initial state and counter
		uint64_t state[2];
		uint64_t counter[2];
		uint64_t byte_counter;

		// input buffer
		uint8_t buffer[16];
		size_t pos;

	public:
		CompactHash_streaming(const uint64_t seed = 0) noexcept
		{
			state[0] = SEED_MIX_0 ^ seed;
			state[1] = SEED_MIX_1 + seed;
			counter[0] = INITIAL_COUNTER_0 ^ seed;
			counter[1] = INITIAL_COUNTER_1 + seed;
			byte_counter = 0;
			pos = 0;
		}

		inline CompactHash_streaming& insert(const void* v, size_t size) noexcept {
			assert(v != nullptr || size == 0);

			const uint8_t* data = reinterpret_cast<const uint8_t*>(v);
			byte_counter += size;

			// Handle case where buffer is not empty
			if (pos != 0) {
				size_t M = std::min(16 - pos, size);
				memcpy(buffer + pos, data, M);
				data += M;
				size -= M;

				if ((pos += M) == 16) {
					update_state(buffer);
					pos = 0;
				}
			}

			// If we used up all the data when writing to the buffer, we
			// can exit early.
			if (size == 0)return *this;

			// If we get here, size>0. That means the buffer was filled
			// and processed, and pos got set to 0. Since pos is 0, we can
			// bypass the buffer, inserting directly out of data.
			const uint8_t* d = data;
			const size_t NBlocks = size / 16;
			for (size_t i = 0; i < NBlocks; i++) {
				update_state(d + 16 * i);
			}
			data += NBlocks * 16;
			size -= NBlocks * 16;

			// Copy any remaining bytes to the buffer
			if (size > 0) {
				memcpy(buffer, data, size);
				pos = size;
			}

			return *this;
		}

		inline std::array<uint64_t, 2> finalize() const noexcept {
			// Working copies (const method → no mutation of object state)
			uint64_t s0 = state[0];
			uint64_t s1 = state[1];
			uint64_t c0 = counter[0];
			uint64_t c1 = counter[1];

			// Tail handling: pad with 0x80 byte + inject folded length into the **other** lane
			// → ensures length never overwrites the padding bit
			uint64_t tail[2] = { 0, 0 };
			if (pos > 0)
				memcpy(tail, buffer, pos);
			((uint8_t*)tail)[pos] = 0x80;		// marker bit

			uint64_t L = byte_counter * 0xbea225f9eb34556dull;
			L ^= L >> 29;
			tail[pos < 8] ^= L;  // folds length into the proper word

			// Mix tail into state
			s0 = wymix_protected(s0, tail[0] ^ (c0 += COUNTER_INCREMENT_0));
			s1 = wymix_protected(s1, tail[1] ^ (c1 += COUNTER_INCREMENT_1));

			// Final cross-lane mix (extra counter steps match original)
			s0 = wymix_protected(s0, s1 ^ (c0 += COUNTER_INCREMENT_0));
			s1 = wymix_protected(s1, s0 ^ (c1 += COUNTER_INCREMENT_1));

			return { s0, s1 };
		}

	private:
		inline void update_state(const uint8_t* p) noexcept {
			uint64_t x[2];
			memcpy(x, p, 16);
			state[0] = wymix_protected(state[0], x[0] ^ (counter[0] += COUNTER_INCREMENT_0));
			state[1] = wymix_protected(state[1], x[1] ^ (counter[1] += COUNTER_INCREMENT_1));
		}
	}; // class CompactHash_streaming

	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const int8_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const int16_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const int32_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const int64_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const uint8_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const uint16_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const uint32_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const uint64_t value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const float value) noexcept {
		return ch.insert(&value, sizeof(value));
	}
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const double value) noexcept {
		return ch.insert(&value, sizeof(value));
	}

	// Forward declarations of container processing

	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, std::string_view str) noexcept;
	template <class T, std::size_t N>
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const std::array<T, N>& arr) noexcept;
	template <class T, class Alloc>
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, const std::vector<T, Alloc>& vec) noexcept;

	// Definitions

	// hasher << std::string_view
	inline CompactHash_streaming& operator<<(CompactHash_streaming& ch, std::string_view str) noexcept {
		return ch.insert(str.data(), str.size());
	}

	// hasher << std::array<T,N>
	template <class T, std::size_t N>
	inline CompactHash_streaming& operator<<(
		CompactHash_streaming& ch,
		const std::array<T, N>& arr) noexcept
	{
		if constexpr (std::is_trivially_copyable_v<T>) {
			// Safe: std::array has guaranteed contiguous storage
			// and T is trivially copyable → object representation = value representation
			ch.insert(arr.data(), N * sizeof(T));
		}
		else {
			for (const T& x : arr) {
				ch << x;
			}
		}

		return ch;
	}

	// hasher << std::vector<T>
	template <class T, class Alloc> // Header is correct
	inline CompactHash_streaming& operator<<(
		CompactHash_streaming& ch,
		const std::vector<T, Alloc>& vec) noexcept
	{
		if constexpr (std::is_trivially_copyable_v<T>) {
			ch.insert(vec.data(), vec.size() * sizeof(T));
		}
		else {
			for (const T& x : vec) {
				ch << x;
			}
		}
		return ch;
	}
}

#endif // NC_TOOLS.H


// To enable this demo, define ENABLE_NCTOOLS_DEMO before including this header file.
#ifdef ENABLE_NCTOOLS_DEMO
#include <iostream>
inline void ncTools_demo() {
	// This 
	uint64_t hash_seed = 42;

	std::cout << "CompactHash of wyrand-created data:\n";
	ncTools::wyrand wy(54321);
	constexpr int N = 5;
	uint64_t data[N];
	for (int i = 0; i < N; ++i) {
		data[i] = wy();
	}
	auto h = ncTools::CompactHash(data, sizeof(data), hash_seed);
	std::cout << "CompactHash{ ";
	for (int i = 0; i < N; ++i)
		std::cout << data[i] << ((i < N - 1) ? ", " : " } ");
	std::cout << "\n\t = {";
	for (int i = 0; i < h.size(); ++i)
		std::cout << h[i] << ((i < h.size() - 1) ? ", " : " }\n");
	std::cout << "\n\n";


	std::cout << "CompactHash of RNG256-created data:\n";
	ncTools::RNG256 rng256(54321);
	for (int i = 0; i < N; ++i) {
		data[i] = rng256();
	}
	h = ncTools::CompactHash(data, sizeof(data), hash_seed);
	std::cout << "CompactHash{ ";
	for (int i = 0; i < N; ++i)
		std::cout << data[i] << ((i < N - 1) ? ", " : " }");
	std::cout << "\n\t = {";
	for (int i = 0; i < h.size(); ++i)
		std::cout << h[i] << ((i < h.size() - 1) ? ", " : " }\n");
	std::cout << "\n\n";

	std::cout << "CompactHash_streaming of RNG256-created data:\n";
	h = ncTools::CompactHash_streaming(hash_seed).insert(data, sizeof(data)).finalize();
	std::cout << "CompactHash_streaming = {";
	for (int i = 0; i < h.size(); ++i)
		std::cout << h[i] << ((i < h.size() - 1) ? ", " : " }\n");
	std::cout << "\n\n";

	std::cout << "Compare hashers for different small-seed inputs.\n";
	std::array<int, 9> sizes = { 0, 1, 5, 8, 10, 15, 16, 17, 20 };
	for (int i = 0; i < sizes.size(); i++) {
		std::cout << "size = " << sizes[i] << "\n";
		auto h1 = ncTools::CompactHash(data, sizes[i], hash_seed);
		auto h2 = ncTools::CompactHash_streaming(hash_seed).insert(data, sizes[i]).finalize();
		std::cout << "CompactHash           = {";
		for (int j = 0; j < h.size(); ++j) {
			std::cout << h1[j] << ((j < h1.size() - 1) ? ", " : " }\n");
		}
		std::cout << "CompactHash_streaming = {";
		for (int j = 0; j < h.size(); ++j) {
			std::cout << h2[j] << ((j < h2.size() - 1) ? ", " : " }\n");
		}
	}

	std::cout << "\n\nTest bulk input:\n";
	{
		std::vector<uint8_t> big(1024);
		ncTools::wyrand wy(123);
		for (auto& b : big) b = static_cast<uint8_t>(wy());
		auto h1 = ncTools::CompactHash(big.data(), big.size(), 42);
		auto h2 = ncTools::CompactHash_streaming(42).insert(big.data(), big.size()).finalize();
		std::cout << "CompactHash           = {";
		for (int i = 0; i < h1.size(); ++i) {
			std::cout << h1[i] << ((i < h1.size() - 1) ? ", " : " }\n");
		}
		std::cout << "CompactHash_streaming = {";
		for (int i = 0; i < h2.size(); ++i) {
			std::cout << h2[i] << ((i < h2.size() - 1) ? ", " : " }\n");
		}
	}
}
#endif
