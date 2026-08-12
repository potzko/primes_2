# primes_2

[![crates.io](https://img.shields.io/crates/v/primes_2.svg)](https://crates.io/crates/primes_2)

A multi-threaded segmented sieve of Eratosthenes over odd numbers, using bit-packed
storage and a 3/5/7 presieve wheel. Yields every prime up to `u32::MAX`. No dependencies.

## Install

```bash
cargo add primes_2
```

## Usage

Two entry points.

`stream()` — a lazy `Iterator<Item = u32>` over the primes in order. Worker threads
spawn on first use past the bootstrap segment and shut down when the iterator drops,
so taking a small prefix stays cheap.

```rust
use primes_2::stream;

let first_ten: Vec<u32> = stream().take(10).collect();
assert_eq!(first_ten, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);

// The 100_000th prime.
assert_eq!(stream().nth(99_999).unwrap(), 1_299_709);
```

`PrimeSieve::up_to(limit)` — every prime `<= limit`, collected into a `Vec<u64>`.

```rust
use primes_2::PrimeSieve;

assert_eq!(PrimeSieve::up_to(30), vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
assert_eq!(PrimeSieve::up_to(1_000_000).len(), 78498);
```

### Range

The sieve's ceiling is `u32::MAX` (4_294_967_295). `stream()` ends there, and a `limit`
above it returns the same result as `u32::MAX`. For a larger range, see the
single-threaded `u64` variant below.

## Repository layout

- `src/lib.rs` — the library: `stream()`, `PrimeSieve`, and the unit tests
- `src/bin/bench.rs` — prints the `100_000_001`st prime and the time taken
- `src/bin/single_u64.rs` — standalone single-threaded `u64` variant, not part of the
  library. Its ceiling is the square of the largest bootstrap prime (~17.18B).

## Build, run, test

```bash
cargo build --release
cargo run --release --bin bench
cargo run --release --bin single_u64
cargo test
```

Requires a Rust toolchain supporting edition 2024 (1.85+).

## License

Apache-2.0. See [LICENSE](LICENSE).
