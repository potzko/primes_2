# primes_2

`primes_2` is a Rust project that generates prime numbers with a segmented sieve of Eratosthenes over odd numbers, using bit-packed storage and a 3/5/7 presieve wheel.

## What is in this repository

- `src/main.rs` — multi-threaded prime stream (`u32` output range) plus unit tests
- `src/bin/single_u64.rs` — single-threaded `u64` variant

Both binaries currently compute and print the `100_000_001`st prime.

## Requirements

- Rust toolchain (stable)
- Cargo

## Build

```bash
cargo build
```

## Run

Run the default binary:

```bash
cargo run --release
```

Run the single-threaded `u64` variant:

```bash
cargo run --release --bin single_u64
```

## Test

```bash
cargo test
```
