fn main() {
    let start = std::time::Instant::now();
    let p = primes_2::stream()
        .nth(100_000_000)
        .expect("stream exhausted before the 100_000_001st prime");
    let end = std::time::Instant::now();
    println!("Time taken: {:.2?}", end - start);
    println!("100_000_001st prime: {}", p);
}
