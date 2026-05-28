// Segmented sieve of Eratosthenes over odd numbers, packed into bits.
// Each bit b in a segment with first odd number `seg_lo` represents `seg_lo + 2*b`.
const SEGSIZE: usize = 1 << 16;
const NUMBER_WORDS: usize = SEGSIZE / 64;

// The 3-5-7 presieve wheel: length lcm(3*64, 5*64, 7*64) = 6720 bits = 105 u64 words.
const PRESIEVE_WORDS: usize = 3 * 5 * 7;
// Primes 2, 3, 5, 7 are not in the per-segment marking loop:
// 2 has no odd-bit slot, and 3/5/7 are handled by the presieve copy.
const SKIP_PRIMES: usize = 4;

fn isqrt(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u64;
    while x > 0 && x.saturating_mul(x) > n {
        x -= 1;
    }
    while (x + 1).saturating_mul(x + 1) <= n {
        x += 1;
    }
    x
}

// Bit g in odd-only space represents the odd number 2g+1.
// 2g+1 is divisible by 3 iff g % 3 == 1, by 5 iff g % 5 == 2, by 7 iff g % 7 == 3.
const fn build_presieve() -> [u64; PRESIEVE_WORDS] {
    let mut p = [0u64; PRESIEVE_WORDS];
    let mut g = 0;
    while g < PRESIEVE_WORDS * 64 {
        if g % 3 == 1 || g % 5 == 2 || g % 7 == 3 {
            p[g / 64] |= 1u64 << (g % 64);
        }
        g += 1;
    }
    p
}

const PRESIEVE: [u64; PRESIEVE_WORDS] = build_presieve();

// First mark of {offset, offset + step, offset + 2*step, ...} that lands at or past SEGSIZE.
fn next_mark_past_segment(offset: usize, step: usize) -> usize {
    offset + ((SEGSIZE - 1 - offset) / step + 1) * step
}

struct BoolArray {
    data: Vec<u64>,
}

impl BoolArray {
    fn new() -> Self {
        Self {
            data: vec![0; NUMBER_WORDS],
        }
    }

    // Marks bits at {offset, offset + step, ...} that fall within [0, SEGSIZE).
    // Returns the index of the first mark at or past SEGSIZE so callers can keep
    // a "next mark" cursor across segments without recomputing it.
    fn mark_multiples(&mut self, offset: usize, step: usize) -> usize {
        if offset >= SEGSIZE {
            return offset;
        }
        if step < 64 {
            self.mark_with_wheel(offset, step);
            next_mark_past_segment(offset, step)
        } else {
            self.mark_one_by_one(offset, step)
        }
    }

    // Bit-by-bit marking. Cheap for large strides.
    fn mark_one_by_one(&mut self, offset: usize, step: usize) -> usize {
        let mut i = offset;
        while i < SEGSIZE {
            self.data[i / 64] |= 1u64 << (i % 64);
            i += step;
        }
        i
    }

    // Build a step-word wheel encoding one period of the mark pattern, then OR it
    // into the segment cyclically. The first segment word gets masked so wraparound
    // bits below `offset % 64` (which belong to the previous period) are not set.
    fn mark_with_wheel(&mut self, offset: usize, step: usize) {
        let bit_off = offset % 64;
        let period_bits = step * 64;
        let mut wheel = vec![0u64; step];
        for k in 0..64 {
            let pos = (bit_off + k * step) % period_bits;
            wheel[pos / 64] |= 1u64 << (pos % 64);
        }

        let start = offset / 64;
        let first_mask = !((1u64 << bit_off) - 1);
        self.data[start] |= wheel[0] & first_mask;

        let mut j = 1 % step;
        for w in (start + 1)..NUMBER_WORDS {
            self.data[w] |= wheel[j];
            j += 1;
            if j == step {
                j = 0;
            }
        }
    }

    // Push the actual odd numbers that remain unmarked. `seg_lo` is the first odd
    // number represented (bit 0). Stops early once values would exceed u32::MAX.
    fn push_unmarked(&self, seg_lo: u64, out: &mut Vec<u32>) {
        for (i, &word) in self.data.iter().enumerate() {
            let mut bits = !word;
            let chunk_base = seg_lo + 2 * (i as u64) * 64;

            // Fast path: every value in this word fits in u32, so we skip per-bit checks.
            if chunk_base + 126 <= u32::MAX as u64 {
                while bits != 0 {
                    let b = bits.trailing_zeros() as u64;
                    out.push((chunk_base + 2 * b) as u32);
                    bits &= bits - 1;
                }
            } else {
                while bits != 0 {
                    let b = bits.trailing_zeros() as u64;
                    let val = chunk_base + 2 * b;
                    if val > u32::MAX as u64 {
                        return;
                    }
                    out.push(val as u32);
                    bits &= bits - 1;
                }
            }
        }
    }
}

// Number of post-bootstrap segments needed to cover [2*SEGSIZE+1, u32::MAX].
// Post-segment N covers odd numbers up to (2N+4)*SEGSIZE - 1.
const TOTAL_POST_SEGMENTS: usize = (1usize << 32) / (2 * SEGSIZE) - 1;

// Segments per work unit handed to a worker thread. Bigger = fewer activations
// per chunk (the only per-chunk division cost), smaller = better load balancing.
const CHUNK_SIZE: usize = 64;
const TOTAL_CHUNKS: usize = TOTAL_POST_SEGMENTS / CHUNK_SIZE
    + if TOTAL_POST_SEGMENTS % CHUNK_SIZE != 0 { 1 } else { 0 };

// Backpressure on workers: at most this many completed chunks waiting in the channel.
const CHANNEL_BOUND: usize = 16;

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::Arc;
use std::thread::{self, JoinHandle};

// Sieve a single segment. Used by both bootstrap and the worker loop.
fn sieve_one_segment(
    seg: &mut BoolArray,
    seg_lo: u64,
    presieve_offset: usize,
    small_primes: &[u32],
    next_bit: &mut Vec<usize>,
    out: &mut Vec<u32>,
) {
    // Presieve copy: lay down the 3-5-7 multiples pattern (also clears the segment).
    let mut wp = presieve_offset;
    for i in 0..NUMBER_WORDS {
        seg.data[i] = PRESIEVE[wp];
        wp += 1;
        if wp == PRESIEVE_WORDS {
            wp = 0;
        }
    }

    let seg_hi = seg_lo + 2 * SEGSIZE as u64;
    let sqrt_end = isqrt(seg_hi - 2);

    // Activate any new primes whose square has entered the sieving range.
    while next_bit.len() + SKIP_PRIMES < small_primes.len() {
        let p = small_primes[next_bit.len() + SKIP_PRIMES] as u64;
        if p > sqrt_end {
            break;
        }
        let lo = seg_lo.max(p * p);
        let mut m = ((lo + p - 1) / p) * p;
        if m % 2 == 0 {
            m += p;
        }
        let bit = ((m - seg_lo) / 2) as usize;
        next_bit.push(bit);
    }

    // Mark with active primes; cached cursor advances by mark_multiples' return.
    for i in 0..next_bit.len() {
        let p = small_primes[i + SKIP_PRIMES] as usize;
        let next = seg.mark_multiples(next_bit[i], p);
        next_bit[i] = next - SEGSIZE;
    }

    seg.push_unmarked(seg_lo, out);
}

// Run the bootstrap segment (odd numbers 1..2*SEGSIZE). Returns the primes found,
// starting with 2. Single-threaded -- workers depend on this finishing first.
fn run_bootstrap() -> Vec<u32> {
    let mut seg = BoolArray::new();
    let mut primes: Vec<u32> = vec![2];
    seg.data[0] |= 1u64; // bit 0 = number 1, not prime

    let largest = 2 * SEGSIZE - 1;
    let mut bit = 1usize;
    while (2 * bit + 1).pow(2) <= largest {
        let p = 2 * bit + 1;
        if seg.data[bit / 64] & (1u64 << (bit % 64)) == 0 {
            seg.mark_multiples((p * p - 1) / 2, p);
        }
        bit += 1;
    }
    seg.push_unmarked(1, &mut primes);
    primes
}

// Worker thread: pull chunks off the shared counter, sieve them, ship results.
// Each worker keeps its own segment buffer and next_bit cache, so the cache
// benefit is preserved within a chunk's contiguous segments.
fn worker_loop(
    small_primes: Arc<Vec<u32>>,
    counter: Arc<AtomicUsize>,
    stop: Arc<AtomicBool>,
    tx: SyncSender<(usize, Vec<u32>)>,
) {
    let mut seg = BoolArray::new();
    let mut next_bit: Vec<usize> = Vec::new();

    loop {
        if stop.load(Ordering::Relaxed) {
            break;
        }
        let chunk_idx = counter.fetch_add(1, Ordering::Relaxed);
        if chunk_idx >= TOTAL_CHUNKS {
            break;
        }

        let start_seg = chunk_idx * CHUNK_SIZE;
        let mut seg_lo = (2 + 2 * start_seg as u64) * SEGSIZE as u64 + 1;
        let mut presieve_offset = ((start_seg + 1) * NUMBER_WORDS) % PRESIEVE_WORDS;
        next_bit.clear();

        let mut primes = Vec::new();
        for _ in 0..CHUNK_SIZE {
            if seg_lo > u32::MAX as u64 {
                break;
            }
            sieve_one_segment(
                &mut seg,
                seg_lo,
                presieve_offset,
                &small_primes,
                &mut next_bit,
                &mut primes,
            );
            presieve_offset = (presieve_offset + NUMBER_WORDS) % PRESIEVE_WORDS;
            seg_lo += 2 * SEGSIZE as u64;
        }

        if tx.send((chunk_idx, primes)).is_err() {
            break;
        }
    }
}

struct PrimeStream {
    // Bootstrap primes are yielded first, in order.
    bootstrap: Vec<u32>,
    boot_cursor: usize,

    // Lazy worker spawn: only fire up threads if the consumer pulls past bootstrap.
    small_primes: Arc<Vec<u32>>,
    workers: Option<WorkerPool>,

    // Current in-order chunk being yielded, plus the cursor into it.
    current: Vec<u32>,
    cur_cursor: usize,
    next_chunk: usize,

    // Chunks that arrived ahead of `next_chunk` go here.
    pending: BTreeMap<usize, Vec<u32>>,
}

struct WorkerPool {
    rx: Receiver<(usize, Vec<u32>)>,
    stop: Arc<AtomicBool>,
    handles: Vec<JoinHandle<()>>,
}

impl PrimeStream {
    fn new() -> Self {
        let bootstrap = run_bootstrap();
        Self {
            small_primes: Arc::new(bootstrap.clone()),
            bootstrap,
            boot_cursor: 0,
            workers: None,
            current: Vec::new(),
            cur_cursor: 0,
            next_chunk: 0,
            pending: BTreeMap::new(),
        }
    }

    fn spawn_workers(&mut self) {
        let counter = Arc::new(AtomicUsize::new(0));
        let stop = Arc::new(AtomicBool::new(false));
        let (tx, rx) = sync_channel(CHANNEL_BOUND);

        let num_workers = thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
        let mut handles = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let sp = Arc::clone(&self.small_primes);
            let ct = Arc::clone(&counter);
            let st = Arc::clone(&stop);
            let tx = tx.clone();
            handles.push(thread::spawn(move || worker_loop(sp, ct, st, tx)));
        }
        drop(tx); // so rx sees Err once all workers exit

        self.workers = Some(WorkerPool { rx, stop, handles });
    }

    // Pull (and load) the next chunk in order, or None if all chunks are done.
    fn pull_next_chunk(&mut self) -> Option<Vec<u32>> {
        if self.next_chunk >= TOTAL_CHUNKS {
            return None;
        }
        if let Some(chunk) = self.pending.remove(&self.next_chunk) {
            self.next_chunk += 1;
            return Some(chunk);
        }
        let rx = &self.workers.as_ref().unwrap().rx;
        loop {
            match rx.recv() {
                Ok((idx, chunk)) => {
                    if idx == self.next_chunk {
                        self.next_chunk += 1;
                        return Some(chunk);
                    }
                    self.pending.insert(idx, chunk);
                }
                Err(_) => return None,
            }
        }
    }
}

impl Iterator for PrimeStream {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        // Phase 1: bootstrap primes.
        if self.boot_cursor < self.bootstrap.len() {
            let p = self.bootstrap[self.boot_cursor];
            self.boot_cursor += 1;
            return Some(p);
        }

        // Phase 2: chunked primes from worker threads. Spawn workers on first need.
        if self.workers.is_none() {
            self.spawn_workers();
        }

        loop {
            if self.cur_cursor < self.current.len() {
                let p = self.current[self.cur_cursor];
                self.cur_cursor += 1;
                return Some(p);
            }
            match self.pull_next_chunk() {
                Some(chunk) => {
                    self.current = chunk;
                    self.cur_cursor = 0;
                }
                None => return None,
            }
        }
    }
}

impl Drop for PrimeStream {
    fn drop(&mut self) {
        // Tell workers to stop at the next chunk boundary, then wait for them.
        if let Some(pool) = self.workers.take() {
            pool.stop.store(true, Ordering::Relaxed);
            drop(pool.rx); // closes channel; workers blocked on send will unblock with Err
            for h in pool.handles {
                let _ = h.join();
            }
        }
    }
}

fn stream() -> impl Iterator<Item = u32> {
    PrimeStream::new()
}

struct PrimeSieve;
impl PrimeSieve {
    fn up_to(limit: u64) -> Vec<u64> {
        stream()
            .take_while(|&p| p as u64 <= limit)
            .map(|p| p as u64)
            .collect()
    }
}

fn main() {
    let start = std::time::Instant::now();
    let p = stream()
        .nth(100_000_000)
        .expect("stream exhausted before the 100_000_001st prime");
    let end = std::time::Instant::now();
    println!("Time taken: {:.2?}", end - start);
    println!("100_000_001st prime: {}", p);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wheel_matches_simple() {
        for step in 1..64 {
            for offset in 0..200 {
                let mut a = BoolArray::new();
                a.mark_multiples(offset, step);

                let mut b = BoolArray::new();
                let mut i = offset;
                while i < SEGSIZE {
                    b.data[i / 64] |= 1u64 << (i % 64);
                    i += step;
                }
                assert_eq!(a.data, b.data, "offset={} step={}", offset, step);
            }
        }
    }

    #[test]
    fn sieve_small() {
        assert_eq!(
            PrimeSieve::up_to(30),
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        );
    }

    #[test]
    fn sieve_edge_cases() {
        assert_eq!(PrimeSieve::up_to(0), Vec::<u64>::new());
        assert_eq!(PrimeSieve::up_to(1), Vec::<u64>::new());
        assert_eq!(PrimeSieve::up_to(2), vec![2]);
        assert_eq!(PrimeSieve::up_to(3), vec![2, 3]);
        assert_eq!(PrimeSieve::up_to(4), vec![2, 3]);
        assert_eq!(PrimeSieve::up_to(5), vec![2, 3, 5]);
    }

    #[test]
    fn sieve_pi_million() {
        // pi(10^6) = 78498
        assert_eq!(PrimeSieve::up_to(1_000_000).len(), 78498);
    }

    #[test]
    fn sieve_pi_ten_million() {
        // pi(10^7) = 664579 -- exercises several segments past bootstrap
        assert_eq!(PrimeSieve::up_to(10_000_000).len(), 664579);
    }

    #[test]
    fn stream_first_10() {
        let first: Vec<u32> = stream().take(10).collect();
        assert_eq!(first, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn stream_lazy() {
        // Should not exhaust the iterator -- just confirm it yields lazily.
        let p100k = stream().nth(99_999).unwrap();
        assert_eq!(p100k, 1_299_709); // the 100000th prime
    }
}
