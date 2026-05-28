// Single-threaded u64 variant of the segmented sieve. Cap on iteration is the
// square of the largest bootstrap prime (~17.18B with SEGSIZE = 1 << 16,
// i.e. ~720M primes). Bump SEGSIZE for more range.

const SEGSIZE: usize = 1 << 16;
const NUMBER_WORDS: usize = SEGSIZE / 64;

const PRESIEVE_WORDS: usize = 3 * 5 * 7;
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

    fn mark_one_by_one(&mut self, offset: usize, step: usize) -> usize {
        let mut i = offset;
        while i < SEGSIZE {
            self.data[i / 64] |= 1u64 << (i % 64);
            i += step;
        }
        i
    }

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

    fn push_unmarked(&self, seg_lo: u64, out: &mut Vec<u64>) {
        for (i, &word) in self.data.iter().enumerate() {
            let mut bits = !word;
            let chunk_base = seg_lo + 2 * (i as u64) * 64;
            while bits != 0 {
                let b = bits.trailing_zeros() as u64;
                out.push(chunk_base + 2 * b);
                bits &= bits - 1;
            }
        }
    }
}

struct PrimeStream {
    primes: Vec<u64>,
    cursor: usize,
    seg: BoolArray,
    seg_lo: u64,
    done: bool,
    next_bit: Vec<usize>,
    presieve_offset: usize,
    // max_sieve = last_bootstrap_prime^2. Numbers up to this are correctly sieved.
    max_sieve: u64,
}

impl PrimeStream {
    fn new() -> Self {
        let mut s = Self {
            primes: Vec::new(),
            cursor: 0,
            seg: BoolArray::new(),
            seg_lo: 2 * SEGSIZE as u64 + 1,
            done: false,
            next_bit: Vec::new(),
            presieve_offset: 0,
            max_sieve: 0,
        };
        s.bootstrap();
        let last = *s.primes.last().expect("bootstrap finds at least 2");
        s.max_sieve = last * last;
        s
    }

    fn bootstrap(&mut self) {
        self.primes.push(2);
        self.seg.data[0] |= 1u64;

        let largest = 2 * SEGSIZE - 1;
        let mut bit = 1usize;
        while (2 * bit + 1).pow(2) <= largest {
            let p = 2 * bit + 1;
            if self.seg.data[bit / 64] & (1u64 << (bit % 64)) == 0 {
                self.seg.mark_multiples((p * p - 1) / 2, p);
            }
            bit += 1;
        }
        self.seg.push_unmarked(1, &mut self.primes);

        self.presieve_offset = NUMBER_WORDS % PRESIEVE_WORDS;
    }

    fn sieve_next_segment(&mut self) {
        // Stop if the segment's largest number would go past what we can correctly sieve.
        if self.seg_lo + 2 * SEGSIZE as u64 - 1 > self.max_sieve {
            self.done = true;
            return;
        }

        self.lay_down_presieve();

        let seg_hi = self.seg_lo + 2 * SEGSIZE as u64;
        let sqrt_end = isqrt(seg_hi - 2);

        self.activate_new_primes(sqrt_end);
        self.mark_with_active_primes();

        self.seg.push_unmarked(self.seg_lo, &mut self.primes);
        self.seg_lo = seg_hi;
    }

    fn lay_down_presieve(&mut self) {
        let mut wp = self.presieve_offset;
        for i in 0..NUMBER_WORDS {
            self.seg.data[i] = PRESIEVE[wp];
            wp += 1;
            if wp == PRESIEVE_WORDS {
                wp = 0;
            }
        }
        self.presieve_offset = (self.presieve_offset + NUMBER_WORDS) % PRESIEVE_WORDS;
    }

    fn activate_new_primes(&mut self, sqrt_end: u64) {
        while self.next_bit.len() + SKIP_PRIMES < self.primes.len() {
            let p = self.primes[self.next_bit.len() + SKIP_PRIMES];
            if p > sqrt_end {
                break;
            }
            let lo = self.seg_lo.max(p * p);
            let mut m = ((lo + p - 1) / p) * p;
            if m % 2 == 0 {
                m += p;
            }
            let bit = ((m - self.seg_lo) / 2) as usize;
            self.next_bit.push(bit);
        }
    }

    fn mark_with_active_primes(&mut self) {
        for i in 0..self.next_bit.len() {
            let p = self.primes[i + SKIP_PRIMES] as usize;
            let next = self.seg.mark_multiples(self.next_bit[i], p);
            self.next_bit[i] = next - SEGSIZE;
        }
    }
}

impl Iterator for PrimeStream {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        loop {
            if self.cursor < self.primes.len() {
                let p = self.primes[self.cursor];
                self.cursor += 1;
                return Some(p);
            }
            if self.done {
                return None;
            }
            self.sieve_next_segment();
        }
    }
}

fn stream() -> impl Iterator<Item = u64> {
    PrimeStream::new()
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
