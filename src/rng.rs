use rand::{
    Rng, RngCore, SeedableRng,
    distr::{
        Distribution, StandardUniform,
        uniform::{SampleRange, SampleUniform},
    },
    rngs::SmallRng,
};
use std::cell::RefCell;

thread_local! {
    static THREAD_RNG: RefCell<SmallRng> = RefCell::new(SmallRng::from_rng(&mut rand::rng()));
}

#[derive(Clone, Copy, Debug)]
pub struct FastRng;

pub fn fast_rng() -> FastRng {
    FastRng
}

pub fn fast_random<T>() -> T
where
    StandardUniform: Distribution<T>,
{
    fast_rng().random()
}

pub fn fast_random_range<T, R>(range: R) -> T
where
    T: SampleUniform,
    R: SampleRange<T>,
{
    fast_rng().random_range(range)
}

impl RngCore for FastRng {
    fn next_u32(&mut self) -> u32 {
        THREAD_RNG.with(|rng| rng.borrow_mut().next_u32())
    }

    fn next_u64(&mut self) -> u64 {
        THREAD_RNG.with(|rng| rng.borrow_mut().next_u64())
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        THREAD_RNG.with(|rng| rng.borrow_mut().fill_bytes(dest))
    }
}
