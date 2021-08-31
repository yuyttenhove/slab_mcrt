use crate::slab::{GridLessSlab, Slab};

mod slab;

fn main() {
    let slab = GridLessSlab::new(2., 0.5, 0.);
    let result = slab.run_mcrt(1e9 as u64, 10, 500);
    print!("{{");
    for v in result {
        print!("{}, ", v);
    }
    println!("}}");
}
