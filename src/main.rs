use crate::slab::{GridLessSlab, Slab};

mod slab;
mod vector;

fn main() {
    let tau_max = 10.0;
    let albedo = 0.5;
    let g = 0.0;
    let n_bins = 10;
    let log_n_photons = 9;

    let slab = GridLessSlab::new(tau_max, albedo, g);
    GridLessSlab::save_intensities_to_file(
        slab.run_mcrt(u64::pow(10, log_n_photons), n_bins, 500),
        &format!(
            "output/gridless_t{}_a{}_g{}_n{}_p{}.txt",
            tau_max as u64,
            10 * albedo as u64,
            10 * g as u64,
            n_bins,
            log_n_photons
        )
    ).unwrap();
}
