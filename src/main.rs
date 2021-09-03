use crate::slab::{GridLessSlab, RegularGridSlab, Slab};
use crate::voronoi_slab::VoronoiSlab;

mod slab;
mod vector;
mod gridless_slab;
mod regular_grid_slab;
mod voronoi_slab;

fn main() {
    let tau_max = 10.0;
    let albedo = 0.5;
    let g = 0.0;
    let n_bins = 10;
    let log_n_photons = 4;

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

    let regular_grid_slab = RegularGridSlab::new(tau_max, albedo, g, 10);
    RegularGridSlab::save_intensities_to_file(
        regular_grid_slab.run_mcrt(u64::pow(10, log_n_photons), n_bins, 500),
        &format!(
            "output/regular_t{}_a{}_g{}_n{}_p{}.txt",
            tau_max as u64,
            10 * albedo as u64,
            10 * g as u64,
            n_bins,
            log_n_photons
        )
    ).unwrap();

    let voronoi_slab = VoronoiSlab::new(tau_max, albedo, g, 10, 0.5);
    voronoi_slab.save_grid("vor_test.txt");
}
