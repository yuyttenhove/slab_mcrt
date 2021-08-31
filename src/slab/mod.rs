mod gridless_slab;
mod voronoi_slab;
mod regular_grid_slab;

use rayon::prelude::{ParallelIterator, IntoParallelIterator};

// Re-exports
pub use gridless_slab::GridLessSlab;

pub trait Slab: Sync {
    fn raytrace_photon_packets(&self, number_of_packets: u64, number_of_detection_bins: usize) -> Vec<f64>;

    fn run_mcrt(&self, number_of_packets: u64, number_of_detection_bins: usize, number_of_jobs: u64) -> Vec<f64> {
        let number_of_packets_per_job = number_of_packets / number_of_jobs;
        let all_results: Vec<Vec<f64>> = (0..number_of_jobs).into_par_iter().map(|_| {
            self.raytrace_photon_packets(number_of_packets_per_job, number_of_detection_bins)
        }).collect();

        let mut result = vec![0.; number_of_detection_bins];
        for thread_result in all_results {
            for (i, v) in thread_result.iter().enumerate() {
                result[i] += v;
            }
        }
        Self::convert_photon_count_to_intensity(result, number_of_packets)
    }

    fn convert_photon_count_to_intensity(counts: Vec<f64>, number_of_packets: u64) -> Vec<f64> {
        let n_bins = counts.len();
        let mut intensities = vec![0.0; n_bins];
        for j in 0..n_bins {
            let mu = (j as f64 + 0.5) / n_bins as f64;
            intensities[j] = counts[j] * (n_bins as f64) / (number_of_packets as f64 * mu * 2.);
        }
        intensities
    }
}
