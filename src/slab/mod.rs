mod gridless_slab;
mod voronoi_slab;
mod regular_grid_slab;

use rayon::prelude::{ParallelIterator, IntoParallelIterator};
use indicatif::{ProgressBar, ParallelProgressIterator, ProgressStyle};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

// Re-exports
pub use gridless_slab::GridLessSlab;

pub trait Slab: Sync {
    fn raytrace_photon_packets(&self, number_of_packets: u64, number_of_detection_bins: usize) -> Vec<f64>;

    fn run_mcrt(&self, number_of_packets: u64, number_of_detection_bins: usize, number_of_jobs: u64) -> Vec<f64> {
        let number_of_packets_per_job = number_of_packets / number_of_jobs;
        let pb = ProgressBar::new(number_of_jobs);
        pb.set_style(ProgressStyle::default_bar().template("[{elapsed_precise}] [{bar:40}] [{percent}%]").progress_chars("##-"));
        pb.enable_steady_tick(1000);
        pb.println("Launching photon packets!");
        let all_results: Vec<Vec<f64>> = (0..number_of_jobs).into_par_iter().progress_with(pb).map(|_| {
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

    fn save_intensities_to_file(intensities: Vec<f64>, filename: &str) -> std::io::Result<()> {
        // Create a path to the desired file
        let path = Path::new(filename);
        println!("Saving intensities to {}", path.display());
        let mut file = File::create(&path)?;
        file.write(b"#mu\tI\n")?;
        let n_bins = intensities.len();
        for (bin, intensity) in intensities.iter().enumerate() {
            let mu = (bin as f64 + 0.5) / n_bins as f64;
            let line = format!("{}\t{}\n", mu, intensity);
            file.write(line.as_ref())?;
        }
        Ok(())
    }
}
