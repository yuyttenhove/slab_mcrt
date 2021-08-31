use super::Slab;
use rand::Rng;

pub struct GridLessSlab {
    tau_max: f64,
    albedo: f64,
    g: f64,
}

impl GridLessSlab {
    pub fn new(tau_max: f64, albedo: f64, g: f64) -> GridLessSlab {
        GridLessSlab { tau_max, albedo, g }
    }
}

impl Slab for GridLessSlab {
    fn raytrace_photon_packets(&self, number_of_packets: u64, number_of_detection_bins: usize) -> Vec<f64> {
        let mut result = vec![0.; number_of_detection_bins];
        let mut rng = rand::thread_rng();

        for _ in 0..number_of_packets {
            let mut tau = self.tau_max;
            let mut mu = f64::sqrt(rng.gen_range(0.0..1.0));
            let mut weight = 1.0;

            loop {
                // calculate the optical depth along the path
                let tau_path = if mu > 0. { tau / mu } else { (tau - self.tau_max) / mu };

                // determine a random optical depth, optionally using path length stretching
                let tau_random = -f64::ln_1p(-rng.gen_range(0.0..1.0));

                // if the photon package leaves the slab, terminate its life cycle
                if tau_random > tau_path {
                    // if the photon package emerges at the bottom, detect it in the appropriate bin
                    if mu > 0. {
                        let bin = f64::floor(number_of_detection_bins as f64 * mu) as usize;
                        result[bin] += weight;
                    }
                    break;
                }

                // Correct the photon package weight for absorption
                weight *= self.albedo;

                // calculate the interaction location
                tau -= tau_random * mu;

                // generate a new propagation direction
                mu = Self::scatter_photon(mu, self.g, &mut rng);
            }
        }
        result
    }
}