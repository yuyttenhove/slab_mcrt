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
                let taupath = if mu > 0. { tau / mu } else { (tau - self.tau_max) / mu };

                // determine a random optical depth, optionally using path length stretching
                let taurandom = -f64::ln_1p(-rng.gen_range(0.0..1.0));

                // if the photon package leaves the slab, terminate its life cycle
                if taurandom > taupath {
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
                tau -= taurandom * mu;

                // generate a new propagation direction
                if f64::abs(self.g) < 1e-6 {
                    // isotropic scattering
                    mu = rng.gen_range(-1.0..1.0);
                } else {
                    // anisotropic scattering; use 3D Henyey-Greenstein phase function
                    let f = ((1.0 - self.g) * (1.0 + self.g)) / (1.0 - self.g + self.g * rng.gen_range(0.0..2.0));
                    let cos_theta = (1.0 + self.g * self.g - f * f) / (2.0 * self.g);
                    let sin_theta = f64::sqrt(f64::abs((1.0 - cos_theta) * (1.0 + cos_theta)));
                    let cos_phi = f64::cos(std::f64::consts::PI * rng.gen_range(0.0..2.0));
                    let nu = f64::sqrt(f64::abs((1.0 - mu) * (1.0 + mu)));
                    mu = nu * sin_theta * cos_phi + mu * cos_theta;
                }
            }
        }
        result
    }
}