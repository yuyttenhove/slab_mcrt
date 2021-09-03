mod voronoi;
mod delaunay;
mod geometry;
mod voronoi_generator;

use crate::vector::Vec2;
use super::Slab;
use voronoi_generator::VoronoiGenerator;
use voronoi::VoronoiGrid;
use delaunay::DelaunayTriangulation;
use rand::{thread_rng, Rng};


#[derive(Copy, Clone, Debug)]
enum NeighbourDirection {
    Up,
    Right,
    Down,
    Left,
}

pub struct VoronoiSlab {
    anchor: Vec2<f64>,
    sides: Vec2<f64>,
    generators: Vec<VoronoiGenerator>,
    vor_tess: VoronoiGrid,
    resolution: i64,
}

impl VoronoiSlab {
    pub fn new(tau_max: f64, albedo: f64, g: f64, resolution: i64, perturbation: f64) -> VoronoiSlab {
        let anchor = Vec2::new(0., 0.);
        let sides = Vec2::new(1., 1.);
        let mut generators = Vec::new();

        let mut rng = thread_rng();
        for row in 0..resolution {
            for col in 0..resolution {
                let x = (col as f64 + 0.5) / resolution as f64 + (perturbation * rng.gen_range(-0.5..0.5) / resolution as f64);
                let y = (row as f64 + 0.5) / resolution as f64 + (perturbation * rng.gen_range(-0.5..0.5) / resolution as f64);
                generators.push(VoronoiGenerator::new(tau_max, albedo, g, x, y));
            }
        }

        let del_tess = DelaunayTriangulation::from_generators(&generators, anchor, sides, true);
        let vor_tess = VoronoiGrid::from_delaunay_triangulation(&del_tess);

        VoronoiSlab { anchor, sides, generators, vor_tess, resolution }
    }

    pub fn save_grid(&self, filename: &str) {
        self.vor_tess.to_file(filename);
    }
}

impl Slab for VoronoiSlab {
    fn raytrace_photon_packets(&self, number_of_packets: u64, number_of_detection_bins: usize) -> Vec<f64> {
        let mut result = vec![0.; number_of_detection_bins];
        let mut rng = rand::thread_rng();

        for _ in 0..number_of_packets {
            let mut x = 0.5;
            let mut y = 1.;
            let mut cos_theta = f64::sqrt(rng.gen_range(0.0..1.0));
            let mut weight = 1.0;

            /* Determine start cell */
            let mut cell_idx = 0;
            let mut min_dist = 1.;
            for i in (self.resolution * (self.resolution - 1)) as usize..(self.resolution * self.resolution) as usize {
                let cur_dist = (Vec2::new(x, y) - self.generators[i].position()).norm();
                if cur_dist < min_dist {
                    cell_idx = i;
                    min_dist = cur_dist;
                }
            }

            'outer:
            loop {
                // Determine sin_theta, only correct for isotropic scattering!
                let sign = if rng.gen_bool(0.5) { 1. } else { -1. };
                let sin_theta = sign * f64::sqrt((1. - cos_theta) * (1. + cos_theta));

                // determine a random path length, optionally using path length stretching
                let density = self.generators[cell_idx].density();
                let mut path_length_random = -f64::ln_1p(-rng.gen_range(0.0..1.0)) / density;

                let mut intersection_result= self.vor_tess.raytrace_cell(cell_idx, x, y, cos_theta, sin_theta, -1).unwrap();

                // Handle case of packet leaving cell or domain
                while path_length_random > intersection_result.distance {
                    // update location
                    x = intersection_result.intersection.x();
                    y = intersection_result.intersection.y();
                    path_length_random -= intersection_result.distance;

                    let exit_face = intersection_result.exit_face;
                    cell_idx = if exit_face.adjacent_cells[0] == cell_idx as i32 {
                        exit_face.adjacent_cells[1] as usize
                    } else {
                        exit_face.adjacent_cells[0] as usize
                    };

                    if let Some(neighbour_direction) = exit_face.neighbour_direction {
                        match neighbour_direction {
                            NeighbourDirection::Left => x += self.sides.x(),
                            NeighbourDirection::Right => x -= self.sides.x(),
                            NeighbourDirection::Up => break 'outer,
                            NeighbourDirection::Down => {
                                let bin = f64::floor(number_of_detection_bins as f64 * cos_theta) as usize;
                                result[bin] += weight;
                                break 'outer;
                            }
                        }
                    }

                    let previous_face = intersection_result.exit_face_idx;
                    match self.vor_tess.raytrace_cell(cell_idx, x, y, cos_theta, sin_theta, previous_face) {
                        Ok(res) => intersection_result = res,
                        Err(e) => {
                            println!("Encountered error! Skipping current photon packet... \n{}", e);
                            break 'outer;
                        }
                    }
                }

                // update location
                x += path_length_random * sin_theta;
                y -= path_length_random * cos_theta;

                // Correct the photon package weight for absorption
                weight *= self.generators[cell_idx].albedo();

                // generate a new propagation direction
                cos_theta = Self::scatter_photon(cos_theta, self.generators[cell_idx].g(), &mut rng);
            }
        }

        result
    }
}

