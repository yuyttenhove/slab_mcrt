mod grid_cell;

use grid_cell::GridCell;
use crate::vector::Vec2;
use crate::slab::Slab;
use rand::Rng;

pub struct RegularGridSlab {
    anchor: Vec2<f64>,
    sides: Vec2<f64>,
    resolution: i64,
    cells: Vec<GridCell>,
}

impl RegularGridSlab {
    pub fn new(tau_max: f64, albedo: f64, g: f64, resolution: i64) -> RegularGridSlab {
        // We assume a box of size 1 in all dimensions, of homogeneous density and with kappa = 1.
        // Tau_max is then equal to z_max * rho = rho.
        let anchor = Vec2::new(0., 0.);
        let sides = Vec2::new(1., 1.);
        let cell_sides = sides / resolution as f64;

        let mut cells = Vec::new();
        for row in 0..resolution {
            for col in 0..resolution {
                cells.push(
                    GridCell::new(
                        tau_max,
                        albedo,
                        g,
                        anchor + Vec2::new(col as f64 * cell_sides.x(), row as f64 * cell_sides.y()),
                        cell_sides
                    )
                )
            }
        }
        RegularGridSlab {anchor, sides, resolution, cells}
    }

    fn get_cell_at(&self, row: i64, col: i64) -> &GridCell {
        &self.cells[row as usize * self.resolution as usize + col as usize]
    }

    fn path_length_cell(cell: &GridCell, cos_theta: f64, sin_theta: f64, x: f64, y: f64) -> (f64, f64) {
        let path_length_x = if sin_theta > 0. {
            (cell.anchor().x() + cell.sides().x() - x) / sin_theta
        } else {
            (cell.anchor().x() - x) / sin_theta
        };
        let path_length_y = if cos_theta > 0. {
            (y - cell.anchor().y()) / cos_theta
        } else {
            (y - cell.anchor().y() - cell.sides().y()) / cos_theta
        };
        (path_length_x, path_length_y)
    }
}


impl Slab for RegularGridSlab {
    fn raytrace_photon_packets(&self, number_of_packets: u64, number_of_detection_bins: usize) -> Vec<f64> {
        let mut result = vec![0.; number_of_detection_bins];
        let mut rng = rand::thread_rng();

        for _ in 0..number_of_packets {
            let mut x = rng.gen_range(0.0..1.0);
            let mut y = 1.0;
            let mut cos_theta = f64::sqrt(rng.gen_range(0.0..1.0));
            let mut weight = 1.0;

            /* Determine start cell */
            let mut row = self.resolution as i64 - 1;
            let mut col = f64::floor(x * self.resolution as f64) as i64;
            let mut cell = self.get_cell_at(row, col);

            'outer:
            loop {
                // Determine sin_theta, only correct for isotropic scattering!
                let sign = if rng.gen_bool(0.5) { 1. } else { -1. };
                let sin_theta = sign * f64::sqrt((1. - cos_theta) * (1. + cos_theta));

                // determine a random path length, optionally using path length stretching
                let density = cell.density();
                let mut path_length_random = -f64::ln_1p(-rng.gen_range(0.0..1.0)) / density;

                // calculate the path length within the current cell along the path
                let mut path_lengths_cell = Self::path_length_cell(cell, cos_theta, sin_theta, x, y);
                let mut path_length_cell = f64::min(path_lengths_cell.0, path_lengths_cell.1);

                // Handle case of packet leaving cell or domain
                while path_length_random > path_length_cell {
                    // update location
                    x += path_length_cell * sin_theta;
                    y -= path_length_cell * cos_theta;
                    path_length_random -= path_length_cell;
                    if path_lengths_cell.0 < path_lengths_cell.1 {
                        col = if sin_theta > 0. { col + 1 } else { col - 1 };
                        if col >= self.resolution {
                            col -= self.resolution;
                            x -= self.sides.x();
                        } else if col < 0 {
                            col += self.resolution;
                            x += self.sides.x();
                        }
                    } else {
                        row = if cos_theta < 0. { row + 1 } else { row - 1 };
                        if row >= self.resolution {
                            break 'outer;
                        } else if row < 0 {
                            let bin = f64::floor(number_of_detection_bins as f64 * cos_theta) as usize;
                            result[bin] += weight;
                            break 'outer;
                        }
                    }
                    cell = self.get_cell_at(row, col);
                    path_lengths_cell = Self::path_length_cell(cell, cos_theta, sin_theta, x, y);
                    path_length_cell = f64::min(path_lengths_cell.0, path_lengths_cell.1);
                }

                // update location
                x += path_length_random * sin_theta;
                y -= path_length_random * cos_theta;

                // Correct the photon package weight for absorption
                weight *= cell.albedo();

                // generate a new propagation direction
                cos_theta = Self::scatter_photon(cos_theta, cell.g(), &mut rng);
            }
        }

        result
    }
}