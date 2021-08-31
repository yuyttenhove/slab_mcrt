mod grid_cell;

use grid_cell::GridCell;
use crate::vector::Vec2;

pub struct RegularGridSlab {
    anchor: Vec2<f64>,
    sides: Vec2<f64>,
    resolution: u64,
    cells: Vec<GridCell>,
}

impl RegularGridSlab {
    pub fn new(tau_max: f64, albedo: f64, g: f64, resolution: u64) -> RegularGridSlab {
        // We assume a box of size 1 in all dimensions, of homogeneous density and with kappa = 1.
        // Tau_max is then equal to z_max * rho = rho.
        let anchor = Vec2::new(0., 0.);
        let sides = Vec2::new(0., 0.);
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
}