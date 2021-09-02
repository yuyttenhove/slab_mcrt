mod voronoi;
mod delaunay;
mod geometry;
mod voronoi_generator;

use crate::vector::Vec2;
use voronoi_generator::VoronoiGenerator;
use voronoi::VoronoiGrid;
use delaunay::DelaunayTriangulation;
use rand::{thread_rng, Rng};


pub struct VoronoiSlab {
    anchor: Vec2<f64>,
    sides: Vec2<f64>,
    generators: Vec<VoronoiGenerator>,
    vor_tess: VoronoiGrid,
}

impl VoronoiSlab {
    pub fn new(tau_max: f64, albedo: f64, g: f64, resolution: i64, perturbation: f64) -> VoronoiSlab {
        let anchor = Vec2::new(0., 0.);
        let sides = Vec2::new(1., 1.);
        let mut generators = Vec::new();

        let mut rng = thread_rng();
        for row in 0..resolution {
            for col in 0..resolution {
                let x = (col as f64 + 0.5) / resolution as f64 * (1.0 + perturbation * rng.gen_range(-0.5..0.5));
                let y = (row as f64 + 0.5) / resolution as f64 * (1.0 + perturbation * rng.gen_range(-0.5..0.5));
                generators.push(VoronoiGenerator::new(tau_max, albedo, g, x, y));
            }
        }

        let del_tess = DelaunayTriangulation::from_generators(&generators, anchor, sides, true);
        let vor_tess = VoronoiGrid::from_delaunay_triangulation(&del_tess);

        VoronoiSlab{anchor, sides, generators, vor_tess}
    }

    pub fn save_grid(&self, filename: &str) {
        self.vor_tess.to_file(filename);
    }
}

