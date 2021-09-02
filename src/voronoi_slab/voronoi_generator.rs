use crate::vector::Vec2;

pub struct VoronoiGenerator {
    density: f64,
    albedo: f64,
    g: f64,
    position: Vec2<f64>,
}

impl VoronoiGenerator {
    pub fn density(&self) -> f64 {
        self.density
    }
    pub fn albedo(&self) -> f64 {
        self.albedo
    }
    pub fn g(&self) -> f64 {
        self.g
    }
    pub fn position(&self) -> Vec2<f64> {
        self.position
    }

    pub fn new(density: f64, albedo: f64, g: f64, x: f64, y: f64) -> VoronoiGenerator {
        VoronoiGenerator { density, albedo, g, position: Vec2::new(x, y) }
    }

    pub fn x(&self) -> f64 {
        self.position.x()
    }

    pub fn y(&self) -> f64 {
        self.position.y()
    }
}