use crate::vector::Vec2;

pub struct GridCell {
    density: f64,
    albedo: f64,
    g: f64,
    anchor: Vec2<f64>,
    sides: Vec2<f64>
}

impl GridCell {
    pub fn new(density: f64, albedo: f64, g: f64, anchor: Vec2<f64>, sides: Vec2<f64>) -> GridCell {
        GridCell {density, albedo, g, anchor, sides}
    }

    pub fn anchor(&self) -> &Vec2<f64> {
        &self.anchor
    }

    pub fn sides(&self) -> &Vec2<f64> {
        &self.sides
    }

    pub fn density(&self) -> f64 {
        self.density
    }

    pub fn albedo(&self) -> f64 {
        self.albedo
    }

    pub fn g(&self) -> f64 {
        self.g
    }
}