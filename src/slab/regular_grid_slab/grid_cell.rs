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
}