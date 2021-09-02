#[derive(Debug, Default, Copy, Clone)]
pub struct SimulationDomain {
    anchor: [f64; 2],
    sides: [f64; 2],
}

impl SimulationDomain {
    pub fn new(a: [f64; 2], s: [f64; 2]) -> SimulationDomain {
        SimulationDomain {anchor: a, sides: s}
    }

    pub fn anchor(&self) -> [f64; 2] {
        self.anchor.clone()
    }

    pub fn sides(&self) -> [f64; 2] {
        self.sides.clone()
    }
}