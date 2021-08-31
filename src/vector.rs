use std::ops::{Add, Sub, Mul, Div};

#[derive(Copy, Clone)]
pub struct Vec2<T> where T: Copy {
    x: T,
    y: T,
}

impl<T> Vec2<T> where T: Copy {
    pub fn new(x: T, y: T) -> Vec2<T> {
        Vec2 {x, y}
    }

    pub fn x(&self) -> &T {
        &self.x
    }

    pub fn y(&self) -> &T {
        &self.y
    }
}

impl<T> Vec2<T> where T: Mul<Output = T> + Add<Output = T> + Copy {
    pub fn dot(&self, other: Self) -> T {
        self.x.mul(other.x).add(self.y.mul(other.y))
    }
}

impl<T> Add for Vec2<T> where T: Add<Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<T> Sub for Vec2<T> where T: Sub<Output = T> + Copy {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<T> Mul<T> for Vec2<T> where T: Mul<Output = T> + Copy {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Vec2::new(self.x * rhs, self.y * rhs)
    }
}

impl<T> Div<T> for Vec2<T> where T: Div<Output = T> + Copy {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output{
        Vec2::new(self.x / rhs, self.y / rhs)
    }
}