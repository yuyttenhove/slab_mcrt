use std::ops::{Add, Sub, Mul, Div, AddAssign, DivAssign};
use num_traits::Float;

#[derive(Copy, Clone, Debug, Default)]
pub struct Vec2<T> where T: Copy {
    x: T,
    y: T,
}

impl<T> Vec2<T> where T: Copy {
    pub fn new(x: T, y: T) -> Vec2<T> {
        Vec2 {x, y}
    }

    pub fn x(&self) -> T {
        self.x
    }

    pub fn y(&self) -> T {
        self.y
    }
}

impl<T> Vec2<T> where T: Mul<Output = T> + Add<Output = T> + Copy {
    pub fn dot(&self, other: &Self) -> T {
        self.x.mul(other.x).add(self.y.mul(other.y))
    }

    pub fn norm2(&self) -> T {
        self.dot(&self)
    }
}

impl<T> Vec2<T> where T: Mul<Output = T> + Add<Output = T> + Copy + Float {
    pub fn norm(&self) -> T {
        self.norm2().sqrt()
    }
}

impl<T> Add for Vec2<T> where T: Add<Output = T> + Copy {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<T> AddAssign for Vec2<T> where T: AddAssign + Copy {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
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

impl<T> DivAssign<T> for Vec2<T> where T: DivAssign + Copy {
    fn div_assign(&mut self, rhs: T) {
        self.x /= rhs;
        self.y /= rhs;
    }
}