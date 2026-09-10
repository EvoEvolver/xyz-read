use std::ops::{Add, AddAssign, Div, Mul, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    pub const fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    pub fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn length(self) -> f32 {
        self.dot(self).sqrt()
    }

    pub fn normalized(self) -> Self {
        let length = self.length();
        if length == 0.0 {
            Self::default()
        } else {
            self / length
        }
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f32> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Div<f32> for Vec3 {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

pub fn rotate_euler(point: Vec3, x_degrees: f32, y_degrees: f32, z_degrees: f32) -> Vec3 {
    let (sin_x, cos_x) = x_degrees.to_radians().sin_cos();
    let (sin_y, cos_y) = y_degrees.to_radians().sin_cos();
    let (sin_z, cos_z) = z_degrees.to_radians().sin_cos();

    let around_x = Vec3::new(
        point.x,
        point.y * cos_x - point.z * sin_x,
        point.y * sin_x + point.z * cos_x,
    );
    let around_y = Vec3::new(
        around_x.x * cos_y + around_x.z * sin_y,
        around_x.y,
        -around_x.x * sin_y + around_x.z * cos_y,
    );
    Vec3::new(
        around_y.x * cos_z - around_y.y * sin_z,
        around_y.x * sin_z + around_y.y * cos_z,
        around_y.z,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotates_around_z() {
        let rotated = rotate_euler(Vec3::new(1.0, 0.0, 0.0), 0.0, 0.0, 90.0);
        assert!(rotated.x.abs() < 0.0001);
        assert!((rotated.y - 1.0).abs() < 0.0001);
    }
}
