///
/// math.rs
/// Contains definitons and functions for the following structs:
/// Point
/// Vector
/// Matrix4D
/// CubicSplineInterpolation
///

use std::f64;

//
// Point type
//
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    //w: f64,
}

impl Point {
    pub fn new_blank() -> Point {
        Point {x: 0., y: 0., z: 0.}
    }
    pub fn new(xn: f64,yn: f64,zn: f64) -> Point {
        Point {x: xn, y: yn, z: zn}
    }
    pub fn set_points(&mut self, xn: f64,yn: f64,zn: f64) {
        self.x=xn;
        self.y=yn;
        self.z=zn;
    }
}

//
// Vector type
//
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Vector {
    pub point_a: Point,
    pub point_b: Point,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new_blank() -> Vector {
        Vector {point_a: Point::new_blank(), point_b: Point::new_blank(), x: 1., y: 1., z: 1.}
    }
    pub fn new(a: Point, b: Point, xn: f64,yn: f64,zn: f64) -> Vector {
        Vector {point_a: a, point_b: b, x: xn, y: yn, z: zn}
    }
    pub fn dot_product(&self, other: &Vector) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn length(&self) -> f64 {
        f64::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }
    pub fn cross_product(&self, other: &Vector) -> Vector {
        let xn = self.y * other.z - self.z * other.y;
        let yn = self.z * other.x - self.x * other.z;
        let zn = self.x * other.y - self.y * other.x;
        let pt_a = self.point_a.clone();
        let pt_b = Point::new(pt_a.x + xn, pt_a.y + yn, pt_a.z + zn);
        Vector::new(pt_a, pt_b, xn, yn , zn)
    }
    pub fn normalize(&mut self) {
        let length = self.length();
        self.x = self.x/length;
        self.y = self.y/length;
        self.z = self.z/length;
        self.point_b.set_points(
            self.point_a.x+self.x,
            self.point_a.y+self.y,
            self.point_a.z+self.z);
    }
    pub fn multiply_scalar(&mut self, scalar: f64) {
        self.x = self.x*scalar;
        self.y = self.y*scalar;
        self.z = self.z*scalar;
        self.point_b.set_points(
            self.point_a.x+self.x,
            self.point_a.y+self.y,
            self.point_a.z+self.z);
    }
    pub fn angle(&mut self, other: &mut Vector) -> f64 {
        self.normalize();
        other.normalize();

        f64::acos(self.dot_product(other)) //Check if this is correct. TODO
    }
    pub fn angle2(&self, other: &Vector) -> f64 {
        let l1 = self.length();
        let l2 = other.length();
        let ca = (self.dot_product(other)) / ( l1 * l2 );

        f64::acos(ca) //Check if this is correct. TODO
    }
}

//
// CubicSplineInterpolation type
//
#[derive(PartialEq, Clone, Debug)]
pub struct CubicSplineInterpolation {
    sd: Box<[f64]>, //second derivative has been abberviated because it appears in equations...
    //...And a long variable name is distracting.
    x: Box<[f64]>,
    y: Box<[f64]>,
}

impl CubicSplineInterpolation {
    //This did not make sense to me as a non-static function in the original code
    pub fn create(x: Box<[f64]>, y: Box<[f64]>) -> CubicSplineInterpolation {
        let mut u = vec![0.;x.len()];
        let mut d2u = vec![0.;x.len()]; //second derivative

        for i in 1..x.len() {
            let sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
            let p = sig * d2u[i - 1] + 2.;
            d2u[i] = (sig - 1.) / p;
            u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
            u[i] = (6. * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
        }

        let qn = 0.;
        let un = 0.;

        d2u[x.len()-1] =  (un - qn * u[x.len() - 1]) / (qn * d2u[x.len() - 1] + 1.);

        for i in (0..x.len()-2).rev() {
            d2u[i] = d2u[i] * d2u[i+1] + u[i];
        }

        CubicSplineInterpolation {sd: d2u.into_boxed_slice(), x: x, y: y}
    }

    pub fn interpolate(&self, x: f64) -> f64 {
        let mut klo = 0;
        let mut khi = self.x.len();

        while khi - klo > 1 {
            let k = (khi + klo) >> 1;
            if self.x[k]>x {
                khi = k;
            }
            else {
                klo = k;
            }
        }

        let h = self.x[khi] - self.x[klo];
        assert!(h == 0.,"Interpolation Error");
        let a = (self.x[khi] - x) / h;
        let b = (x - self.x[klo]) / h;

        a * self.y[klo] + b * self.y[khi]
        + ((a * a * a - a) * self.sd[klo]
        + (b * b * b - b) * self.sd[khi])
        * (h * h) / 6.
    }
}

//
// Matrix4D type
//
const ROWS: usize = 4;
const COLS: usize = 4;
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Matrix4D {
    pub data: [[f64; COLS]; ROWS],
}

impl Matrix4D {
    pub fn new() -> Matrix4D {
        let mut m4d = Matrix4D {data: [[0f64; COLS]; ROWS]};
        m4d.identity();
        m4d
    }
    // pub fn new_from_data(data: [f64; COLS*ROWS]) {
    //     //Check if this is used anywhere
    // }
    pub fn identity(&mut self) {
        let mut i = 0;
        let mut j = 0;
        while i<COLS && j<ROWS {
            self.data[i][j]=1.;
            i+=1;
            j+=1;
        }
    }
    pub fn x_rotation(angle: f64) -> Matrix4D {
        let ca = f64::cos(angle);
        let sa = f64::sin(angle);

        let r1 = [1., 0.,  0.,  0.];
        let r2 = [0., ca, sa, 0.];
        let r3 = [0., -sa, ca, 0.];
        let r4 = [0., 0., 0., 1.];

        Matrix4D {data: [r1,r2,r3,r4]}
    }
    pub fn y_rotation(angle: f64) -> Matrix4D {
        let ca = f64::cos(angle);
        let sa = f64::sin(angle);

        let r1 = [ca, 0., -sa,  0.];
        let r2 = [0., 1., sa, 0.];
        let r3 = [sa, 0., ca, 0.];
        let r4 = [0., 0., 0., 1.];

        Matrix4D {data: [r1,r2,r3,r4]}
    }
    pub fn z_rotation(angle: f64) -> Matrix4D {
        let ca = f64::cos(angle);
        let sa = f64::sin(angle);

        let r1 = [ca, sa, 0., 0.];
        let r2 = [-sa, ca, 0., 0.];
        let r3 = [0., 0., 1., 0.];
        let r4 = [0., 0., 0., 1.];

        Matrix4D {data: [r1,r2,r3,r4]}
    }
    pub fn translation(dx: f64, dy: f64, dz: f64) -> Matrix4D {
        let r1 = [1., 0., 0., dx];
        let r2 = [0., 1., 0., dy];
        let r3 = [0., 0., 1., dz];
        let r4 = [0., 0., 0., 1.];

        Matrix4D {data: [r1,r2,r3,r4]}
    }
    pub fn scaling(sx: f64, sy: f64, sz: f64) -> Matrix4D {
        let r1 = [sx, 0., 0., 0.];
        let r2 = [0., sy, 0., 0.];
        let r3 = [0., 0., sz, 0.];
        let r4 = [0., 0., 0., 1.];

        Matrix4D {data: [r1,r2,r3,r4]}
    }
    pub fn multiply(&self, other: &Matrix4D) -> Matrix4D {
        let mut m = Matrix4D {data: [[0f64; COLS]; ROWS]};
        // let mut i = 0;
        // let mut j = 0;
        // let mut k = 0;
        for i in 0..ROWS {
            for j in 0..COLS {
                let mut sum = 0.;

                for k in 0..ROWS {
                    sum += self.data[i][k] * other.data[k][i];
                }

                m.data[i][j]=sum;
            }
        }
        return m; //The "m" alone is hard to see
    }
    pub fn multiply_vector(&self, vector: &Vector) -> Vector {
        //TODO Something about this function seems wrong
        //The vector is only assigned the three scalars - and is left with the two default points
        //Is this intended?
        let mut v = Vector::new_blank();

        v.x = self.data[0][0] * vector.x + self.data[0][1] * vector.y + self.data[0][2] * vector.z;
        v.y = self.data[1][0] * vector.x + self.data[1][1] * vector.y + self.data[1][2] * vector.z;
        v.z = self.data[2][0] * vector.x + self.data[2][1] * vector.y + self.data[2][2] * vector.z;

        return v; //The "m" alone is hard to see
    }
    pub fn multiply_point(&self, point: Point) -> Point {
        let x = self.data[0][0] * point.x + self.data[0][1] * point.y + self.data[0][2] * point.z + self.data[0][3] * 1.;// point.w;
        let y = self.data[1][0] * point.x + self.data[1][1] * point.y + self.data[1][2] * point.z + self.data[1][3] * 1.;// point.w;
        let z = self.data[2][0] * point.x + self.data[2][1] * point.y + self.data[2][2] * point.z + self.data[2][3] * 1.;// point.w;
        //let w = self.data[3][0] * point.x + self.data[3][1] * point.y + self.data[3][2] * point.z + self.data[3][3] * point.w;

        Point {x: x, y: y, z: z}
    }
    pub fn transpose(&mut self) {
        //TODO do not use
        panic!();
    }
    pub fn to_array(&self) {
        //TODO do not use
        panic!();
    }
}
