///
/// math.rs
/// Contains definitons and functions for the following structs:
/// Point
/// Vector
/// Matrix4D
/// CubicSplineInterpolation
///

extern crate nalgebra;

use std::f64;
use nalgebra::{Inverse,DMatrix3};

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
    pub p: Point,
    // pub p_b: Point,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new_blank() -> Vector {
        Vector {p: Point::new_blank(), x: 1., y: 1., z: 1.}
    }
    pub fn new_from_points(a: Point, b: Point) -> Vector {
        Vector {p: a, x: b.x-a.x, y: b.y-a.y, z: b.z-a.z}
    }
    pub fn new(a: Point, xn: f64,yn: f64,zn: f64) -> Vector {
        Vector {p: a, x: xn, y: yn, z: zn}
    }
    pub fn dot_product(&self, other: &Vector) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross_product(&self, other: &Vector) -> Vector {
        let xn = self.y * other.z - self.z * other.y;
        let yn = self.z * other.x - self.x * other.z;
        let zn = self.x * other.y - self.y * other.x;
        let pt_a = self.p.clone();
        //let pt_b = Point::new(pt_a.x + xn, pt_a.y + yn, pt_a.z + zn);
        Vector::new(pt_a, xn, yn , zn)
    }
    pub fn length(&self) -> f64 {
        f64::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }
    pub fn normalize(&mut self) {
        let length = self.length();
        self.x = self.x/length;
        self.y = self.y/length;
        self.z = self.z/length;
        // self.p_b.set_points(
        //     self.p.x+self.x,
        //     self.p.y+self.y,
        //     self.p.z+self.z);
    }
    pub fn intersect(&self,other: &Vector) -> Point {
        //Assumes intersection happens, does not check if actually does.
        let mix = Vector::new_from_points(self.p, other.p);
        let left = self.cross_product(&other);
        let right = mix.cross_product(&other);
        // let t = right.length()/left.length();
        let t = right.dot_product(&left)/(left.length().powi(2));
        Point {x: self.p.x+self.x*t, y: self.p.y+self.y*t, z: self.p.z+self.z*t}
    }
    pub fn reflection(&self, normal: &Vector, point: Point) -> Vector {
        let dot=self.dot_product(&normal);
        let x = self.x-2.*dot*normal.x;
        let y = self.y-2.*dot*normal.y;
        let z = self.z-2.*dot*normal.z;
        let mut v = Vector {p: point, x: x, y: y, z: z};
        v.normalize();
        return v;
    }
    pub fn point_at_t_equals(&self, t: f64) -> Point {
        Point{x: self.p.x+self.x*t, y: self.p.y+self.y*t, z: self.p.z+self.z*t}
    }
    pub fn multiply_scalar(&mut self, scalar: f64) {
        self.x = self.x*scalar;
        self.y = self.y*scalar;
        self.z = self.z*scalar;
        // self.p_b.set_points(
        //     self.p.x+self.x,
        //     self.p.y+self.y,
        //     self.p.z+self.z);
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
        unimplemented!()
    }
    pub fn to_array(&self) {
        //TODO do not use
        unimplemented!()
    }
}


//
// Plane type
// Represented as a+(b-a)*n+(c-a)*m
//
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Plane {
    pub a: Point,
    pub b: Point,
    pub c: Point,
}

impl Plane {
    //An intersection would be v.p+v(x,y,z)*t=a+(b-a)*n+(c-a)*m
    pub fn intersect(&self,v: Vector) -> Option<Point> {
        let vals = vec![-v.x,self.b.x-self.a.x,self.c.x-self.a.x,
                        -v.y,self.b.y-self.a.y,self.c.y-self.a.y,
                        -v.z,self.b.z-self.a.z,self.c.z-self.a.z];
        let to_mult = vec![v.p.x-self.a.x,
                            v.p.y-self.a.y,
                            v.p.z-self.a.z];
        let mat = DMatrix3::from_row_vector(3,3,&vals);
        let mat_mult = DMatrix3::from_row_vector(3,1,&to_mult);
        let inv = mat.inverse();
        let result: DMatrix3<f64>;
        match inv {
            Some(inverted) => {
                result = inverted*mat_mult;
            },
            None => {
                return None;
            }
        }
        let t = result[(0,0)];
        if t<=0.+f64::EPSILON {
            return None;
        }
        return Some(Point {x: v.p.x+v.x*t,
                        y: v.p.y+v.y*t,
                        z: v.p.z+v.z*t});
    }
    pub fn normal(&self, p: Point) -> Vector {
        let va = Vector::new(Point::new_blank(),self.b.x-self.a.x,self.b.y-self.a.y,self.b.z-self.a.z);
        let vb = Vector::new(Point::new_blank(),self.c.x-self.a.x,self.c.y-self.a.y,self.c.z-self.a.z);
        let mut result = va.cross_product(&vb);
        result.normalize();
        result.p = p;
        return result;
    }
    pub fn move_along_vector(&mut self, v: Vector, t: f64) {
        self.a.x=self.a.x+v.x*t;
        self.a.y=self.a.y+v.y*t;
        self.a.z=self.a.z+v.z*t;

        self.b.x=self.b.x+v.x*t;
        self.b.y=self.b.y+v.y*t;
        self.b.z=self.b.z+v.z*t;

        self.c.x=self.c.x+v.x*t;
        self.c.y=self.c.y+v.y*t;
        self.c.z=self.c.z+v.z*t;
    }
}

#[test]
fn plane_basic_collision() {
    let plane = Plane{a: Point{x: 1.,y: 1.,z: 0.}, b: Point{x: 1.,y: 0.,z: 0.}, c: Point{x: 0.,y: 1.,z: 0.}};
    let ray = Vector {p: Point{x: 0.,y: 0.,z: 5.},x: 0.,y: 0.,z: -1.};
    let point = plane.intersect(ray);
    match point {
        Some(p) => {
            println!("{:?}", p );
            assert_eq!(p.x,0.);
            assert_eq!(p.y,0.);
            assert_eq!(p.z,0.);
        },
        _ => {
            panic!();
        },
    }
}
