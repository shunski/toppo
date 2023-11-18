use alg::lin_alg::{Matrix, Vector, SparseMatrix};
use alg::commutative::{PID, Zero, One};

pub trait Function: PartialEq + Clone {
    type Domain;
    type Codomain;
    fn eval(&self, x: Self::Domain) -> Self::Codomain ;
}

trait SmoothFn: Function
    where Self::Codomain: PID
{
    fn jacobian(&self) -> Matrix<Self::Codomain>;
}

impl std::ops::Add for AtomicSmoothFn {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        [self, rhs].into_iter().sum()
    }
}

impl std::ops::AddAssign for AtomicSmoothFn {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl std::iter::Sum for AtomicSmoothFn {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut out: Vector<_, _> = iter
            .filter(|s| s != &Self::Zero)
            .map(|s| {
                    if let Self::Sum(v) = s {
                        v
                    } else {
                        Vector::from(s)
                    }
                } )
            .sum();

        // if 'out' has only one element and if that element has coefficient 1.0, then unpack it
        if out.len() == 1 {
            let (coeff, basis) = out.remove(0);
            basis * coeff
        } 
        // if 'out' is empty, return 'Self::Zero'
        else if out.len() == 0 {
            Self::Zero
        } else {
            Self::Sum(out)
        }
    }
}

impl std::ops::Sub for AtomicSmoothFn {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let out = if let Self::Sum(v) = self {
            if let Self::Sum(w) = rhs {
                v - w
            } else {
                v - rhs
            }
        } else if let Self::Sum(w) = rhs {
            -w+self
        } else {
            Vector::from(self)-rhs
        };
        Self::Sum(out)
    }
}

impl std::ops::SubAssign for AtomicSmoothFn {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl std::ops::Neg for AtomicSmoothFn {
    type Output = Self;
    fn neg(self) -> Self {
        let v = if let Self::Sum(x) = self {
            -x
        } else {
            -Vector::from(self)
        };
        Self::Sum(v)
    }
}


impl std::ops::Mul for AtomicSmoothFn {
    type Output = Self;
    fn mul(mut self, mut rhs: Self) -> Self::Output {
        // If either factor is 'Self::Zero', return 'Self::Zero'
        if self == Self::Zero || rhs==Self::Zero {return Self::Zero}
        
        // if either factor is 'Self::One', return the other factor.
        if self == Self::One {
            return rhs;
        } else if rhs == Self::One {
            return self;
        }

        // if either factor is 'Self::Sum' with one element, decompose the Sum structure 
        if let Self::Sum(v) = &mut self {
            if v.len() == 1 {
                let (coeff, basis) = v.remove(0);
                return (basis * rhs) * coeff
            }
        } else if let Self::Sum(v) = &mut rhs {
            if v.len() == 1 {
                let (coeff, basis) = v.remove(0);
                return (self * basis) * coeff
            }
        }

        // The rest of the cases: 
        [self, rhs].into_iter().product()
    }
}

impl std::ops::MulAssign for AtomicSmoothFn {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

impl std::iter::Product for AtomicSmoothFn {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut out: Vector<_, _> = iter
            .filter(|s| s != &Self::One)
            .map(|s| {
                    if let Self::Product(v) = s {
                        v
                    } else {
                        Vector::from(s)
                    }
                } )
            .sum();

        if out.iter().any(|(_,f)| f == &Self::Zero ) {
            return Self::Zero
        }

        // if 'out' has only one element and if that element has coefficient 1.0, then unpack it
        if out.len() == 1 {
            let (coeff, basis) = out.remove(0);
            basis.powi(coeff)
        } 
        // if 'out' is empty, return 'Self::One'
        else if out.len() == 0 {
            Self::One
        } else {
            Self::Product(out)
        }
    }
}

impl std::ops::Div for AtomicSmoothFn {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.powi(-1)
    }
}

impl std::ops::DivAssign for AtomicSmoothFn {
    fn div_assign(&mut self, rhs: Self) {
        *self = self.clone() / rhs;
    }
}

// scalar multiplication
impl std::ops::Mul<f64> for AtomicSmoothFn {
    type Output = Self;
    fn mul(self, other: f64) -> Self {
        if self == Self::Zero {
            self
        } else if other == 1.0 {
            self
        } else if let Self::Sum( mut v ) = self {
            v.iter_mut().for_each(|(a,_)| *a *= other );
            Self::Sum(v)
        } else {
            Self::Sum(Vector::from(self) * other)
        }
    }
}


// scalar division
impl std::ops::Div<f64> for AtomicSmoothFn {
    type Output = Self;
    fn div(self, other: f64) -> Self {
        if other == 0.0 {
            panic!("cannot devide a function by zero!");
        } else if self == Self::Zero {
            self
        } else if other == 1.0 {
            self
        } else if let Self::Sum( mut v ) = self {
            v.iter_mut().for_each(|(a,_)| *a /= other );
            Self::Sum(v)
        } else {
            Self::Sum(Vector::from(self) / other)
        }
    }
}

impl Zero for AtomicSmoothFn {
    fn zero() -> Self {
        Self::Zero
    }
}


impl One for AtomicSmoothFn {
    fn one() -> Self  {
        Self::One
    }
}


#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum AtomicSmoothFn {
    Zero,
    One,

    // Projection map to the n-th coordinate
    Projection(usize),
    

    RealPower(Box<AtomicSmoothFn>, f64),
    Sin(Box<AtomicSmoothFn>),
    Cos(Box<AtomicSmoothFn>),

    // summands in 'Sum' must be sorted by 'PartialOrd' of 'AtomicSmoothFn' at any time.
    // 'Sum' cannot directly contain another 'Sum' in it.
    // 'Product' cannot directly contain 'Zero' in it.
    Sum(Vector<f64, AtomicSmoothFn>),

    // factors in 'Product' must be sorted by 'PartialOrd' of 'AtomicSmoothFn' at any time.
    // 'Product' cannot directly contain another 'Product' in it.
    // 'Product' cannot directly contain 'Zero' or 'One' in it.
    Product(Vector<i32, AtomicSmoothFn>),
}

impl Function for AtomicSmoothFn {
    type Domain = Matrix<f64>;
    type Codomain = f64;

    fn eval(&self, x: Self::Domain) -> Self::Codomain {
        match self {
            Self::Zero => 0.0,
            Self::One => 1.0,
            Self::Projection(c) => x[(*c, 0)],
            Self::RealPower(f, a) => (f.eval(x)).powf(*a),
            Self::Sin(f) => (f.eval(x)).sin(),
            Self::Cos(f) => (f.eval(x)).cos(),
            Self::Sum(v) => v.iter().map(|(a,f)| a * f.eval(x.clone()) ).sum(),
            Self::Product(v) => v.iter().map(|(a, f)| f.eval(x.clone()).powi(*a) ).product()
        }
    }
}

impl alg::commutative::Symbol for AtomicSmoothFn {
    fn symbol()->String {
        "C^{infinity}".to_string()
    }
}

impl PID for AtomicSmoothFn {
    fn divides(self, _: Self) -> bool {
        true
    }

    fn euclid_div(self, d: Self) -> (Self, Self) {
        (self / d, Self::Zero)
    }

    fn is_prime(self) -> bool {
        false
    }

    fn n_prime_factors(self) -> usize {
        1
    }

    fn simplify(self) -> (Self, Self) {
        (self, Self::One)
    }
}

impl AtomicSmoothFn {
    pub fn partial_derivative(self, n: usize) -> AtomicSmoothFn {
        match self {
            Self::Zero => Self::Zero,
            Self::One => Self::Zero,
            Self::Projection(m) => if n==m {Self::One} else {Self::Zero} 
            Self::Sin(f) => Self::Cos(f.clone()) * f.partial_derivative(n),
            Self::Cos(f) => -Self::Sin(f.clone()) * f.partial_derivative(n),
            Self::RealPower(f, a) => Self::RealPower(f.clone(), a-1.0) * f.partial_derivative(n) * (a as f64),
            Self::Sum(summands) => {
                let out = summands.into_iter()
                    .map(|(a, f)| 
                            (a, f.partial_derivative(n))
                        )
                    .filter(|(_, f)| f != &Self::Zero )
                    .map(|(a,f)| if let Self::Sum(g) = f {
                            g * a
                        } else {
                            Vector::from(f) * a
                        })
                    .sum();

                // before returning, add 'Self::Zero' to reduce some redundant structure
                Self::Sum(out) + Self::Zero
            }
            Self::Product(factors) => {
                (0..factors.len()).map(|i| {
                    let mut clone = factors.clone();
                    let (e, f) = clone.remove(i);
                    let dfdn = f.clone().partial_derivative(n);
                    if dfdn == Self::Zero {
                        Self::Zero
                    } else if e == 1 && dfdn == Self::One && clone.len() == 1 {
                        let (coeff, basis) = clone.remove(0);
                        basis.powi(coeff)
                    } else {
                        Self::Product(clone) * f.powi(e-1) * e as f64 * dfdn
                    }
                })
                .filter(|f| f != &Self::Zero)
                .sum()
            },
        }
    }

    pub fn jacobian(self, dim_domain: usize) -> Matrix<AtomicSmoothFn> {
        let mut out = Matrix::zero(1, dim_domain);

        for i in 0..out.size().1 {
            out[(0,i)] = self.clone().partial_derivative(i);
        }
        out
    }

    pub fn powi(self, exponent: i32) -> Self {
        if exponent == 1 {
            return self;
        } else if exponent == 0 {
            return Self::One;
        }

        if let Self::Product(v) = self {
            let v = v * exponent;
            Self::Product(v)
        } else {
            Self::Product(Vector::from(self) * exponent)
        }
    }

    pub fn powf(self, exponent: f64) -> Self {
        if exponent != 1.0 {
            Self::RealPower(Box::new(self), exponent)
        } else {
            self
        }
    }

    pub fn sin(self) -> Self {
        Self::Sin(Box::new(self))
    }

    pub fn cos(self) -> Self {
        Self::Cos(Box::new(self))
    }
}


impl std::fmt::Display for AtomicSmoothFn {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Zero => write!(f, "0")?,
            Self::One => write!(f, "1")?,
            Self::Projection(m) => write!(f, "x[{m}]")?,
            Self::Sin(g) => write!(f, "sin({g})")?,
            Self::Cos(g) => write!(f, "cos({g})")?,
            Self::RealPower(g, a) => write!(f, "({g})^{a}")?,
            Self::Sum(summands) => write!(f, "{summands}")?,
            Self::Product(factors) => {
                for (e, factor) in factors.iter() {
                    write!(f, " {factor}^{e}")?
                }
            }
        };
        write!(f, "")
    }
}

pub struct AtomicSmoothFnBuilder {
    handle: AtomicSmoothFn,
    dim_domain: usize,
}

impl AtomicSmoothFnBuilder {
    pub fn new (dim_domain: usize) -> Self {
        Self { handle: AtomicSmoothFn::Zero, dim_domain }
    }

    pub fn set(mut self, setter: impl Fn(&mut AtomicSmoothFn, Vec<Coordinate>)) -> AtomicSmoothFn {
        setter(&mut self.handle, (0..self.dim_domain).map(|c| Coordinate(c)).collect::<Vec<_>>());
        self.handle
    }
}


#[derive(Clone, PartialEq)]
struct Var {}

impl Function for Var {
    type Domain = f64;
    type Codomain = f64;
    fn eval(&self, val:Self::Domain ) -> f64 {
        val
    }
}

#[derive(Clone, PartialEq)]
pub struct MultiVarFn {
    fns: Vec<AtomicSmoothFn>,
    dim_domain: usize,
}

impl MultiVarFn {
    pub fn new((dim_codomain, dim_domain): (usize, usize)) -> MultiVarFnBuilder {
        MultiVarFnBuilder{ 
            handle: Self{ fns: vec![AtomicSmoothFn::Zero; dim_codomain], dim_domain }
        }
    }

    pub fn jacobian(self) -> Matrix<AtomicSmoothFn> {
        let n = self.fns.len();
        let mut out = Matrix::zero(n, self.dim_domain);
        for (i, f) in self.fns.into_iter().enumerate() {
            out[(i,..)].write(&f.jacobian(self.dim_domain));
        }
        out
    }
}

impl Function for MultiVarFn {
    type Domain = Matrix<f64>;
    type Codomain = Matrix<f64>;

    fn eval(&self, x: Self::Domain) -> Self::Codomain {
        let mut out = Matrix::zero(self.fns.len(), 1);
        for (i, f) in self.fns.iter().enumerate() {
            out[(i,0)] = f.eval(x.clone());
        }
        out
    }
}

pub struct MultiVarFnBuilder {
    handle: MultiVarFn,
}

impl MultiVarFnBuilder {
    pub fn set(mut self, setter: impl Fn(usize, &mut AtomicSmoothFn, Vec<Coordinate>)) -> MultiVarFn {
        for (i, f) in self.handle.fns.iter_mut().enumerate() {
            setter(i, f, (0..self.handle.dim_domain).map(|c| Coordinate(c)).collect::<Vec<_>>())
        }
        self.handle
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Coordinate(usize);

pub mod coordinate_impl {
    use crate::{AtomicSmoothFn, Coordinate};

    macro_rules! atomic_fn_interface_impl {
        ($(($fn: ident)), *) => {
            impl Coordinate {
                $(
                    pub fn $fn(self) -> AtomicSmoothFn {
                        AtomicSmoothFn::Projection(self.0).$fn()
                    }
                )*
            }
        };

        ($(($fn: ident, $input: ty)), *) => {
            impl Coordinate {
                $(
                    pub fn $fn(self, input: $input) -> AtomicSmoothFn {
                        AtomicSmoothFn::Projection(self.0).$fn(input)
                    }
                )*
            }
        };
    }

    atomic_fn_interface_impl!{ (sin),  (cos) }
    atomic_fn_interface_impl!{ 
        (powi, i32), 
        (powf, f64)
    }

    macro_rules! ops_impl {
        ($(($op: ident, $op_fn: ident, $symbol: tt)), *) => {$(
            // operation between two 'Coordinate's
            impl std::ops::$op for Coordinate {
                type Output = AtomicSmoothFn;
                fn $op_fn(self, rhs: Self) -> Self::Output {
                    AtomicSmoothFn::Projection(self.0) $symbol AtomicSmoothFn::Projection(rhs.0)
                }
            }

            // operation between 'Coordinate' and 'AtomicSmoothFn'
            impl std::ops::$op<AtomicSmoothFn> for Coordinate {
                type Output = AtomicSmoothFn;
                fn $op_fn(self, rhs: AtomicSmoothFn) -> Self::Output {
                    AtomicSmoothFn::Projection(self.0) $symbol rhs
                }
            }

            // operation between 'AtomicSmoothFn' and 'Coordinate'
            impl std::ops::$op<Coordinate> for AtomicSmoothFn {
                type Output = AtomicSmoothFn;
                fn $op_fn(self, rhs: Coordinate) -> Self::Output {
                    self $symbol AtomicSmoothFn::Projection(rhs.0)
                }
            }
        )*};
    }

    ops_impl!{(Add, add, +), (Sub, sub, -), (Mul, mul, *), (Div, div, /)}

    impl std::ops::Mul<f64> for Coordinate {
        type Output = AtomicSmoothFn;
        fn mul(self, rhs: f64) -> Self::Output {
            AtomicSmoothFn::Projection(self.0) * rhs
        }
    }

    impl std::ops::Div<f64> for Coordinate {
        type Output = AtomicSmoothFn;
        fn div(self, rhs: f64) -> Self::Output {
            AtomicSmoothFn::Projection(self.0) / rhs
        }
    }
}

#[cfg(test)]
mod test {
    use alg::lin_alg::Matrix;

    use crate::{
        MultiVarFn,
        AtomicSmoothFn, AtomicSmoothFnBuilder
    };

    #[test]
    fn partial_derivative() {
        let f = AtomicSmoothFnBuilder::new(2).set(|f, x| 
            *f = - x[0].cos() * x[0] * x[1]
        );


        let dfdx0 = f.clone().partial_derivative(0);
        let ans = AtomicSmoothFnBuilder::new(2).set(|f, x| 
            *f = x[0].sin() * x[0] * x[1] - x[0].cos() * x[1]
        );
        assert!(
            dfdx0 == ans,
            "dfdx0={:?} \nnot equals \nans={:?}", dfdx0, ans,
        );

        let dfdx1 = f.partial_derivative(1);
        let ans = AtomicSmoothFnBuilder::new(2).set(|f, x| 
            *f = - x[0].cos() * x[0]
        );
        assert!(
            dfdx1 == ans,
            "dfdx0={:?} \nnot equals \nans={:?}", dfdx1, ans,
        );
    }

    #[test]
    fn functions_of_multi_variable() {
        let f = MultiVarFn::new((1, 3))
            .set( |_, f, x| {
                *f = -x[0].sin() * x[1].powi(3) + x[0]*x[2];
            });
        
        let j = f.jacobian();
        let j_ans = {
            let mut out = Matrix::<AtomicSmoothFn>::zero(1, 3);
            out[(0,0)] = MultiVarFn::new((1, 3))
                .set( |_, f, x| {
                    *f = -x[0].cos() * x[1].powi(3) + x[2];
                }).fns[0].clone();
            out[(0,1)] = MultiVarFn::new((1, 3))
                .set( |_, f, x| {
                    *f = - x[0].sin() * x[1].powi(2) * 3.0;
                }).fns[0].clone();
            out[(0,2)] = MultiVarFn::new((1, 3))
                .set( |_, f, x| {
                    *f = x[0] * 1.0;
                }).fns[0].clone();
            out
        };

        assert!(j==j_ans, "j={} \nnot equals \nj_ans={}", &j, &j_ans);
    }
}


impl Function for Matrix<AtomicSmoothFn> {
    type Codomain = Matrix<f64>;
    type Domain = Matrix<f64>;
    fn eval(&self, other: Self::Domain) -> Self::Codomain {
        let mut out = Matrix::zero(self.size().0, self.size().1);
        for i in 0..self.size().0 {
            for j in 0..self.size().1 {
                out[(i,j)] = self[(i,j)].clone().eval(other.clone());
            }
        }
        out
    }
}

impl Function for SparseMatrix<AtomicSmoothFn> {
    type Codomain = SparseMatrix<f64>;
    type Domain = Matrix<f64>;
    fn eval(&self, other: Self::Domain) -> Self::Codomain {
        let mut out = Vec::new();
        for i in 0..self.size().0 {
            for (j, f) in self.row_iter(i) {
                out.push(((i, j), f.eval(other.clone())));
            }
        }

        let initializer = |(i,j)| {
            let result = out.binary_search_by(|(idx, _)| idx.cmp(&(i,j)) );
            match result {
                Ok(x) => Some(out[x].1),
                Err(_) => None
            }
        };

        SparseMatrix::from(((self.size().0, self.size().1), initializer))
    }
}