use crate::field::*;


pub struct Matrix <CoeffT> 
    where CoeffT: PID
{
    dim: (usize, usize),
    val: Vec<Vec<CoeffT>>,
    rank: Option<usize>,
}

impl<CoeffT> Matrix 
    where CoeffT: Field
{
    pub fn from(val: &[CoeffT]) -> Matrix<CoeffT>{
        Matrix {
            val: vec![val],
            dim: (val.0, val.1),
            rank: self.rank,
        }
    }
}

#[cfg(test)]
mod matrix_test {
    #[test]
    fn init_test() {
        let array = [[1, 2], [3, 4], [5, 6]];
        array: &[[_]] = array.iter().map(|&x| x.iter().map(|elem| elem.rational())).collect();
        let m = Matrix::from();
    }
}
