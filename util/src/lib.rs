#![allow(unused)]
pub mod iter;

pub mod bitwise {
    pub fn count_ones(mut data: u64) ->  usize {
        let mut count = 0;
        while data != 0 {
            count+=1;
            data &= data-1;
        }
        return count;
    }
}