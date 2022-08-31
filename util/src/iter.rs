pub struct Counter {
    count: usize,
    max: usize,
}

impl Counter {
    pub fn new() -> Counter {
        Counter { 
            count: 0,
            max: std::usize::MAX,
        }
    }

    pub fn new_with_max(max_in: usize) -> Counter {
        Counter { 
            count: 0,
            max: max_in,
        }
    }

    pub fn new_with_min(min_in: usize) -> Counter {
        Counter { 
            count: min_in,
            max: std::usize::MAX,
        }
    }

    pub fn new_with_min_and_max(min_in: usize, max_in: usize) -> Counter {
        Counter { 
            count: min_in,
            max: max_in,
        }
    }
}

impl Iterator for Counter {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.count < self.max {
            self.count += 1;
            Some(self.count - 1)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod counter_test {
    #[test]
    fn initialization(){
        use crate::iter::Counter;

        let v: Vec<_> = Counter::new_with_max(5).collect();
        assert_eq!(v, vec![0,1,2,3,4]);

        let v: Vec<_> = Counter::new_with_min_and_max(11, 18).collect();
        assert_eq!(v, vec![11,12,13,14,15,16,17]);
    }
}