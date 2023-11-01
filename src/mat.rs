//! Define the type of the witness matrices used in the library.
//! 
//! # Example
//! 
//! 
use core::convert::From;

#[derive(Debug, Clone)]
pub struct Mat<T> {
    pub id: String,
    pub shape: (usize, usize),
    pub data: Vec<(usize, usize, T)>,
}

impl< T: From<u64> > Mat<T> {
    pub fn new(
        id: &str, 
        shape: (usize, usize)) -> Self {
        let id = String::from(id);
        let data: Vec<(usize, usize, T)> = Vec::new();
        Self {id, shape, data}
    }

    pub fn new_from_data_vec(
        id: &str, 
        shape: (usize, usize),
        data: Vec<(usize, usize, T)>) -> Self {
        let id = String::from(id);
        Self {id, shape, data}
    }

    pub fn new_from_u64_vec(
        id: &str, 
        shape: (usize, usize),
        data: Vec<(usize, usize, u64)>) -> Self {
        let id = String::from(id);
        let data = data.iter().map(
            |(row, col, val)| (*row, *col, T::from(*val))
        ).collect();
        
        Self {id, shape, data}
    }

    pub fn push(&mut self, row: usize, col: usize, val: T) {
        self.data.push((row, col, val));
    }
}

impl<T: PartialEq + Clone> PartialEq for Mat<T> {

    fn eq(&self, other: &Self) -> bool {
        if self.shape != other.shape {
            return false;
        }

        let mut self_data = self.data.clone();
        let mut other_data = other.data.clone();

        self_data.sort_by(|a, b| a.0.cmp(&b.0));
        self_data.sort_by(|a, b| a.1.cmp(&b.1));

        other_data.sort_by(|a, b| a.0.cmp(&b.0));
        other_data.sort_by(|a, b| a.1.cmp(&b.1));

        self_data == other_data
    }
}

impl<T: PartialEq + Clone> Eq for Mat<T> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::GtElement;

    #[test]
    fn test_mat() {
        let mut mat_1:Mat<GtElement>  = Mat::new("test_mat", (2, 2));
        assert_eq!(mat_1.id, "test_mat");
        assert_eq!(mat_1.shape, (2, 2));

        mat_1.push(1, 1, GtElement::from(3));
        mat_1.push(1, 2, GtElement::from(1));

        assert_eq!(mat_1.data[1], (1, 2, GtElement::from(1)));

        let mut mat_2:Mat<GtElement>  = Mat::new("test_mat_2", (2, 2));
        mat_2.push(1, 2, GtElement::from(1));
        mat_2.push(1, 1, GtElement::from(3));
        
        assert_eq!(mat_1, mat_2);
    }
}