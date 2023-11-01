//! Define the type of the witness matrices used in the library.
//! 
//! # Example
//! 
#[derive(Debug, Clone)]
struct Mat<'a, T> {
    pub id: &'a str,
    pub shape: (usize, usize),
    pub data: Vec<(usize, usize, T)>,
}

impl<'a, T> Mat<'a, T> {
    pub fn new(
        id: &'a str, 
        shape: (usize, usize)) -> Self {
        let data: Vec<(usize, usize, T)> = Vec::new();
        Self {
            id,
            shape,
            data,
        }
    }

    pub fn push(&mut self, row: usize, col: usize, val: T) {
        self.data.push((row, col, val));
    }
}

impl<T: PartialEq + Clone> PartialEq for Mat<'_, T> {

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

impl<T: PartialEq + Clone> Eq for Mat<'_, T> {}

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