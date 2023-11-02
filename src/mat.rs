//! Define the type of the witness matrices used in the library.
//! 
//! # Example
//! 
//! 
use core::convert::From;
use std::fs::File;
use std::io::Write;

use bincode;
use serde::{Serialize, Deserialize};

use crate::config::{DATA_DIR_PUBLIC, DATA_DIR_PRIVATE};

#[derive(Debug, Clone, Serialize, Deserialize)]
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

trait ToFile {
    fn to_file(&self, file_name: &str, public: bool) -> std::io::Result<()>;
    fn from_file(file_name: &str, public:bool) -> std::io::Result<Self> 
    where
        Self: Sized;
}  

impl<T: Serialize + for<'de> Deserialize<'de> > ToFile for Mat<T> {
    fn to_file(&self, file_name: &str, public: bool) -> std::io::Result<()> {
        
        let dir = if public {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let mut file = File::create(format!("{}{}", dir, file_name))?;
        
        let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
        file.write_all(&encoded)?;
        Ok(())
    }

    fn from_file(file_name: &str, public: bool) -> std::io::Result<Self> {
        let dir = if public {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let file = File::open(format!("{}{}", dir, file_name))?;
        let decoded: Self = bincode::deserialize_from(file).unwrap();
        Ok(decoded)
    }
}

impl<T: Serialize + for<'de> Deserialize<'de> > ToFile for Vec<T> {
    fn to_file(&self, file_name: &str, public: bool) -> std::io::Result<()> {
        let dir = if public {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let mut file = File::create(format!("{}{}", dir, file_name))?;

        let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
        file.write_all(&encoded)?;
        Ok(())
    }

    fn from_file(file_name: &str, public: bool) -> std::io::Result<Self> {
        let dir = if public {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let file = File::open(format!("{}{}", dir, file_name))?;
        let decoded: Self = bincode::deserialize_from(file).unwrap();
        Ok(decoded)
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
    use std::vec;

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

        mat_1.to_file("mat_test.dat", true).unwrap();
        let mat_1_read: Mat<GtElement> = Mat::from_file("mat_test.dat", true)
        .unwrap();
        assert_eq!(mat_1, mat_1_read);

        let vec_1: Vec<GtElement> = vec![GtElement::from(1), GtElement::from(2)];
        vec_1.to_file("vec_test.dat", false).unwrap();
        let vec_1_read: Vec<GtElement> = Vec::from_file( "vec_test.dat", false)
        .unwrap();
        assert_eq!(vec_1, vec_1_read);
    }
}