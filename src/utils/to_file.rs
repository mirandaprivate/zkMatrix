use std::fs::File;
use std::io::Write;

use bincode;
use serde::{Serialize, Deserialize};

use crate::config::{DATA_DIR_PUBLIC, DATA_DIR_PRIVATE};

pub trait FileIO {
    
    fn to_file(&self, file_name: String, public_data: bool) -> std::io::Result<()>
    where
        Self: Serialize,
    {
        let dir = if public_data {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let mut file = File::create(format!("{}{}.dat", dir, file_name))?;
        
        let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
        file.write_all(&encoded)?;
        Ok(())
    }

    fn from_file(file_name: String, public_data:bool) -> std::io::Result<Self> 
    where
        Self: Sized + for<'de> Deserialize<'de>,
    {
        let dir = if public_data {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let file = File::open(format!("{}{}.dat", dir, file_name))?;
        let decoded: Self = bincode::deserialize_from(file)
        .map_err(
            |e| 
            std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!("Error: {:?} when reading {}", e, file_name))
            )?;
        Ok(decoded)
    }
}  