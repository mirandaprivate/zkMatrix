//! Utility functions for file IO operations
// 
// 
use std::fs::File;
use std::io::{Read, Write};

use bincode;
use serde::{Serialize, Deserialize};

use crate::config::{DATA_DIR_PUBLIC, DATA_DIR_PRIVATE};

pub trait FileIO {
    
    fn to_file(
        &self, file_name: String, public_data: bool
    ) -> std::io::Result<()>
    where
        Self: Serialize,
    {
        let dir = 
            if public_data {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let mut file = File::create(
            format!("{}{}", dir, file_name)
        )?;
        
        let encoded: Vec<u8> = bincode::serialize(&self).unwrap();
        file.write_all(&encoded)?;

        // println!(" -----------------------------------------------------  ");
        println!(
            " ===== File {}{} written; Size {} Bytes ", 
            dir, 
            file_name,
            encoded.len()
        );
        println!(" -----------------------------------------------------  ");

        Ok(())
    }

    fn from_file(
        file_name: String, public_data:bool
    ) -> std::io::Result<Self> 
    where
        Self: Sized + for<'de> Deserialize<'de>,
    {
        let dir = 
            if public_data {DATA_DIR_PUBLIC} else {DATA_DIR_PRIVATE};
        let mut file = File::open(
            format!("{}{}", dir, file_name)
        ).map_err(
            |e| 
            std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!(
                    "Error: {:?} when opening file {}", e, file_name
                )
            )
        )?;
        
        let mut encoded: Vec<u8> = Vec::new();
        file.read_to_end(&mut encoded).map_err(
            |e| 
            std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!(
                    "Error: {:?} when reading u8 in file {}", e, file_name
                )
            )
        )?;

        let decoded: Self = 
            bincode::deserialize(&encoded)
            .map_err(
                |e| 
                std::io::Error::new(
                    std::io::ErrorKind::Other, 
                    format!(
                        "Error: {:?} when deserilizing file {}", e, file_name
                    )
                )
            )?;
        Ok(decoded)
    }
}  


impl<T: Serialize + for<'de> Deserialize<'de> > FileIO for Vec<T> {}