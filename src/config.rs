pub const NUM_THREADS: usize = 8;

pub const LOG_DIM: usize = 12;

pub const Q: usize = 2usize.pow(LOG_DIM as u32+1);

pub const SQRT_MATRIX_DIM: usize = 2usize.pow(LOG_DIM as u32/2);
pub const MATRIX_DIM: usize = SQRT_MATRIX_DIM * SQRT_MATRIX_DIM;

pub const LOG_DIM_TEST: usize = 6;
pub const Q_TEST: usize = 2usize.pow(LOG_DIM_TEST as u32+1);

pub const SQRT_MATRIX_DIM_TEST: usize = 2usize.pow(3);
pub const MATRIX_DIM_TEST: usize = SQRT_MATRIX_DIM_TEST * SQRT_MATRIX_DIM_TEST;

pub const DATA_DIR_PUBLIC: &str = "data/public/";
pub const DATA_DIR_PRIVATE: &str = "data/private/";