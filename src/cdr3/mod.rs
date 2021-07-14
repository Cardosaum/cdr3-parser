pub mod cdr3;

pub mod proteins;

pub mod prelude {
    pub use crate::cdr3::cdr3::{
        aa_groups, aromaticity, build_cdr3_struct, create_cdr3_dict, extract_cdr3,
        get_cdr3_sequences_attributes, molecular_weight, write_cdr3_attributes, CDR3Prop,
    };
    pub use crate::cdr3::proteins::IsoelectricPoint;
}
