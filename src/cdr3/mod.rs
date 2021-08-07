pub mod cdr3;

// pub mod isoelectric_point;

pub mod protein_analysis;

pub mod prelude {
    pub use crate::cdr3::cdr3::{
        aa_groups, aromaticity, create_cdr3_dict, extract_cdr3, molecular_weight, pipeline_cdr3,
        write_cdr3_attributes, CDR3Prop,
    };
    // pub use crate::cdr3::isoelectric_point::IsoelectricPoint;
    pub use crate::cdr3::protein_analysis::{ProteinAnalysis, PROTEIN_WEIGHTS};
}
