//
// All code here is just a port from Biopython project
// all rights reserved to the original authors.
//

use lazy_static::lazy_static;
use num::pow;
use std::collections::{HashMap, HashSet};

lazy_static! {
    // Normalized flexibility parameters (B-values), average
    // Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).
    static ref FLEXIBILITY_TABLE: HashMap<char, f64> = {
    [
    ("A", 0.984), ("C", 0.906),
    ("E", 1.094), ("D", 1.068),
    ("G", 1.031), ("F", 0.915),
    ("I", 0.927), ("H", 0.950),
    ("K", 1.102), ("M", 0.952),
    ("L", 0.935), ("N", 1.048),
    ("Q", 1.037), ("P", 1.049),
    ("S", 1.046), ("R", 1.008),
    ("T", 0.997), ("W", 0.904),
    ("V", 0.931), ("Y", 0.929),
    ]
        .iter()
        .map(|(s,v)| {(s.to_string().chars().next().expect("failed to get char for amino acid"), (*v) as f64)})
        .collect()
};

    pub static ref MONOISOTOPIC_PROTEIN_WEIGHTS:  HashMap<char, f64> = {
    [
    ("A", 89.047678), ("C", 121.019749),
    ("D", 133.037508), ("E", 147.053158),
    ("F", 165.078979), ("G", 75.032028),
    ("H", 155.069477), ("I", 131.094629),
    ("K", 146.105528), ("L", 131.094629),
    ("M", 149.051049), ("N", 132.053492),
    ("O", 255.158292), ("P", 115.063329),
    ("Q", 146.069142), ("R", 174.111676),
    ("S", 105.042593), ("T", 119.058243),
    ("U", 168.964203), ("V", 117.078979),
    ("W", 204.089878), ("Y", 181.073893),
    ]
        .iter()
        .map(|(s,v)| {(s.to_string().chars().next().expect("failed to get char for amino acid"), (*v) as f64)})
        .collect()
};

    pub static ref PROTEIN_LETTERS: String = String::from("ACDEFGHIKLMNPQRSTVWY");
}

#[derive(Default, Debug)]
pub struct ProteinAnalysis {
    sequence: String,
    amino_acids_content: HashMap<char, usize>,
    amino_acids_percent: HashMap<char, f64>,
    length: usize,
    monoisotopic: bool,
}

impl ProteinAnalysis {
    pub fn new(sequence: &str) -> Self {
        let mut PA = ProteinAnalysis {
            sequence: sequence.to_owned(),
            monoisotopic: true,
            length: sequence.chars().count(),
            ..Default::default()
        };
        PA
    }

    fn count_amino_acids(&self) -> HashMap<char, usize> {
        let mut aa: HashMap<char, usize> = HashMap::with_capacity(PROTEIN_LETTERS.chars().count());
        for c in self.sequence.chars() {
            *aa.entry(c).or_insert(0) += 1;
        }
        for c in PROTEIN_LETTERS.chars() {
            if !aa.contains_key(&c) {
                aa.insert(c, 0);
            }
        }
        aa
    }

    fn get_amino_acids_percent(&self) -> HashMap<char, f64> {
        let mut aa_percent: HashMap<char, f64> =
            HashMap::with_capacity(PROTEIN_LETTERS.chars().count());
        for (aa, v) in &self.amino_acids_content {
            aa_percent.insert(
                *aa,
                *self
                    .amino_acids_content
                    .get(aa)
                    .expect("Failed to get amino acid content") as f64
                    / self.length as f64,
            );
        }
        aa_percent
    }

    fn molecular_weight(&self) -> f64 {
        let mut weigth: f64 = 0.0;
        let water: f64 = 18.010565;
        for c in self.sequence.chars() {
            weigth += MONOISOTOPIC_PROTEIN_WEIGHTS
                .get(&c)
                .expect("failed to retriece monoisotopic protein weight")
        }
        weigth -= (&self.length - 1) as f64 * water;
        weigth
    }

    fn aromaticity(&self) -> f64 {
        let aromatic_aas = "YWF";
        let aromaticity: f64 = 0.0;
        let a: Vec<f64> = self
            .amino_acids_percent
            .iter()
            .map(|(k, v)| *v)
            .collect();
        a.iter().sum()
    }

    // fn flexibility() -> i32 {}
}
