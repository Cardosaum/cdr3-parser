//
// All code here is just a port from Biopython project
// all rights reserved to the original authors.
//

use crate::cdr3::isoelectric_point::IsoelectricPoint;
use lazy_static::lazy_static;
use num::pow;
use std::collections::{HashMap, HashSet};

lazy_static! {
    // Normalized flexibility parameters (B-values), average
    // Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).
    static ref FLEXIBILITY_TABLE: HashMap<char, f64> = {
    [
    ('A', 0.984), ('C', 0.906),
    ('E', 1.094), ('D', 1.068),
    ('G', 1.031), ('F', 0.915),
    ('I', 0.927), ('H', 0.950),
    ('K', 1.102), ('M', 0.952),
    ('L', 0.935), ('N', 1.048),
    ('Q', 1.037), ('P', 1.049),
    ('S', 1.046), ('R', 1.008),
    ('T', 0.997), ('W', 0.904),
    ('V', 0.931), ('Y', 0.929),
    ]
        .iter()
        .map(|(c,f)| {(*c, (*f) as f64)})
        .collect()
};

    pub static ref MONOISOTOPIC_PROTEIN_WEIGHTS:  HashMap<char, f64> = {
    [
    ('A', 89.047678),  ('C', 121.019749),
    ('D', 133.037508), ('E', 147.053158),
    ('F', 165.078979), ('G', 75.032028),
    ('H', 155.069477), ('I', 131.094629),
    ('K', 146.105528), ('L', 131.094629),
    ('M', 149.051049), ('N', 132.053492),
    ('O', 255.158292), ('P', 115.063329),
    ('Q', 146.069142), ('R', 174.111676),
    ('S', 105.042593), ('T', 119.058243),
    ('U', 168.964203), ('V', 117.078979),
    ('W', 204.089878), ('Y', 181.073893),
    ]
        .iter()
        .map(|(c,f)| {(*c, (*f) as f64)})
        .collect()
};

    pub static ref PROTEIN_LETTERS: String = String::from("ACDEFGHIKLMNPQRSTVWY");

    pub static ref DIWV: HashMap<char, HashMap<char, f64>> = {
        // thanks to @kroltan to help with this map.
        let H: HashMap<char, HashMap<char, f64>>  =
        [
        ('A', [('A', 1.0), ('C', 44.94), ('E', 1.0), ('D', -7.49),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', -7.49),
              ('K', 1.0), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 1.0), ('P', 20.26), ('S', 1.0), ('R', 1.0),
              ('T', 1.0), ('W', 1.0), ('V', 1.0), ('Y', 1.0)]),
        ('C', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 20.26),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', 33.60),
              ('K', 1.0), ('M', 33.60), ('L', 20.26), ('N', 1.0),
              ('Q', -6.54), ('P', 20.26), ('S', 1.0), ('R', 1.0),
              ('T', 33.60), ('W', 24.68), ('V', -6.54), ('Y', 1.0)]),
        ('E', [('A', 1.0), ('C', 44.94), ('E', 33.60), ('D', 20.26),
              ('G', 1.0), ('F', 1.0), ('I', 20.26), ('H', -6.54),
              ('K', 1.0), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 20.26), ('P', 20.26), ('S', 20.26), ('R', 1.0),
              ('T', 1.0), ('W', -14.03), ('V', 1.0), ('Y', 1.0)]),
        ('D', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', 1.0), ('F', -6.54), ('I', 1.0), ('H', 1.0),
              ('K', -7.49), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 1.0), ('P', 1.0), ('S', 20.26), ('R', -6.54),
              ('T', -14.03), ('W', 1.0), ('V', 1.0), ('Y', 1.0)]),
        ('G', [('A', -7.49), ('C', 1.0), ('E', -6.54), ('D', 1.0),
              ('G', 13.34), ('F', 1.0), ('I', -7.49), ('H', 1.0),
              ('K', -7.49), ('M', 1.0), ('L', 1.0), ('N', -7.49),
              ('Q', 1.0), ('P', 1.0), ('S', 1.0), ('R', 1.0),
              ('T', -7.49), ('W', 13.34), ('V', 1.0), ('Y', -7.49)]),
        ('F', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 13.34),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', 1.0),
              ('K', -14.03), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 1.0), ('P', 20.26), ('S', 1.0), ('R', 1.0),
              ('T', 1.0), ('W', 1.0), ('V', 1.0), ('Y', 33.601)]),
        ('I', [('A', 1.0), ('C', 1.0), ('E', 44.94), ('D', 1.0),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', 13.34),
              ('K', -7.49), ('M', 1.0), ('L', 20.26), ('N', 1.0),
              ('Q', 1.0), ('P', -1.88), ('S', 1.0), ('R', 1.0),
              ('T', 1.0), ('W', 1.0), ('V', -7.49), ('Y', 1.0)]),
        ('H', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', -9.37), ('F', -9.37), ('I', 44.94), ('H', 1.0),
              ('K', 24.68), ('M', 1.0), ('L', 1.0), ('N', 24.68),
              ('Q', 1.0), ('P', -1.88), ('S', 1.0), ('R', 1.0),
              ('T', -6.54), ('W', -1.88), ('V', 1.0), ('Y', 44.94)]),
        ('K', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', -7.49), ('F', 1.0), ('I', -7.49), ('H', 1.0),
              ('K', 1.0), ('M', 33.60), ('L', -7.49), ('N', 1.0),
              ('Q', 24.64), ('P', -6.54), ('S', 1.0), ('R', 33.60),
              ('T', 1.0), ('W', 1.0), ('V', -7.49), ('Y', 1.0)]),
        ('M', [('A', 13.34), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', 58.28),
              ('K', 1.0), ('M', -1.88), ('L', 1.0), ('N', 1.0),
              ('Q', -6.54), ('P', 44.94), ('S', 44.94), ('R', -6.54),
              ('T', -1.88), ('W', 1.0), ('V', 1.0), ('Y', 24.68)]),
        ('L', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', 1.0),
              ('K', -7.49), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 33.60), ('P', 20.26), ('S', 1.0), ('R', 20.26),
              ('T', 1.0), ('W', 24.68), ('V', 1.0), ('Y', 1.0)]),
        ('N', [('A', 1.0), ('C', -1.88), ('E', 1.0), ('D', 1.0),
              ('G', -14.03), ('F', -14.03), ('I', 44.94), ('H', 1.0),
              ('K', 24.68), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', -6.54), ('P', -1.88), ('S', 1.0), ('R', 1.0),
              ('T', -7.49), ('W', -9.37), ('V', 1.0), ('Y', 1.0)]),
        ('Q', [('A', 1.0), ('C', -6.54), ('E', 20.26), ('D', 20.26),
              ('G', 1.0), ('F', -6.54), ('I', 1.0), ('H', 1.0),
              ('K', 1.0), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 20.26), ('P', 20.26), ('S', 44.94), ('R', 1.0),
              ('T', 1.0), ('W', 1.0), ('V', -6.54), ('Y', -6.54)]),
        ('P', [('A', 20.26), ('C', -6.54), ('E', 18.38), ('D', -6.54),
              ('G', 1.0), ('F', 20.26), ('I', 1.0), ('H', 1.0),
              ('K', 1.0), ('M', -6.54), ('L', 1.0), ('N', 1.0),
              ('Q', 20.26), ('P', 20.26), ('S', 20.26), ('R', -6.54),
              ('T', 1.0), ('W', -1.88), ('V', 20.26), ('Y', 1.0)]),
        ('S', [('A', 1.0), ('C', 33.60), ('E', 20.26), ('D', 1.0),
              ('G', 1.0), ('F', 1.0), ('I', 1.0), ('H', 1.0),
              ('K', 1.0), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 20.26), ('P', 44.94), ('S', 20.26), ('R', 20.26),
              ('T', 1.0), ('W', 1.0), ('V', 1.0), ('Y', 1.0)]),
        ('R', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', -7.49), ('F', 1.0), ('I', 1.0), ('H', 20.26),
              ('K', 1.0), ('M', 1.0), ('L', 1.0), ('N', 13.34),
              ('Q', 20.26), ('P', 20.26), ('S', 44.94), ('R', 58.28),
              ('T', 1.0), ('W', 58.28), ('V', 1.0), ('Y', -6.54)]),
        ('T', [('A', 1.0), ('C', 1.0), ('E', 20.26), ('D', 1.0),
              ('G', -7.49), ('F', 13.34), ('I', 1.0), ('H', 1.0),
              ('K', 1.0), ('M', 1.0), ('L', 1.0), ('N', -14.03),
              ('Q', -6.54), ('P', 1.0), ('S', 1.0), ('R', 1.0),
              ('T', 1.0), ('W', -14.03), ('V', 1.0), ('Y', 1.0)]),
        ('W', [('A', -14.03), ('C', 1.0), ('E', 1.0), ('D', 1.0),
              ('G', -9.37), ('F', 1.0), ('I', 1.0), ('H', 24.68),
              ('K', 1.0), ('M', 24.68), ('L', 13.34), ('N', 13.34),
              ('Q', 1.0), ('P', 1.0), ('S', 1.0), ('R', 1.0),
              ('T', -14.03), ('W', 1.0), ('V', -7.49), ('Y', 1.0)]),
        ('V', [('A', 1.0), ('C', 1.0), ('E', 1.0), ('D', -14.03),
              ('G', -7.49), ('F', 1.0), ('I', 1.0), ('H', 1.0),
              ('K', -1.88), ('M', 1.0), ('L', 1.0), ('N', 1.0),
              ('Q', 1.0), ('P', 20.26), ('S', 1.0), ('R', 1.0),
              ('T', -7.49), ('W', 1.0), ('V', 1.0), ('Y', -6.54)]),
        ('Y', [('A', 24.68), ('C', 1.0), ('E', -6.54), ('D', 24.68),
              ('G', -7.49), ('F', 1.0), ('I', 1.0), ('H', 13.34),
              ('K', 1.0), ('M', 44.94), ('L', 1.0), ('N', 1.0),
              ('Q', 1.0), ('P', 13.34), ('S', 1.0), ('R', -15.91),
              ('T', -7.49), ('W', -9.37), ('V', 1.0), ('Y', 13.34)]),
]
            .iter()
            .map(|(key, data)| (*key, data.iter().map(|&(key, value)| (key, value))
            .collect()))
            .collect();
        H
    };

    pub static ref GRAVY_SCALE: HashMap<char, f64> = {
        [
            ('A', 1.8), ('C', 2.5), ('D', -3.5), ('E', -3.5),
            ('F', 2.8), ('G', -0.4), ('H', -3.2), ('I', 4.5),
            ('K', -3.9), ('L', 3.8), ('M', 1.9), ('N', -3.5),
            ('P', -1.6), ('Q', -3.5), ('R', -4.5), ('S', -0.8),
            ('T', -0.7), ('V', 4.2), ('W', -0.9), ('Y', -1.3),
        ]
            .iter()
            .map(|(c,f)| {(*c, (*f) as f64)})
            .collect()
    };

}

#[derive(Default, Debug)]
pub struct ProteinAnalysis {
    pub sequence: String,
    pub amino_acids_content: HashMap<char, usize>,
    pub amino_acids_percent: HashMap<char, f64>,
    pub length: usize,
    pub monoisotopic: bool,
    pub molecular_weight: f64,
    pub aromaticity: f64,
    pub instability_index: f64,
    pub flexibility: Vec<f64>,
    pub gravy: f64,
    pub isoelectric_point: f64,
    pub charge_at_pH: f64,
    pub secondary_structure_fraction: (f64, f64, f64),
    pub molar_extinction_coefficient: (usize, usize),
}

impl ProteinAnalysis {
    pub fn new(sequence: &str) -> Self {
        let mut PA = ProteinAnalysis {
            sequence: sequence.to_owned(),
            monoisotopic: true,
            length: sequence.chars().count(),
            ..Default::default()
        };
        PA.count_amino_acids();
        PA.get_amino_acids_percent();
        PA.aromaticity();
        PA.molecular_weight();
        PA.instability_index();
        PA.flexibility();
        PA.gravy();
        PA.isoelectric_point();
        PA.charge_at_pH(None);
        PA.secondary_structure_fraction();
        PA.molar_extinction_coefficient();
        PA
    }

    fn count_amino_acids(&mut self) -> HashMap<char, usize> {
        if self.amino_acids_content.is_empty() {
            let mut aa: HashMap<char, usize> =
                HashMap::with_capacity(PROTEIN_LETTERS.chars().count());
            for c in self.sequence.chars() {
                *aa.entry(c).or_insert(0) += 1;
            }
            for c in PROTEIN_LETTERS.chars() {
                if !aa.contains_key(&c) {
                    aa.insert(c, 0);
                }
            }
            self.amino_acids_content = aa;
        }
        return self.amino_acids_content.clone();
    }

    fn get_amino_acids_percent(&mut self) -> HashMap<char, f64> {
        if self.amino_acids_percent.is_empty() {
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
            self.amino_acids_percent = aa_percent;
        }
        return self.amino_acids_percent.clone();
    }

    fn molecular_weight(&mut self) -> f64 {
        if self.molecular_weight == 0.0 && !self.sequence.is_empty() {
            let mut weight: f64 = 0.0;
            let water: f64 = 18.010565;
            for c in self.sequence.chars() {
                weight += MONOISOTOPIC_PROTEIN_WEIGHTS
                    .get(&c)
                    .expect("failed to retriece monoisotopic protein weight")
            }
            weight -= (&self.length - 1) as f64 * water;
            self.molecular_weight = weight;
        }
        return self.molecular_weight;
    }

    fn aromaticity(&mut self) -> f64 {
        let aromatic_aas = "YWF";
        let aromaticity: f64 = 0.0;
        let a: Vec<f64> = self
            .amino_acids_percent
            .iter()
            .filter(|(k, v)| **k == 'Y' || **k == 'W' || **k == 'F')
            .map(|(k, v)| *v)
            .collect();
        self.aromaticity = a.iter().sum();
        return self.aromaticity;
    }

    fn instability_index(&mut self) -> f64 {
        let index = &DIWV;
        let mut score = 0.0_f64;

        for i in 0..(self.length - 1) {
            let this = self.sequence.chars().nth(i).unwrap();
            let next = self.sequence.chars().nth(i + 1).unwrap();
            let dipeptide_value = index[&this][&next];
            score += dipeptide_value;
        }
        self.instability_index = (10.0_f64 / self.length as f64) * score;
        return self.instability_index;
    }

    fn flexibility(&mut self) -> Vec<f64> {
        let window_size: usize = 9;
        if self.flexibility.is_empty() && self.length >= window_size {
            let flexibilities = &FLEXIBILITY_TABLE;
            let weights: Vec<f64> = vec![0.25, 0.4375, 0.625, 0.8125, 1.0];
            let mut scores: Vec<f64> = Vec::with_capacity(self.length - window_size);

            for i in 0..(self.length - window_size) {
                let subsequence: &str = self.sequence.get(i..(i + window_size)).unwrap();
                let mut score: f64 = 0.0;

                for j in 0..(window_size / 2) {
                    let front: char = subsequence.chars().nth(j).unwrap();
                    let back: char = subsequence.chars().nth(window_size - j - 1).unwrap();
                    score += (flexibilities[&front] + flexibilities[&back]) * weights[j];
                }

                let middle: char = subsequence.chars().nth(window_size / 2 + 1).unwrap();
                score += flexibilities[&middle];
                scores.push(score / 5.25);
            }

            self.flexibility = scores;
        }
        return self.flexibility.clone();
    }

    fn gravy(&mut self) -> f64 {
        let scale = &GRAVY_SCALE;
        let g: f64 = self.sequence.chars().into_iter().map(|c| scale[&c]).sum();
        self.gravy = g / self.length as f64;
        return self.gravy;
    }

    fn isoelectric_point(&mut self) -> f64 {
        self.isoelectric_point = IsoelectricPoint::new(&self.sequence, None).isoeletric_point;
        return self.isoelectric_point;
    }

    fn charge_at_pH(&mut self, pH: Option<f64>) -> f64 {
        self.charge_at_pH = IsoelectricPoint::new(&self.sequence, pH).charge_at_pH;
        return self.charge_at_pH;
    }

    fn secondary_structure_fraction(&mut self) -> (f64, f64, f64) {
        let helix: f64 = "VIYFWL".chars().map(|c| self.amino_acids_percent[&c]).sum();
        let turn: f64 = "NPGS".chars().map(|c| self.amino_acids_percent[&c]).sum();
        let sheet: f64 = "EMAL".chars().map(|c| self.amino_acids_percent[&c]).sum();
        self.secondary_structure_fraction = (helix, turn, sheet);
        return self.secondary_structure_fraction;
    }

    fn molar_extinction_coefficient(&mut self) -> (usize, usize) {
        let num_aa = &self.count_amino_acids();
        let mec_reduced: usize = num_aa[&'W'] * 5500 + num_aa[&'Y'] * 1490;
        let mec_cystines: usize = mec_reduced + (num_aa[&'C'] / 2) * 125;
        self.molar_extinction_coefficient = (mec_reduced, mec_cystines);
        return self.molar_extinction_coefficient;
    }
}
