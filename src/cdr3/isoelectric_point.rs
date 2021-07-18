//
// port of IsoelectricPoint.py from Biopython project.
// all copyrights reserved to the original authors.
//
//
// # Copyright 2003 Yair Benita.  All rights reserved.
// # Revisions copyright 2020 by Tianyi Shi.  All rights reserved.
// # This file is part of the Biopython distribution and governed by your
// # choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
// # Please see the LICENSE file that should have been included as part of this
// # package.
// """Calculate isoelectric points of polypeptides using methods of Bjellqvist.
//
// pK values and the methos are taken from::
//
//     * Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F.,
//     Sanchez, J.-Ch., Frutiger, S. & Hochstrasser, D.F.
//     The focusing positions of polypeptides in immobilized pH gradients can be
//     predicted from their amino acid sequences. Electrophoresis 1993, 14,
//     1023-1031.
//
//     * Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
//     Reference points for comparisons of two-dimensional maps of proteins from
//     different human cell types defined in a pH scale where isoelectric points
//     correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
//
// I designed the algorithm according to a note by David L. Tabb, available at:
// http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf
// """
//

use lazy_static::lazy_static;
use num::pow;
use std::collections::{HashMap, HashSet};

use crate::cdr3::protein_analysis::PROTEIN_LETTERS;

lazy_static! {
    static ref POSITIVE_PKS: HashMap<String, f64> = {
        [("Nterm", 7.5), ("K", 10.0), ("R", 12.0), ("H", 5.98)]
            .iter()
            .map(|(s, v)| (s.to_string(), *v))
            .collect()
    };
    static ref NEGATIVE_PKS: HashMap<String, f64> = {
        [
            ("Cterm", 3.55),
            ("D", 4.05),
            ("E", 4.45),
            ("C", 9.0),
            ("Y", 10.0),
        ]
        .iter()
        .map(|(s, v)| (s.to_string(), *v))
        .collect()
    };
    static ref PKCTERMINAL: HashMap<String, f64> = {
        [("D", 4.55), ("E", 4.75)]
            .iter()
            .map(|(s, v)| (s.to_string(), *v))
            .collect()
    };
    static ref PKNTERMINAL: HashMap<String, f64> = {
        [
            ("A", 7.59),
            ("M", 7.0),
            ("S", 6.93),
            ("P", 8.36),
            ("T", 6.82),
            ("V", 7.44),
            ("E", 7.7),
        ]
        .iter()
        .map(|(s, v)| (s.to_string(), *v))
        .collect()
    };
    static ref CHARGED_AAS: HashSet<char> = {
        ["K", "R", "H", "D", "E", "C", "Y"]
            .iter()
            .map(|s| {
                s.to_string()
                    .chars()
                    .next()
                    .expect("Failed to generate Hashset")
            })
            .collect()
    };
    static ref DEFAULT_PH: f64 = 7.775;
}

#[derive(Default, Debug)]
pub struct IsoelectricPoint {
    sequence: String,
    aa_content: HashMap<char, usize>,
    pos_pKs: HashMap<String, f64>,
    neg_pKs: HashMap<String, f64>,
    charged_aas_content: HashMap<String, f64>,
    pub isoeletric_point: Option<f64>,
    pub charge_at_pH: Option<f64>,
}

impl IsoelectricPoint {
    pub fn new(sequence: &str, pH: Option<f64>) -> Self {
        let mut tmp = IsoelectricPoint {
            sequence: sequence.to_ascii_uppercase(),
            ..Default::default()
        };
        tmp.count_amino_acids();
        tmp.pK_table();
        tmp.select_charged();
        tmp.charge_at_pH(pH);
        tmp.pi(pH, None, None);
        tmp
    }

    fn count_amino_acids(&mut self) -> HashMap<char, usize> {
        if self.aa_content.is_empty() {
            let mut aa: HashMap<char, usize> = HashMap::with_capacity(20); // There is only 20 valid amino acids
            for c in self.sequence.chars() {
                *aa.entry(c).or_insert(0) += 1;
            }
            for c in PROTEIN_LETTERS.chars() {
                if !aa.contains_key(&c) {
                    aa.insert(c, 0);
                }
            }
            self.aa_content = aa;
        }
        return self.aa_content.clone();
    }

    fn pK_table(&mut self) -> (HashMap<String, f64>, HashMap<String, f64>) {
        if self.pos_pKs.is_empty() || self.neg_pKs.is_empty() {
            let mut pos_pKs: HashMap<String, f64> = POSITIVE_PKS.clone();
            let mut neg_pKs: HashMap<String, f64> = NEGATIVE_PKS.clone();
            let nterm = self
                .sequence
                .chars()
                .last()
                .expect("Sequence was empty")
                .to_string();
            let cterm = self
                .sequence
                .chars()
                .rev()
                .last()
                .expect("Sequence was empty")
                .to_string();
            let mut pKnterminal = PKNTERMINAL.clone();
            let mut pKcterminal = PKCTERMINAL.clone();

            if let Some(n) = pKnterminal.get_mut(&nterm) {}
            if let Some(n) = pKcterminal.get_mut(&cterm) {
                *n = neg_pKs.get("Cterm").expect("Failed to get 'Cterm'").clone();
            }
            self.pos_pKs = pos_pKs;
            self.neg_pKs = neg_pKs;
        }
        return (self.pos_pKs.clone(), self.neg_pKs.clone());
    }

    fn select_charged(&mut self) -> HashMap<String, f64> {
        if self.charged_aas_content.is_empty() {
            let mut charged: HashMap<String, f64> = HashMap::new();
            for aa in CHARGED_AAS.iter() {
                if let Some(v) = self.aa_content.get(&aa) {
                    charged.insert(aa.to_string(), (*v) as f64);
                }
            }
            charged.insert("Nterm".to_owned(), 1.0);
            charged.insert("Cterm".to_owned(), 1.0);
            self.charged_aas_content = charged;
        }
        return self.charged_aas_content.clone();
    }

    pub fn charge_at_pH(&mut self, mut pH: Option<f64>) -> f64 {
        if self.charge_at_pH.is_none() {
            pH = match pH {
                Some(p) => Some(p),
                None => Some(*DEFAULT_PH),
            };
            let mut positive_charge: f64 = 0.0;
            for (aa, pK) in &self.pos_pKs {
                let partial_charge: f64 =
                    1.0_f64 / (10.0_f64.powf(pH.unwrap_or(*DEFAULT_PH) - pK) + 1.0_f64);
                positive_charge += self
                    .charged_aas_content
                    .get(aa)
                    .expect("Failed to retrive amino acid charge")
                    * partial_charge;
            }

            let mut negative_charge: f64 = 0.0;
            for (aa, pK) in &self.neg_pKs {
                let partial_charge: f64 =
                    1.0_f64 / (10.0_f64.powf(pK - pH.unwrap_or(*DEFAULT_PH)) + 1.0_f64);
                negative_charge += self
                    .charged_aas_content
                    .get(aa)
                    .expect("Failed to retrive amino acid charge")
                    * partial_charge;
            }

            self.charge_at_pH = Some(positive_charge - negative_charge);
        }

        return self.charge_at_pH.unwrap();
    }

    pub fn pi(&mut self, mut pH: Option<f64>, mut min_: Option<f64>, mut max_: Option<f64>) -> f64 {
        if self.isoeletric_point.is_none() {
            let mut charge = self.charge_at_pH(pH);
            while max_.unwrap_or(12.0) - min_.unwrap_or(4.05) > 0.0001_f64 {
                if charge > 0.0 {
                    min_ = Some(pH.unwrap_or(*DEFAULT_PH));
                } else {
                    max_ = Some(pH.unwrap_or(*DEFAULT_PH));
                }
                pH = Some((max_.unwrap_or(12.0) + min_.unwrap_or(4.05)) / 2_f64);
                charge = self.charge_at_pH(pH);
            }
            self.isoeletric_point = Some(pH.expect("Failed to compute Isoeletric Point"));
        }
        return self.isoeletric_point.unwrap();
    }
}
