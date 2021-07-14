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
    static ref PROTEIN_LETTERS: String = String::from("ACDEFGHIKLMNPQRSTVWY");
}

#[derive(Debug)]
pub struct IsoelectricPoint {
    sequence: String,
    aa_content: HashMap<char, usize>,
    pos_pKs: HashMap<String, f64>,
    neg_pKs: HashMap<String, f64>,
    charged_ass_content: HashMap<String, f64>,
    pub isoeletric_point: f64,
}

impl IsoelectricPoint {
    pub fn new(sequence: &str, pH: Option<f64>) -> Self {
        let s = sequence.to_ascii_uppercase();
        let a = IsoelectricPoint::count_amino_acids(&s);
        let p = IsoelectricPoint::pK_table(&s);
        let c = IsoelectricPoint::select_charged(&a);
        let mut tmp = IsoelectricPoint {
            sequence: s,
            aa_content: a,
            pos_pKs: p.0,
            neg_pKs: p.1,
            charged_ass_content: c,
            isoeletric_point: 0.0,
        };
        let i = IsoelectricPoint::pi(&tmp, pH, None, None);
        tmp.isoeletric_point = i;
        tmp
    }

    fn count_amino_acids(sequence: &String) -> HashMap<char, usize> {
        let mut aa: HashMap<char, usize> = HashMap::with_capacity(20); // There is only 20 valid amino acids
        for c in sequence.chars() {
            *aa.entry(c).or_insert(0) += 1;
        }
        for c in PROTEIN_LETTERS.chars() {
            if !aa.contains_key(&c) {
                aa.insert(c, 0);
            }
        }
        aa
    }

    fn pK_table(sequence: &String) -> (HashMap<String, f64>, HashMap<String, f64>) {
        let mut pos_pKs: HashMap<String, f64> = POSITIVE_PKS.clone();
        let mut neg_pKs: HashMap<String, f64> = NEGATIVE_PKS.clone();
        let nterm = sequence
            .chars()
            .last()
            .expect("Sequence was empty")
            .to_string();
        let cterm = sequence
            .chars()
            .rev()
            .last()
            .expect("Sequence was empty")
            .to_string();
        let mut pKnterminal = PKNTERMINAL.clone();
        let mut pKcterminal = PKCTERMINAL.clone();

        if let Some(n) = pKnterminal.get_mut(&nterm) {
            *n = pos_pKs.get("Nterm").expect("Failed to get 'Nterm'").clone();
        }
        if let Some(n) = pKcterminal.get_mut(&cterm) {
            *n = neg_pKs.get("Cterm").expect("Failed to get 'Cterm'").clone();
        }

        return (pos_pKs, neg_pKs);
    }

    fn select_charged(aa_content: &HashMap<char, usize>) -> HashMap<String, f64> {
        let mut charged: HashMap<String, f64> = HashMap::new();
        for aa in CHARGED_AAS.iter() {
            if let Some(v) = aa_content.get(&aa) {
                charged.insert(aa.to_string(), (*v) as f64);
            }
        }
        charged.insert("Nterm".to_owned(), 1.0);
        charged.insert("Cterm".to_owned(), 1.0);
        charged
    }

    fn charge_at_pH(&self, pH: f64) -> f64 {
        let mut positive_charge: f64 = 0.0;
        for (aa, pK) in &self.pos_pKs {
            let partial_charge: f64 = 1.0_f64 / (10.0_f64.powf(pH - pK) + 1.0_f64);
            positive_charge += self
                .charged_ass_content
                .get(aa)
                .expect("Failed to retrive amino acid charge")
                * partial_charge;
        }

        let mut negative_charge: f64 = 0.0;
        for (aa, pK) in &self.neg_pKs {
            let partial_charge: f64 = 1.0_f64 / (10.0_f64.powf(pK - pH) + 1.0_f64);
            negative_charge += self
                .charged_ass_content
                .get(aa)
                .expect("Failed to retrive amino acid charge")
                * partial_charge;
        }
        return positive_charge - negative_charge;
    }

    fn pi(&self, mut pH: Option<f64>, mut min_: Option<f64>, mut max_: Option<f64>) -> f64 {
        let mut charge = IsoelectricPoint::charge_at_pH(&self, pH.unwrap_or(7.775));
        while max_.unwrap_or(12.0) - min_.unwrap_or(4.05) > 0.0001_f64 {
            if charge > 0.0 {
                min_ = Some(pH.unwrap_or(7.775));
            } else {
                max_ = Some(pH.unwrap_or(7.775));
            }
            pH = Some((max_.unwrap_or(12.0) + min_.unwrap_or(4.05)) / 2_f64);
            charge = IsoelectricPoint::charge_at_pH(&self, pH.unwrap_or(7.775));
        }
        return pH.expect("Failed to compute Isoeletric Point");
    }
}
