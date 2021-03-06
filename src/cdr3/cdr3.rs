use lazy_static::lazy_static;
use rayon::prelude::*;
use regex::Regex;
use serde::Serialize;
// use serde_json::Result;
use serde_with::serde_as;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::iter::FromIterator;
use std::path::PathBuf;
use std::str;
use velcro::hash_map;

use crate::app::{InputFormat, OutputFormat};
use crate::cdr3::prelude::*;

#[serde_as]
#[derive(Debug, Serialize)]
pub struct CDR3Prop {
    pub cdr3: String,
    pub quantity: usize,
    pub length: usize,
    pub MW: f64,
    pub AV: f64,
    IP: f64,
    // flex: f64,
    // gravy: f64,
    // SSF_Helix: f64,
    // SSF_Turn: f64,
    // SSF_Sheet: f64,
    #[serde_as(as = "Vec<(_, _)>")]
    pub aa_quantity: HashMap<char, u64>, // n_A: usize,
                                         // n_C: usize,
                                         // n_D: usize,
                                         // n_E: usize,
                                         // n_F: usize,
                                         // n_G: usize,
                                         // n_H: usize,
                                         // n_I: usize,
                                         // n_K: usize,
                                         // n_L: usize,
                                         // n_M: usize,
                                         // n_N: usize,
                                         // n_P: usize,
                                         // n_Q: usize,
                                         // n_R: usize,
                                         // n_S: usize,
                                         // n_T: usize,
                                         // n_V: usize,
                                         // n_W: usize,
                                         // n_Y: usize,
                                         // aliphatic: usize,
                                         // aromatic: usize,
                                         // neutral: usize,
                                         // positive: usize,
                                         // negative: usize,
                                         // invalid: usize
}

pub fn molecular_weight(sequence: &str) -> f64 {
    let monoisotopic_protein_weights: HashMap<&str, f64> = hash_map! {
        "A": 89.047678,
        "C": 121.019749,
        "D": 133.037508,
        "E": 147.053158,
        "F": 165.078979,
        "G": 75.032028,
        "H": 155.069477,
        "I": 131.094629,
        "K": 146.105528,
        "L": 131.094629,
        "M": 149.051049,
        "N": 132.053492,
        "O": 255.158292,
        "P": 115.063329,
        "Q": 146.069142,
        "R": 174.111676,
        "S": 105.042593,
        "T": 119.058243,
        "U": 168.964203,
        "V": 117.078979,
        "W": 204.089878,
        "Y": 181.073893,
    };

    let mut mw: f64 = 0.0;

    sequence.bytes().for_each(|c| {
        match monoisotopic_protein_weights
            .get_key_value(&str::from_utf8(&[c]).expect("molecular wieght @ cdr3"))
        {
            Some(aa) => {
                mw += aa.1;
            }
            None => {}
        };
    });

    return mw;
}

pub fn aromaticity(sequence: &str) -> f64 {
    return ((sequence.matches("Y").count() as f64)
        + (sequence.matches("W").count() as f64)
        + (sequence.matches("F").count() as f64))
        / (sequence.chars().count() as f64);
}

pub fn create_cdr3_dict(cdr3_sequence: Vec<String>) -> BTreeMap<String, u64> {
    let mut cdr3_dict = BTreeMap::new();

    for cdr3 in cdr3_sequence {
        let cdr3_key = cdr3_dict.entry(cdr3).or_insert(0);
        *cdr3_key += 1;
    }

    return cdr3_dict;
}

pub fn aa_groups<'seq>(sequence: String) -> BTreeMap<String, HashSet<&'seq str>> {
    let mut group = BTreeMap::new();
    group.insert(
        "aliphatic".to_string(),
        HashSet::from_iter(["G", "A", "P", "V", "L", "I", "M"]),
    );
    group.insert("aromatic".to_string(), HashSet::from_iter(["F", "Y", "W"]));
    group.insert(
        "neutral".to_string(),
        HashSet::from_iter(["S", "T", "C", "N", "Q"]),
    );
    group.insert("positive".to_string(), HashSet::from_iter(["K", "H", "R"]));
    group.insert("negative".to_string(), HashSet::from_iter(["D", "E"]));

    let mut result = BTreeMap::new();
    result.insert("aliphatic".to_string(), 0);
    result.insert("aromatic".to_string(), 0);
    result.insert("neutral".to_string(), 0);
    result.insert("positive".to_string(), 0);
    result.insert("negative".to_string(), 0);
    result.insert("invalid".to_string(), 0);

    // let mut not_listed = HashSet::new();

    // for aa in std::str::FromStr(sequence) {
    //     println!("{}", aa);
    // }

    return group;
}

// pub fn build_cdr3_struct(cdr3: String) -> CDR3Prop {
//     CDR3Prop {
//         cdr3: cdr3.clone(),
//         quantity: 0,
//         length: cdr3.chars().count(),
//         MW: molecular_weight(&cdr3),
//         AV: aromaticity(&cdr3),
//         IP: IsoelectricPoint::new(&cdr3, None).isoeletric_point,
//         // flex: f64,
//         // gravy: f64,
//         // SSF_Helix: f64,
//         // SSF_Turn: f64,
//         // SSF_Sheet: f64,
//         aa_quantity: cdr3.chars().fold(HashMap::<char, u64>::new(), |mut k, c| {
//             *k.entry(c).or_insert(0) += 1;
//             k
//         }),
//         // n_A: cdr3.matches("A").count(),
//         // n_C: cdr3.matches("C").count(),
//         // n_D: cdr3.matches("D").count(),
//         // n_E: cdr3.matches("E").count(),
//         // n_F: cdr3.matches("F").count(),
//         // n_G: cdr3.matches("G").count(),
//         // n_H: cdr3.matches("H").count(),
//         // n_I: cdr3.matches("I").count(),
//         // n_K: cdr3.matches("K").count(),
//         // n_L: cdr3.matches("L").count(),
//         // n_M: cdr3.matches("M").count(),
//         // n_N: cdr3.matches("N").count(),
//         // n_P: cdr3.matches("P").count(),
//         // n_Q: cdr3.matches("Q").count(),
//         // n_R: cdr3.matches("R").count(),
//         // n_S: cdr3.matches("S").count(),
//         // n_T: cdr3.matches("T").count(),
//         // n_V: cdr3.matches("V").count(),
//         // n_W: cdr3.matches("W").count(),
//         // n_Y: cdr3.matches("Y").count(),
//         // aliphatic: u64,
//         // aromatic: u64,
//         // neutral: u64,
//         // positive: u64,
//         // negative: u64,
//         // invalid: u64,
//     }
// }

// pub fn get_cdr3_sequences_attributes(cdr3_dict: BTreeMap<String, u64>) -> Vec<CDR3Prop> {
//     let mut tmp = cdr3_dict
//         .par_iter()
//         .map(|cdr3| build_cdr3_struct(cdr3.0.to_string()))
//         .collect::<Vec<CDR3Prop>>();
//     tmp.par_sort_unstable_by_key(|p| p.cdr3.clone());
//     return tmp;
// }

pub fn extract_cdr3(
    mut sequences: Vec<String>,
    input_format: &InputFormat,
) -> HashMap<String, usize> {
    match input_format {
        InputFormat::Fasta => {
            lazy_static! {
                static ref CDR3_REGEX: Regex =
                    Regex::new(r"((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)").expect("extract cdr3 1");
            }

            sequences.par_iter_mut().for_each(|s| {
                let cdr3 = CDR3_REGEX.captures(&s);
                if cdr3.is_some() {
                    *s = String::from(
                        cdr3.expect("extract cdr3 2")
                            .get(7)
                            .map_or("", |m| m.as_str()),
                    );
                }
            });

            let mut distinct_sequences: HashMap<String, usize> = HashMap::new();

            for i in sequences {
                // distinct_sequences.insert(i);
                *distinct_sequences.entry(i).or_insert(0) += 1;
            }

            return distinct_sequences;
        }
        InputFormat::CdrOnly => {
            let mut distinct_sequences: HashMap<String, usize> = HashMap::new();

            for i in sequences {
                *distinct_sequences.entry(i).or_insert(0) += 1;
            }

            return distinct_sequences;
        }
    }
}

// pub fn write_cdr3_attributes(cdr3_attributes: Vec<CDR3Prop>, output_file: PathBuf) -> Result<()> {
pub fn write_cdr3_header() -> Result<(), Box<dyn std::error::Error>> {
    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.write_record(&[
        "sequence",
        "quantity",
        "length",
        "molecular_weight",
        "aromaticity",
        "charge_at_pH",
        "gravy",
        "instability_index",
        "isoelectric_point",
        "ssf_helix",
        "ssf_turn",
        "ssf_sheet",
        "aliphatic",
        "aromatic",
        "neutral",
        "positive",
        "negative",
    ])
    .expect("failed to write record header");
    wtr.flush().expect("failed to flush header");
    Ok(())
}

pub fn write_cdr3_attributes(
    sequence: &str,
    quantity: usize,
    output_format: &OutputFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut P = ProteinAnalysis::new(sequence, quantity, None);
    if output_format == &OutputFormat::Json {
        let _ = writeln!(std::io::stdout(), "{}", serde_json::json!(P));
    } else {
        let mut wtr = csv::Writer::from_writer(io::stdout());
        if let Err(err) = wtr.write_record(&[
            P.sequence.clone(),
            quantity.to_string(),
            P.length.to_string(),
            P.molecular_weight().to_string(),
            P.aromaticity().to_string(),
            P.charge_at_pH(None).to_string(),
            P.gravy().to_string(),
            P.instability_index().to_string(),
            P.isoelectric_point(None, None, None).to_string(),
            P.secondary_structure_fraction().0.to_string(),
            P.secondary_structure_fraction().1.to_string(),
            P.secondary_structure_fraction().2.to_string(),
            P.aa_aliphatic.to_string(),
            P.aa_aromatic.to_string(),
            P.aa_neutral.to_string(),
            P.aa_positive.to_string(),
            P.aa_negative.to_string(),
        ]) {
            let _ = writeln!(std::io::stderr(), "Error in write record: {}", err);
        }
        if let Err(err) = wtr.flush() {
            if err.kind() == std::io::ErrorKind::BrokenPipe {
                std::process::exit(0);
            }
            let _ = writeln!(std::io::stderr(), "Error while flushing: {}", err);
        }
    }
    Ok(())
}

fn parse_file(input_file: PathBuf, input_format: &InputFormat) -> Vec<String> {
    let input_file_handler = File::open(input_file).expect("Failed to open <INPUT_FILE>");
    let reader = BufReader::new(input_file_handler);
    let mut seqs = Vec::<String>::new();

    match input_format {
        InputFormat::Fasta => {
            for line in reader.lines() {
                let line = match line {
                    Ok(line) => line,
                    Err(error) => {
                        println!("Error reading file: {:?}", error);
                        std::process::exit(4);
                    }
                };

                if !line.starts_with(">") && !line.starts_with("#") && !line.starts_with("*") {
                    seqs.push(line.to_string());
                }
            }
        }
        InputFormat::CdrOnly => {
            for line in reader.lines() {
                let line = match line {
                    Ok(line) => line,
                    Err(error) => {
                        println!("Error reading file: {:?}", error);
                        std::process::exit(4);
                    }
                };

                seqs.push(line.to_string());
            }
        }
    }

    return seqs;
}

fn is_valid_cdr3(sequence: &str) -> bool {
    for c in sequence.chars() {
        match PROTEIN_WEIGHTS.get(&c) {
            Some(_) => (),
            None => return false,
        };
    }
    return true;
}

pub fn pipeline_cdr3(
    input_file: PathBuf,
    input_format: InputFormat,
    output_format: OutputFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    let sequences = extract_cdr3(parse_file(input_file, &input_format), &input_format);

    if output_format == OutputFormat::Csv {
        match write_cdr3_header() {
            Ok(_) => (),
            Err(e) => eprintln!("Error writing header: {:?}", e),
        };
    }

    sequences.par_iter().for_each(|(s, n)| {
        if is_valid_cdr3(s) {
            match write_cdr3_attributes(s, *n, &output_format) {
                Ok(_) => (),
                Err(e) => eprintln!("Error on pipeline: {:?}", e),
            };
            // TODO: create a log flag to save all invalid sequences into a file (?)
            // println!("{:?}", s);
        }
    });

    Ok(())
}
