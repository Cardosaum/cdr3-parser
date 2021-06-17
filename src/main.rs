#![feature(map_first_last)]
extern crate bio;
extern crate clap;
extern crate rayon;
extern crate velcro;

// Import some modules
// use bio::alphabets;
// use bio::data_structures::bwt::{bwt, less, Occ};
// use bio::data_structures::fmindex::{FMIndex, FMIndexable};
// use bio::data_structures::suffix_array::suffix_array;
// use bio::io::fastq;
// use bio::io::fastq::FastqRead;

use clap::{App, Arg, Error, SubCommand};
use csv::Writer;
use lazy_static::lazy_static;
use rayon::prelude::*;
use regex::Regex;
use serde::{Deserialize, Serialize};
use serde_json::Result;
use serde_with::serde_as;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{self, prelude::*, BufReader};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::str;
use velcro::hash_map;

#[serde_as]
#[derive(Debug, Serialize)]
struct CDR3Prop {
    cdr3: String,
    quantity: usize,
    length: usize,
    MW: f64,
    AV: f64,
    // IP: f64,
    // flex: f64,
    // gravy: f64,
    // SSF_Helix: f64,
    // SSF_Turn: f64,
    // SSF_Sheet: f64,
    #[serde_as(as = "Vec<(_, _)>")]
    aa_quantity: HashMap<char, u64>, // n_A: usize,
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

fn main() {
    let cli_args = App::new("cdrparser")
        .version("0.1")
        .author("Matheus C. <cardosaum@pm.me>")
        .about("Convert a `aafreq` file into a tidy csv with chemical properties of the sequences")
        .args_from_usage("<INPUT_FILE>           'Sets the input file to use'")
        .args_from_usage("<OUTPUT_FILE>          'Sets the output file to write to'")
        .args_from_usage("-n, --no-clobber       'Do not overwrite an existing file'")
        .get_matches();

    match cli_args.value_of("INPUT_FILE") {
        Some(input_file) => {
            if !Path::new(input_file).is_file() {
                println!("`{}` must be a valid file.", input_file);
                std::process::exit(1);
            }
        }
        None => {}
    }

    match cli_args.value_of("OUTPUT_FILE") {
        Some(output_file) => {
            if Path::new(output_file).is_file() && cli_args.is_present("no-clobber") {
                println!("File `{}` exists, aborting.", output_file);
                std::process::exit(2);
            }
            if Path::new(output_file).is_dir() {
                println!("`{}` is a directory, aborting.", output_file);
                std::process::exit(3);
            }
        }
        None => {}
    }

    let input_file = PathBuf::from(
        cli_args
            .value_of("INPUT_FILE")
            .get_or_insert("INPUT_FILE")
            .to_string(),
    );
    let output_file = PathBuf::from(
        cli_args
            .value_of("OUTPUT_FILE")
            .get_or_insert("OUTPUT_FILE")
            .to_string(),
    );

    // let cdr3_dict = extract_cdr3(input_file, output_file);
    let sequences = extract_cdr3(parse_file(input_file));

    // for x in &sequences {
    //     println!("{}", x);
    // }

    let mut cdr3_dict = create_cdr3_dict(sequences);
    let cdr3_prop = build_cdr3_struct(cdr3_dict.first_entry().unwrap().key().to_string());
    // println!("{}", molecular_weight(&cdr3_prop.cdr3));
    // println!("{}", aromaticity(&cdr3_prop.cdr3));
    // println!("{:?}", cdr3_prop);

    // println!("{:?}", cdr3_prop.cdr3);
    // println!("{:?}", &cdr3_dict.first_entry().unwrap().key());
    // println!("{:?}", &cdr3_dict);

    let mut cdr3_sequences_attributes: Vec<CDR3Prop> = cdr3_dict
        .par_iter()
        .map(|cdr3| build_cdr3_struct(cdr3.0.to_string()))
        .collect();

    cdr3_sequences_attributes.par_sort_unstable_by_key(|p| p.cdr3.clone());

    // for x in cdr3_sequences_attributes {
    //     println!("{:}", x.cdr3);
    // }

    // cdr3_dict.into_par_iter().map(|cdr3| {build_cdr3_struct(cdr3.0.to_string())}).collect_into_vec(cdr3_sequences_attributes);

    // for cdr3_sequence in cdr3_dict {
    //     // println!("{:?}", cdr3_sequence);
    //     let m = build_cdr3_struct(cdr3_sequence.0);
    //     println!("{:?}", m);
    // }

    // aa_groups("SAGTKJLAS".to_string());

    write_cdr3_attributes(cdr3_sequences_attributes, output_file);
}

fn parse_file(mut input_file: PathBuf) -> Vec<String> {
    let input_file_handler = File::open(input_file).expect("Failed to open <INPUT_FILE>");
    let reader = BufReader::new(input_file_handler);
    let mut seqs = Vec::<String>::new();

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

    return seqs;
}

fn extract_cdr3(mut sequences: Vec<String>) -> Vec<String> {
    lazy_static! {
        static ref CDR3_REGEX: Regex =
            Regex::new(r"((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)").unwrap();
    }

    sequences.par_iter_mut().for_each(|s| {
        let cdr3 = CDR3_REGEX.captures(&s).unwrap();
        *s = String::from(cdr3.get(7).map_or("", |m| m.as_str()));
    });

    return sequences;
}

fn create_cdr3_dict(cdr3_sequence: Vec<String>) -> BTreeMap<String, u64> {
    let mut cdr3_dict = BTreeMap::new();

    for cdr3 in cdr3_sequence {
        let cdr3_key = cdr3_dict.entry(cdr3).or_insert(0);
        *cdr3_key += 1;
    }

    return cdr3_dict;
}

fn aa_groups<'seq>(sequence: String) -> BTreeMap<String, HashSet<&'seq str>> {
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

fn build_cdr3_struct(cdr3: String) -> CDR3Prop {
    CDR3Prop {
        cdr3: cdr3.clone(),
        quantity: 0,
        length: cdr3.chars().count(),
        MW: molecular_weight(&cdr3),
        AV: aromaticity(&cdr3),
        // IP: f64,
        // flex: f64,
        // gravy: f64,
        // SSF_Helix: f64,
        // SSF_Turn: f64,
        // SSF_Sheet: f64,
        aa_quantity: cdr3.chars().fold(HashMap::<char, u64>::new(), |mut k, c| {
            *k.entry(c).or_insert(0) += 1;
            k
        }),
        // n_A: cdr3.matches("A").count(),
        // n_C: cdr3.matches("C").count(),
        // n_D: cdr3.matches("D").count(),
        // n_E: cdr3.matches("E").count(),
        // n_F: cdr3.matches("F").count(),
        // n_G: cdr3.matches("G").count(),
        // n_H: cdr3.matches("H").count(),
        // n_I: cdr3.matches("I").count(),
        // n_K: cdr3.matches("K").count(),
        // n_L: cdr3.matches("L").count(),
        // n_M: cdr3.matches("M").count(),
        // n_N: cdr3.matches("N").count(),
        // n_P: cdr3.matches("P").count(),
        // n_Q: cdr3.matches("Q").count(),
        // n_R: cdr3.matches("R").count(),
        // n_S: cdr3.matches("S").count(),
        // n_T: cdr3.matches("T").count(),
        // n_V: cdr3.matches("V").count(),
        // n_W: cdr3.matches("W").count(),
        // n_Y: cdr3.matches("Y").count(),
        // aliphatic: u64,
        // aromatic: u64,
        // neutral: u64,
        // positive: u64,
        // negative: u64,
        // invalid: u64,
    }
}

fn molecular_weight(sequence: &str) -> f64 {
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
        // println!("{}", c);
        let aa = match monoisotopic_protein_weights.get_key_value(&str::from_utf8(&[c]).unwrap()) {
            Some(aa) => {
                // println!("{:?}", aa);
                mw += aa.1;
            }
            None => {}
        };
    });

    return mw;
}

fn aromaticity(sequence: &str) -> f64 {
    return ((sequence.matches("Y").count() as f64)
        + (sequence.matches("W").count() as f64)
        + (sequence.matches("F").count() as f64))
        / (sequence.chars().count() as f64);
}

fn write_cdr3_attributes(cdr3_attributes: Vec<CDR3Prop>, output_file: PathBuf) -> Result<()> {
    // println!("{:?}", output_file);

    let j = serde_json::to_string(&cdr3_attributes)?;
    println!("{}", j);

    // let mut wtr = Writer::from_path(output_file)?;

    // wtr.serialize(cdr3_attributes)?;
    // println!("{}", output_file);
    // for x in cdr3_attributes {
    //     println!("{:?}", x);
    //     wtr.serialize(x)?;
    // }
    // wtr.flush()?;
    // let m: Vec<String> = cdr3_attributes
    //     .iter()
    //     .map(|c| c.cdr3.clone())
    //     .collect();

    // println!("{:?}", m);

    return Ok(());
}
