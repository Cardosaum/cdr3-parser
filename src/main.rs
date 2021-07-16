// #![feature(map_first_last)]
extern crate bio;
extern crate clap;
extern crate num;
extern crate rayon;
extern crate velcro;

mod app;
mod cdr3;

use cdr3::prelude::*;

use clap::{App, Arg, Error, SubCommand};
use csv::Writer;
use serde_json::json;
use std::fs::{self, File, OpenOptions};
use std::io::{self, prelude::*, BufReader};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::{env, str};

fn main() {
    let matches = app::build_app().get_matches_from(env::args_os());

    match matches.value_of("INPUT_FILE") {
        Some(input_file) => {
            if !Path::new(input_file).is_file() {
                println!("`{}` must be a valid file.", input_file);
                std::process::exit(1);
            }
        }
        None => {}
    }

    match matches.value_of("OUTPUT_FILE") {
        Some(output_file) => {
            if Path::new(output_file).is_file() && matches.is_present("no-clobber") {
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
        matches
            .value_of("INPUT_FILE")
            .get_or_insert("INPUT_FILE")
            .to_string(),
    );
    let output_file = PathBuf::from(
        matches
            .value_of("OUTPUT_FILE")
            .get_or_insert("OUTPUT_FILE")
            .to_string(),
    );

    // // let cdr3_dict = extract_cdr3(input_file, output_file);
    // let sequences = extract_cdr3(parse_file(input_file));

    // // for x in &sequences {
    // //     println!("{}", x);
    // // }

    // let mut cdr3_dict = create_cdr3_dict(sequences);
    // let cdr3_prop = build_cdr3_struct(cdr3_dict.first_entry().unwrap().key().to_string());
    // // println!("{}", molecular_weight(&cdr3_prop.cdr3));
    // // println!("{}", aromaticity(&cdr3_prop.cdr3));
    // // println!("{:?}", cdr3_prop);

    // // println!("{:?}", cdr3_prop.cdr3);
    // // println!("{:?}", &cdr3_dict.first_entry().unwrap().key());
    // // println!("{:?}", &cdr3_dict);

    // let mut cdr3_sequences_attributes = get_cdr3_sequences_attributes(cdr3_dict);

    // // for x in cdr3_sequences_attributes {
    // //     println!("{:}", x.cdr3);
    // // }

    // // cdr3_dict.into_par_iter().map(|cdr3| {build_cdr3_struct(cdr3.0.to_string())}).collect_into_vec(cdr3_sequences_attributes);

    // // for cdr3_sequence in cdr3_dict {
    // //     // println!("{:?}", cdr3_sequence);
    // //     let m = build_cdr3_struct(cdr3_sequence.0);
    // //     println!("{:?}", m);
    // // }

    // // aa_groups("SAGTKJLAS".to_string());

    // write_cdr3_attributes(cdr3_sequences_attributes, output_file);
    // //
    // // println!("{:?}", IsoelectricPoint::new("GATTACA", None));
    // // println!("{:?}", IsoelectricPoint::new("FIVESK", None));

    println!("{:?}", ProteinAnalysis::new("FIVESK").molecular_weight);
    println!("{:?}", ProteinAnalysis::new("GATTACA").molecular_weight);
    println!(
        "{:?}",
        ProteinAnalysis::new("FIVESKVIESLTY").molecular_weight
    );

    // let P = ProteinAnalysis::new("FIVESKVIESLTY");
    // println!("\n{}\n", json!(&P).to_string());

    // let mut wtr = csv::Writer::from_writer(io::stdout());
    // wtr.write_record(&["sequence", "length", "molecular_weight"]);
    // wtr.write_record(&[P.sequence, P.length.to_string(), P.molecular_weight.to_string()]);
    // wtr.flush();
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
