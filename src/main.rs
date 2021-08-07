// #![feature(map_first_last)]
extern crate bio;
extern crate clap;
extern crate num;
extern crate rayon;
extern crate velcro;

mod app;
mod cdr3;

use app::OutputFormat;
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
                println!("\n`{}` must be a valid file.\n", input_file);
                std::process::exit(1);
            }
        }
        None => {}
    }

    // match matches.value_of("OUTPUT_FILE") {
    //     Some(output_file) => {
    //         if Path::new(output_file).is_file() && matches.is_present("no-clobber") {
    //             println!("File `{}` exists, aborting.", output_file);
    //             std::process::exit(2);
    //         }
    //         if Path::new(output_file).is_dir() {
    //             println!("`{}` is a directory, aborting.", output_file);
    //             std::process::exit(3);
    //         }
    //     }
    //     None => {}
    // }

    let input_file = PathBuf::from(
        matches
            .value_of("INPUT_FILE")
            .get_or_insert("INPUT_FILE")
            .to_string(),
    );

    let output_format = match matches.is_present("json") {
        true => OutputFormat::Json,
        false => OutputFormat::Csv,
    };

    // let output_file = PathBuf::from(
    //     matches
    //         .value_of("OUTPUT_FILE")
    //         .get_or_insert("OUTPUT_FILE")
    //         .to_string(),
    // );

    match pipeline_cdr3(input_file, output_format) {
        Ok(_) => (),
        Err(e) => eprintln!("Error while processing input file: {}", e),
    }
}
