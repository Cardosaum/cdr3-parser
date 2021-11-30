extern crate bio;
extern crate clap;
extern crate num;
extern crate rayon;
extern crate velcro;

mod app;
mod cdr3;

use app::{InputFormat, OutputFormat};
use cdr3::prelude::*;

use std::env;
use std::path::{Path, PathBuf};

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

    let input_format = match matches.is_present("cdr-only") {
        true => InputFormat::CdrOnly,
        false => InputFormat::Fasta,
    };

    match pipeline_cdr3(input_file, input_format, output_format) {
        Ok(_) => (),
        Err(e) => eprintln!("Error while processing input file: {}", e),
    }
}
