use clap::{crate_version, App, AppSettings};

#[derive(PartialEq)]
pub enum InputFormat {
    CdrOnly,
    Fasta,
}

#[derive(PartialEq)]
pub enum OutputFormat {
    Csv,
    Json,
}

pub fn build_app() -> App<'static, 'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        AppSettings::ColoredHelp
    } else {
        AppSettings::ColorNever
    };

    let app = App::new("cdr3-parser")
        .version(crate_version!())
        .setting(clap_color_setting)
        .about(
            "Output chemical properties and statistics of CDR3 sequences within a `aafreq` file.\n\
             You may choose either `json` or `csv` as output formats.",
        )
        .args_from_usage("<INPUT_FILE> 'Sets the input file to use'")
        .args_from_usage("--cdr-only   'Informs program that the input contains only CDR3VH sequences, separated by a new line'")
        .args_from_usage("-j, --json   'Select `json` as output format'")
        .args_from_usage("-c, --csv    'Select `csv` as output format [default]'");

    app
}
