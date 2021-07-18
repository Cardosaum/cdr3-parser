use clap::{crate_version, App, AppSettings, Arg};

pub fn build_app() -> App<'static, 'static> {
    let clap_color_setting = if std::env::var_os("NO_COLOR").is_none() {
        AppSettings::ColoredHelp
    } else {
        AppSettings::ColorNever
    };

    let mut app = App::new("cdr3-parser")
        .version(crate_version!())
        .setting(clap_color_setting)
        .about("Convert a `aafreq` file into a tidy csv with chemical properties of the sequences")
        .args_from_usage("<INPUT_FILE>           'Sets the input file to use'")
        // .args_from_usage("<OUTPUT_FILE>          'Sets the output file to write to'")
        .args_from_usage("-n, --no-clobber       'Do not overwrite an existing file'");

    app
}
