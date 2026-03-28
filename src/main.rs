use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::process;

use cpc2_rs::error::Cpc2Error;
use cpc2_rs::{analyze_reader, write_results, AnalyzeOptions};

fn main() {
    match parse_args(env::args().skip(1)) {
        Ok(Command::Help) => {
            print_help();
        }
        Ok(Command::Run(config)) => {
            if let Err(err) = run(config) {
                eprintln!("[ERROR] {err}");
                process::exit(1);
            }
        }
        Err(err) => {
            eprintln!("[ERROR] {err}");
            eprintln!();
            print_help();
            process::exit(1);
        }
    }
}

#[derive(Debug)]
struct CliConfig {
    input: PathBuf,
    output: PathBuf,
    reverse: bool,
    include_orf_start: bool,
    include_peptide: bool,
}

enum Command {
    Help,
    Run(CliConfig),
}

fn run(config: CliConfig) -> Result<(), Cpc2Error> {
    let input = File::open(&config.input)?;
    let reader = BufReader::new(input);
    let options = AnalyzeOptions {
        check_reverse: config.reverse,
        include_orf_start: config.include_orf_start,
        include_peptide: config.include_peptide,
    };
    let results = analyze_reader(reader, &options)?;

    let output = File::create(&config.output)?;
    let mut writer = BufWriter::new(output);
    write_results(&mut writer, &results, &options)?;
    eprintln!("[INFO] Running Done!");
    Ok(())
}

fn parse_args<I>(args: I) -> Result<Command, Cpc2Error>
where
    I: IntoIterator<Item = String>,
{
    let mut args = args.into_iter().peekable();
    if args.peek().is_none() {
        return Ok(Command::Help);
    }

    let mut input = None;
    let mut output = None;
    let mut reverse = false;
    let mut include_orf_start = false;
    let mut include_peptide = false;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => return Ok(Command::Help),
            "-i" | "--input" => {
                let value = args
                    .next()
                    .ok_or_else(|| Cpc2Error::Args("missing value for input".to_string()))?;
                input = Some(PathBuf::from(value));
            }
            "-o" | "--output" => {
                let value = args
                    .next()
                    .ok_or_else(|| Cpc2Error::Args("missing value for output".to_string()))?;
                output = Some(resolve_output_path(value));
            }
            "-r" | "--reverse" => reverse = true,
            "--ORF" | "--orf-start" => include_orf_start = true,
            "--peptide" => include_peptide = true,
            _ => return Err(Cpc2Error::Args(format!("unknown argument: {arg}"))),
        }
    }

    let input = input.ok_or_else(|| Cpc2Error::Args("input FASTA is required".to_string()))?;
    let output = output.unwrap_or_else(|| resolve_output_path("cpc2output".to_string()));

    Ok(Command::Run(CliConfig {
        input,
        output,
        reverse,
        include_orf_start,
        include_peptide,
    }))
}

fn resolve_output_path(path: String) -> PathBuf {
    if path.ends_with(".txt") {
        PathBuf::from(path)
    } else {
        PathBuf::from(format!("{path}.txt"))
    }
}

fn print_help() {
    println!(
        "cpc2-rs\n\
         \n\
         Usage:\n\
           cargo run --release -- -i <input.fasta> [-o <output>] [-r] [--ORF] [--peptide]\n\
         \n\
         Options:\n\
           -i, --input <FILE>    Input sequence in FASTA format\n\
           -o, --output <FILE>   Output file or basename (default: cpc2output.txt)\n\
           -r, --reverse         Also check the reverse strand\n\
           --ORF, --orf-start    Add the 1-based longest ORF start column\n\
           --peptide             Add the putative peptide column\n\
           -h, --help            Show this help message"
    );
}
