use clap::{App, Arg};

pub fn fetchargs() -> Args{
    let matches = App::new("basecounts").
    arg(Arg::with_name("bam").short("b").takes_value(true).multiple(true).required(true)).
    arg(Arg::with_name("fasta").short("f").takes_value(true).multiple(false).required(true)).
    arg(Arg::with_name("loci").short("l").takes_value(true).multiple(false).required(true)).
    get_matches();

    let fasta = matches.value_of("fasta").unwrap();
    let loci = matches.value_of("loci").unwrap();
    let bam:  Vec<_> = matches.values_of("bam").unwrap().collect();

    Args::new(fasta, bam, loci)
}

pub struct Args{
    pub fasta: String,
    pub bam: Vec<String>,
    pub loci: String
}

impl Args{
    pub fn new(fasta: &str, files: Vec<&str>, loci: &str) -> Args{
        let fasta = fasta.to_string();
        let loci = loci.to_string();
        let mut bam: Vec<String> = Vec::new();
        for b in files {
            bam.push(b.to_string());
        }

        Args{fasta, bam, loci}
    }
}