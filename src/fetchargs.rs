use clap::{App, Arg};

pub fn fetchargs() -> Args{
    let matches = App::new("basecounts").
    about("A tool to generate nucleotide counts and VAF table from user specific loci or gentotypes").
    version("0.1.0").
    arg(Arg::with_name("bam").takes_value(true).multiple(true).required(true).index(3).help("One or more BAM files")).
    arg(Arg::with_name("loci").takes_value(true).multiple(false).required(true).index(2).help("A 2 column tsv file with chr\\tpos or a 4 column tsv with Chr\\tStart\\tRef\\tAlt")).
    arg(Arg::with_name("fasta").short("f").takes_value(true).multiple(false).required(false).help("Indexed fasta file. If provided, extracts and adds reference base to the ouput")).
    arg(Arg::with_name("placeholder").short("p").takes_value(true).multiple(false).required(false).help("Optional column id from loci file which should be included in placeholder slot of output")).
    arg(Arg::with_name("verbose").short("v").takes_value(false).multiple(false).required(false).help("Show progress bar")).
    arg(Arg::with_name("vaf").short("a").takes_value(false).multiple(false).required(false).help("Prints allelic fractions (vaf) instead of counts")).
    arg(Arg::with_name("refalt").short("r").takes_value(false).multiple(false).required(false).help("Assumes input loci are in 4 column format. i.e, Chr\\tStart\\tRef\\tAlt")).
    get_matches();

    let loci = matches.value_of("loci").unwrap();
    let bam:  Vec<_> = matches.values_of("bam").unwrap().collect();

    let mut fasta = String::new();
    if matches.is_present("fasta") == true{
        fasta = matches.value_of("fasta").unwrap().to_string();
    }

    let mut placeholder: usize = 0;
    if matches.is_present("placeholder") == true{
        placeholder = matches.value_of("placeholder").unwrap().to_string().parse().unwrap();
        placeholder = placeholder - 1;
    }

    let mut vaf: bool = false;
    if matches.is_present("vaf") == true{
        vaf = true;
    }

    let mut refalt: bool = false;
    if matches.is_present("refalt") == true{
        refalt = true;
    }

    let mut verbose: bool = false;
    if matches.is_present("verbose") == true{
        verbose = true;
    }

    Args::new(&fasta, bam, loci, vaf, refalt, verbose, placeholder)
}

pub struct Args{
    pub fasta: String,
    pub bam: Vec<String>,
    pub loci: String,
    pub vaf: bool,
    pub refalt: bool,
    pub verbose: bool,
    pub placeholder: usize,
}

impl Args{
    pub fn new(fasta: &str, files: Vec<&str>, loci: &str, vaf: bool, refalt: bool, verbose: bool, placeholder: usize) -> Args{
        let fasta = fasta.to_string();
        let loci = loci.to_string();
        let mut bam: Vec<String> = Vec::new();
        for b in files {
            bam.push(b.to_string());
        }

        Args{fasta, bam, loci, vaf, refalt, verbose, placeholder}
    }
}