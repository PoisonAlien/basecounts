use std::fs::File;
use std::process;
use std::str::from_utf8;
use std::convert::TryInto;
use std::collections::HashMap;
use std::io::{BufReader, BufRead, SeekFrom, Seek, Read};
//use std::io::Read;

use rust_htslib::bam;
use rust_htslib::bam::Read as OtherRead;
use rust_htslib::bam::Record;

use clap::{Arg, App};

fn main() {
	let matches = App::new("basecounts").version("1.0").about("Fetch base counts for given loci from BAM files").author("Anand Mayakonda")
	.arg(Arg::with_name("loci").short("l").value_name("FILE").help("tsv/bed file chromosome and position").takes_value(true).required(true).max_values(1))
	.arg(Arg::with_name("bam").short("b").value_name("FILE").required(true).max_values(1).help("Input BAM file").takes_value(true).multiple(false))
	.arg(Arg::with_name("mapq").short("q").required(false).default_value("0").help("Mapping quality threshold (Default 30)"))
	.arg(Arg::with_name("flag").short("F").required(false).default_value("1024").help("Filter reads by SAM flags (Default 1024)"))
	.arg(Arg::with_name("one").short("o").required(false).takes_value(false).help("loci one based (Default zero based)"))
	.arg(Arg::with_name("fasta").short("f").required(false).takes_value(true).help("Indexed fasta file"))
	.get_matches();

	let loci = matches.value_of("loci").unwrap();
	let bam = matches.value_of("bam").unwrap();
	let mq = matches.value_of("mapq").unwrap();
	let flag = matches.value_of("flag").unwrap();
	let mut fasta = String::new();
	
	//println!("{:?}", fasta);
	if matches.is_present("fasta"){
		fasta = matches.value_of("fasta").unwrap().to_string();
	}
	
	let onebased = matches.is_present("one");

	let rf = Readfilt::new(&mq, &flag, onebased);

	parse_loci(&loci, &bam, &rf, &fasta);
}

fn parse_loci(loci_file: &str, bam_file: &str, rf: &Readfilt, fa: &String) {
	let file = File::open(loci_file).expect("Unable to open");
	let reader = BufReader::new(file);

	let bedstartbases = get_fasta_base(fa, loci_file);
	

	println!("Loci\tRef\tA\tT\tG\tC\tOthers\ttotal_depth");
	for line in reader.lines(){
		let line = line.unwrap();

		let line_split: Vec<&str> = line.split('\t').collect();
		let basecount = fetch_bases(bam_file, line_split[0], line_split[1], rf);
		
		let refbase: String = match bedstartbases.get(&basecount.loci){
			None => "NA".to_string(),
			Some(x) => x.to_string(),
		};

		println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", basecount.loci, refbase, &basecount.bcounts[0], &basecount.bcounts[1], &basecount.bcounts[2], &basecount.bcounts[3], &basecount.bcounts[4], &basecount.bcounts[5]);
		//println!("{}\t{}\t{:?}", , , basecount.bcounts);		
	}

}


fn summarize_loci(chr: &str, pos: &i64, hm: &HashMap<String, String>) -> Locisummary{
							 	 //A, T, G, C, INDEL
	let mut depth: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

	for v in hm.values(){
		match v.as_str() {
			"A" => depth[0] = depth[0]+1.0,
			"T" => depth[1] = depth[1]+1.0,
			"G" => depth[2] = depth[2]+1.0,
			"C" => depth[3] = depth[3]+1.0,
			_ => depth[4] = depth[4]+1.0,
		}
	}

	depth[5] = &depth[0] + &depth[1] + &depth[2] + &depth[3] + &depth[4];

	let id: String = chr.to_owned() + ":" + &pos.to_owned().to_string();
	Locisummary::new(&id, &depth)
	//println!("{:?}", tdepth);

	//println!("A\tT\tG\tC\tINDELS");
	//println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", chr, pos, tdepth, &depth[0], &depth[1], &depth[2], &depth[3], &depth[4]);
	//println!("{}\t{}\t{}({:.2})\t{}({:.2})\t{}({:.2})\t{}({:.2})\t{}({:.2})", chr, pos, &depth[0], (&depth[0]/tdepth), 
	//	&depth[1], &depth[1]/tdepth, &depth[2], &depth[2]/tdepth, &depth[3], &depth[3]/tdepth, &depth[4], &depth[4]/tdepth);

}

struct Locisummary {
	loci: String,
	bcounts: Vec<f32>,
}


impl Locisummary {
	fn new(chr: &String, counts: &Vec<f32>) -> Locisummary {
		let loci: String = chr.to_owned().to_string();
		let bcounts: Vec<f32> = counts.clone();
		Locisummary{loci, bcounts}
	}
}


fn fetch_bases(bam_file: &str, chr: &str, pos: &str, rf: &Readfilt) -> Locisummary{
	let mut posi64: i64 = pos.parse().unwrap();
	//println!("{:?}", pos);
	let mut bam = bam::IndexedReader::from_path(bam_file).unwrap();
	let mut region: String = chr.to_owned() + ":" + pos + "-" + pos;
	
	if rf.onebased == true{
		posi64 = posi64 + 1;
		let region: String = chr.to_owned() + ":" + &posi64.to_string() + "-" + &posi64.to_string();
	}
	

	let mut rinfo: HashMap<String, String> = HashMap::new();

	bam.fetch_str(region.as_bytes()).unwrap();	
	for r in bam.records(){
		let r = r.unwrap();
		if r.mapq() >= rf.mq {
			let base = seq_at(&r, &posi64); 	
			rinfo.insert(base.rid, base.base);
		}
	}

	summarize_loci(&chr, &posi64, &rinfo)
}

fn seq_at(rec: &Record, pos: &i64) -> Rbase {

	let mut pos_on_read: usize = 0;
	match rec.cigar().read_pos((*pos).try_into().unwrap(), false,  true).unwrap(){
			Some(x) => pos_on_read = x as usize,
			None => pos_on_read = 0,
	}

	let seq: String = from_utf8(&rec.seq().as_bytes()).unwrap().to_string();
	let rname: String = from_utf8(&rec.qname()).unwrap().to_string();
	//eprintln!("{:?}", from_utf8(&r.qname()).unwrap().to_string());
	//let mut base_at: String = (&seq[pos_on_read..(pos_on_read+1)]).to_string();
	let mut base_at_manual: String = String::from("NA");
	
	
	let mut temp_pos: i64 = rec.pos().clone();

	for c in &rec.cigar(){
		let cig_len: i64 = c.len().try_into().unwrap();
		

		if c.char() == 'M'{
			temp_pos = temp_pos + cig_len;
			if temp_pos > *pos{
				base_at_manual =  (&seq[pos_on_read..(pos_on_read+1)]).to_string();
				break;
			}
		}

		if c.char() == 'I'{
			if temp_pos == *pos{
				let cig_len: usize = c.len().try_into().unwrap();
				base_at_manual =  (&seq[pos_on_read..(pos_on_read+cig_len)]).to_string();
				break;
			}
		}
	}
	Rbase::new(rname, base_at_manual)
}

#[derive(Debug)]
struct Rbase {
	rid: String,
	base: String,
}

impl Rbase {
	fn new(rid: String, rbase: String) -> Rbase{
		let mut base = String::new();
		match rbase.as_str() {
			"A" | "T" | "G" | "C" => base = rbase,
			_ => base = "INDEL".to_string(),
		}
		Rbase{rid, base}
	}
}

#[derive(Debug)]
struct Readfilt {
	mq: u8,
	flag: u16,
	onebased: bool,
}

impl Readfilt {
	fn new(mq: &str, flag: &str, onebased: bool) -> Readfilt{
		let mq: u8 = mq.clone().parse().unwrap();
		let flag: u16 = flag.clone().parse().unwrap();
		let onebased: bool = onebased;
		Readfilt{mq, flag, onebased}
	}
}


#[derive(Debug)]
struct Fai {
	fafile: String,
	fai: HashMap<String, (u64, u8)>,
}

impl Fai {
	fn new(fafile: String) -> Fai{
		let mut fai: HashMap<String, (u64, u8)> = HashMap::new();

		let faifile: String = fafile.to_owned() + &".fai".to_string();
		let faifile = File::open(faifile).expect("Unable to open the fai file");
		let faireader = BufReader::new(faifile);

		for line in faireader.lines(){
			let line = line.unwrap();
			let line_spl: Vec<&str> = line.split('\t').collect();
			let chr: String = line_spl[0].to_string();
			let offset: u64 = line_spl[2].parse().unwrap();
			let linelen: u8 = line_spl[3].parse().unwrap();
			fai.insert(chr, (offset, linelen));
		}
		Fai{fafile, fai}
	}
}

fn get_fasta_base(fa: &str, bed: &str) -> HashMap<String, String>{
	
	let mut refbase: HashMap<String, String> = HashMap::new();

	if fa.is_empty(){
		return refbase;
	}
	
	let fa = Fai::new(fa.to_string());

    let mut fafile = File::open(&fa.fafile).expect("Nope!");
    let mut fareader = BufReader::new(fafile);

    let bedfile = File::open(bed).expect("Nope!");
    let mut bedreader = BufReader::new(bedfile);

    

    for bedlines in bedreader.by_ref().lines(){
    	let bedlines = bedlines.unwrap();
    	let bedspl: Vec<&str> = bedlines.split('\t').collect();
    	let itshere: i64 = get_start_seek_pos(&bedspl, &fa); //position of bed start positon base in ref fasta
    	let current_pos = fareader.by_ref().seek(SeekFrom::Current(0)).expect ("Could not get current position!"); //position at buffer cursor is located
    	let seekthismuch = itshere - current_pos as i64; //instead of seeking from begining, seek by this much from current cursor position 
    	fareader.by_ref().seek(SeekFrom::Current(seekthismuch)).expect ("Could not get current position!") as i64;
    	//get_seq(bedspl, &fareader);
    	//println!("{}\t{}\t{}\t{}", itshere, current_pos, current_pos2, seekthismuch);
    	let bedstart: u64 = bedspl[1].to_owned().parse().unwrap();
		//let bedend: u64 = bedspl[2].to_owned().parse().unwrap();
		let bedend: u64 = bedstart + 1;
		let rlen: usize = (bedend - bedstart) as usize;

		//println!("{}\t{}\t{}", bedstart, bedend, rlen);

		if rlen > 0{
			let mut fastr: String = "".to_string();

			for l in fareader.by_ref().lines(){
				let l = l.unwrap();
				fastr.push_str(&l);
				if fastr.len() >= rlen{ 
					let id: String = bedspl[0].to_string() + ":" + bedspl[1];
					let base: String = fastr[0..rlen].to_string();
					refbase.insert(id, base);
					//println!(">{}:{}-{}", bedspl[0], bedstart, bedend);
					//println!("{}", &fastr[0..rlen]);
					//prettyprint(&fastr[0..rlen], lwd);
					//println!("{:?}", buffpos+linetracker);
		            break;
		        	}
				}
			}
    	}
    	refbase
	}
	

fn get_start_seek_pos(bedline: &Vec<&str>, fa: &Fai) -> i64{

	let bedstart: u64 = bedline[1].parse().unwrap();
	let bedend: u64 = bedstart + 1;
	let rlen: usize = (bedend - bedstart) as usize;
	let chrstart: u64 = fa.fai[bedline[0]].0;
	let linelen: u8 = fa.fai[bedline[0]].1;
	

	let skiplines: u64 = (bedstart as f64 /linelen as f64).floor() as u64;
	let reminder = bedstart as f64 % linelen as f64;
	let bytestoskip = chrstart + skiplines + (skiplines * linelen as u64) + reminder as u64;
	
	bytestoskip as i64

}