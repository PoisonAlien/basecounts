use rust_htslib::bam::Record;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use crate::fetchargs;
use crate::refbases;
use std::convert::TryInto;
use std::str::from_utf8;
use std::collections::BTreeMap;
use std::path::Path;
use indicatif::ProgressBar;

pub fn bambases(){
	let args = fetchargs::fetchargs();
	
	let refbases;
	if args.refalt{
		refbases = refbases::refalt_db(&args.loci);
	}else{
		refbases = refbases::get_fasta_base(&args.fasta, &args.loci);
	}
	
	print_header(&args.bam, args.refalt);

	let pb = ProgressBar::new(refbases.len() as u64);

	for (k, v) in refbases{
		let kspl: Vec<&str> = k.split(':').collect();
		fetch_bases(&args.bam, kspl[0], kspl[1], &v, args.vaf, args.refalt);
		pb.inc(1);
	}
	pb.finish_with_message("done");
}

fn print_header(bam_files: &Vec<String>, refalt: bool){
	
	if refalt{
		print!("chr\tpos\tgenotype");
	}else{
		print!("chr\tpos\trefbase");
	}
	

	for bam in bam_files{
		let path = Path::new(bam);
		let file_stem = path.file_stem().unwrap();
		print!("\t{:?}", file_stem);
	}
	print!("\n");
}

pub fn fetch_bases(bam_files: &Vec<String>, chr: &str, pos: &str, refbase: &str, vaf: bool, ratbl: bool){

    let region: String = chr.to_owned() + ":" + pos + "-" + pos;
	let posi64: i64 = pos.parse().unwrap();

	let mut pos_depth_db:  BTreeMap<String, Vec<u32>> = BTreeMap::new();

    for bam_file in bam_files{
        let mut bam = bam::IndexedReader::from_path(bam_file).unwrap();
		bam.fetch_str(region.as_bytes()).unwrap();	
		let mut countsdb: Vec<u32> = vec!(0; 6);

        for r in bam.records(){
            let r = r.unwrap();
			let rbase = seq_at(&r, &posi64);
			
			if rbase.rcigar == "M".to_string(){
				match rbase.rbase.as_str() {
					"A" => countsdb[0] = countsdb[0]+1,
					"T" => countsdb[1] = countsdb[1]+1,
					"G" => countsdb[2] = countsdb[2]+1,
					_ => countsdb[3] = countsdb[3]+1,
				}
			} else if rbase.rcigar == "I".to_string(){
				countsdb[4] = countsdb[4]+1;
			}else if rbase.rcigar == "D".to_string(){
				countsdb[5] = countsdb[5]+1;
			}
		}

		pos_depth_db.insert(bam_file.to_string(), countsdb);
	}

	print!("{}\t{}\t{}", chr, pos, refbase);
	for bam_file in bam_files{
		for val in pos_depth_db.get(bam_file){
			print!("\t");
			if ratbl{
				print_racounts(val, refbase, vaf);
			}else{
				printcounts(val, vaf);
			}
		}
	}
	print!("\n");
}

fn print_racounts(ip: &Vec<u32>, ra: &str, vaf: bool){
	let raspl: Vec<&str> = ra.split("/").collect();
	let vartype; //String

	let mut t_depth: f32 = 0.0;
	for x in ip{
		t_depth = t_depth + *x as f32;
	}

	if raspl[0].len() > raspl[1].len(){
		vartype = String::from("Del");
	}else if raspl[0].len() < raspl[1].len(){
		vartype = String::from("Ins");
	}else{
		if raspl[0] == "-".to_string() {
			vartype = String::from("Ins");
		}else if raspl[1] == "-".to_string() {
			vartype = String::from("Del");
		}else{
			vartype = String::from("SNP");
		}
	}

	let rcount;
	if vartype == "SNP"{
		rcount = match raspl[0] {
			"A" => ip[0],
			"T" => ip[1],
			"G" => ip[2],
			_ => ip[3],
		}
	}else{
		//If var allele is not a SNP, refallele would be one of the 4 bases (with max readcounts)
		let maxval = ip[0..4].iter().max();
		rcount = match maxval {
			Some(max) => *max,
			None      => 0,
		}
	}

	let acount;
	if vartype == "SNP".to_string(){
		acount = match raspl[1] {
			"A" => ip[0],
			"T" => ip[1],
			"G" => ip[2],
			_ => ip[3],
		}
	}else if vartype == "Ins".to_string(){
		acount = ip[4]
	}else if vartype == "Del".to_string(){
		acount = ip[5]
	}else{
		panic!("Unknown variant type!")
	}

	if vaf{
		let mut vf: f32 = acount as f32 / t_depth ;
		vf = (vf * 1000.0).round() / 1000.0;
		print!("{}", vf);
	}else{
		print!("{}|{}", rcount, acount);
	}

	
}

fn printcounts(ip: &Vec<u32>, fract: bool){
	
	let mut valstr: String = String::new();
	if fract{
		let mut t_depth: f32 = 0.0;
		for x in ip{
			t_depth = t_depth + *x as f32;
		}
		
		for v in ip{
			let mut vf: f32 = *v as f32 /t_depth;
			vf = (vf * 1000.0).round() / 1000.0;
			if valstr.is_empty(){
			valstr = valstr + &vf.to_string();
			}else{
				valstr = valstr + "|" + &vf.to_string();
			}
		}
	}else{
		for v in ip{
			if valstr.is_empty(){
				valstr = valstr + &v.to_string();
			}else{
				valstr = valstr + "|" + &v.to_string();
			}
		}
	}
	
	print!("{}", valstr)
}

pub fn seq_at(rec: &Record, pos: &i64) -> Rbase{

	let pos_on_read: usize = match rec.cigar().read_pos((*pos).try_into().unwrap(), false,  true).unwrap(){
			Some(x) => x as usize,
			None =>  0 as usize,
	};

	let seq: String = from_utf8(&rec.seq().as_bytes()).unwrap().to_string();
	let rname: String = from_utf8(&rec.qname()).unwrap().to_string();
    let mut base_at_manual: String = String::from("NA");
    let mut base_at_cigar: String = String::new();
	
	let mut temp_pos: i64 = rec.pos().clone();

	for c in &rec.cigar(){
		let cig_len: i64 = c.len().try_into().unwrap();
		

		if c.char() == 'M'{
			temp_pos = temp_pos + cig_len;
			if temp_pos > *pos{
                base_at_manual =  (&seq[pos_on_read..(pos_on_read+1)]).to_string();
                base_at_cigar = "M".to_string();
				break;
			}
		}

		if c.char() == 'I'{
			if temp_pos == *pos{
				let cig_len: usize = c.len().try_into().unwrap();
                base_at_manual =  (&seq[pos_on_read..(pos_on_read+cig_len)]).to_string();
                base_at_cigar = "I".to_string();
				break;
			}
        }
        
        if c.char() == 'D'{
			if temp_pos == *pos{
				let cig_len: usize = c.len().try_into().unwrap();
                base_at_manual =  (&seq[pos_on_read..(pos_on_read+cig_len)]).to_string();
                base_at_cigar = "D".to_string();
				break;
			}
		}
	}
	Rbase::new(rname, base_at_manual, base_at_cigar)
}

#[derive(Debug)]
pub struct Rbase {
	rid: String,
	rbase: String,
	rcigar: String
}

impl Rbase {
	fn new(rid: String, rbase: String, rcigar: String) -> Rbase{
		Rbase{rid, rbase, rcigar}
	}
}