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

	let refbases = refbases::get_fasta_base(&args.fasta, &args.loci);
	//println!("{:?}", refbases);
	print_header(&args.bam);

	//refbases.iter().map(|(k, v)|, {fetch_bases(&args.bam, k, v)});
	let pb = ProgressBar::new(refbases.len() as u64);

	for (k, v) in refbases{
		let kspl: Vec<&str> = k.split(':').collect();
		fetch_bases(&args.bam, kspl[0], kspl[1], &v);
		pb.inc(1);
	}
	pb.finish_with_message("done");
}

fn print_header(bam_files: &Vec<String>){
	
	print!("chr\tpos\trefbase");

	for bam in bam_files{
		let path = Path::new(bam);
		let file_stem = path.file_stem().unwrap();
		print!("\t{:?}", file_stem);
	}
	print!("\n");
}

pub fn fetch_bases(bam_files: &Vec<String>, chr: &str, pos: &str, refbase: &str){

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
			}else {
				countsdb[5] = countsdb[5]+1;
			}
		}

		pos_depth_db.insert(bam_file.to_string(), countsdb);
	}

	print!("{}\t{}\t{}", chr, pos, refbase);
	for bam_file in bam_files{
		print!("\t{:?}", pos_depth_db.get(bam_file).unwrap());
	}
	print!("\n");
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