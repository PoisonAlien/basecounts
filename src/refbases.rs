use std::collections::BTreeMap;
use std::io::{BufReader, BufRead, SeekFrom, Seek, Read};
use std::fs::File;
use std::process;

pub fn get_fasta_base(fa: &str, bed: &str, ph: usize) -> BTreeMap<String, String>{
	
	let mut refbase: BTreeMap<String, String> = BTreeMap::new();

	if fa.is_empty(){
		refbase = get_fasta_base_decoy(bed, ph);
		return refbase;
    }
    
	let fa = Fai::new(fa.to_string());

    let fafile = File::open(&fa.fafile).expect("Could not open the file!");
    let mut fareader = BufReader::new(fafile);

    let file = File::open(bed).expect("Cound not open the file!");
    let filer = BufReader::new(file);

    for line in filer.lines(){
        let line = line.unwrap();
        let bedspl: Vec<&str> = line.split('\t').collect();
        let chr = bedspl[0];
        let mut bedstart: u64 = bedspl[1].parse().unwrap();
        bedstart = bedstart - 1;

        let itshere: i64 = get_start_seek_pos(chr, bedstart, &fa); //position of bed start positon base in ref fasta
    	let current_pos = fareader.by_ref().seek(SeekFrom::Current(0)).expect ("Could not get current position!"); //position at buffer cursor is located
        let seekthismuch = itshere - current_pos as i64; //instead of seeking from begining, seek by this much from current cursor position 
        fareader.by_ref().seek(SeekFrom::Current(seekthismuch)).expect ("Could not get current position!") as i64;
        let bedend: u64 = bedstart + 1;

        let rlen: usize = (bedend - bedstart) as usize;		

		if rlen > 0{
			let mut fastr: String = "".to_string();

			for l in fareader.by_ref().lines(){
				let l = l.unwrap();
				fastr.push_str(&l);
				if fastr.len() >= rlen{ 
					let id: String = chr.to_string() + ":" + &bedstart.to_string();
					let mut base: String = fastr[0..rlen].to_string();
					if ph != 0{
						base = base + "/" + &bedspl[ph].to_string();
					}else{
						base = base + "/" + &"-".to_string();
					}
					

					refbase.insert(id, base);
		            break;
		        }
			}
		}
        
    }
    refbase
}

fn get_start_seek_pos(chr: &str, bedstart: u64, fa: &Fai) -> i64{

    if fa.fai.contains_key(chr) == false{
        eprintln!("Chromosome {} could not be found in reference fasta file!", chr);
        process::exit(0);
    }

	let chrstart: u64 = fa.fai[chr].0;
	let linelen: u8 = fa.fai[chr].1;
	
	let skiplines: u64 = (bedstart as f64 /linelen as f64).floor() as u64;
	let reminder = bedstart as f64 % linelen as f64;
	let bytestoskip = chrstart + skiplines + (skiplines * linelen as u64) + reminder as u64;
	
	bytestoskip as i64
}


#[derive(Debug)]
struct Fai {
	fafile: String,
	fai: BTreeMap<String, (u64, u8)>,
}

impl Fai {
	fn new(fafile: String) -> Fai{
		let mut fai: BTreeMap<String, (u64, u8)> = BTreeMap::new();

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


pub fn refalt_db(bed: &str, ph: usize) -> BTreeMap<String, String> {

	let file = File::open(bed).expect("Cound not open the file!");
	let filer = BufReader::new(file);
	
	let mut refbase: BTreeMap<String, String> = BTreeMap::new();

    for line in filer.lines(){
        let line = line.unwrap();
		let bedspl: Vec<&str> = line.split('\t').collect();
        let chr = bedspl[0];
        let mut bedstart: u64 = bedspl[1].parse().unwrap();
		bedstart = bedstart - 1;

		let id: String = chr.to_string() + ":" + &bedstart.to_string();
		
		let rbase: String = bedspl[2].parse().unwrap();
		let abase: String = bedspl[3].parse().unwrap();
		let mut base: String = rbase + "/" + &abase;

		if ph != 0{
			base = base + "/" + &bedspl[ph].to_string();
		}else{
			base = base + "/" + &"-".to_string();
		}

		refbase.insert(id, base);
	}

	//println!("{:?}", refbase);

	refbase
}

//In case fasta file not provided, generate a db with "-" as refbase
pub fn get_fasta_base_decoy(bed: &str, ph: usize) -> BTreeMap<String, String>{

	let file = File::open(bed).expect("Cound not open the file!");
	let filer = BufReader::new(file);

	let mut refbase: BTreeMap<String, String> = BTreeMap::new();

    for line in filer.lines(){
		let line = line.unwrap();
        let bedspl: Vec<&str> = line.split('\t').collect();
        let chr = bedspl[0];
        let mut bedstart: u64 = bedspl[1].parse().unwrap();
		bedstart = bedstart - 1;

		let id: String = chr.to_string() + ":" + &bedstart.to_string();
		let mut rbase: String = "-".to_string();
		if ph != 0{
			rbase = rbase + "/" + &bedspl[ph].to_string();
		}else{
			rbase = rbase + "/" + &"-".to_string();
		}
		refbase.insert(id, rbase);
	}

	refbase
}