#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::collections::HashMap;

use crate::homopolymer::HomopolymerRecord;



pub fn read_homo_pol_file(filename: String) -> Vec<HomopolymerRecord> {

    let mut homos: Vec<HomopolymerRecord> = Vec::new();
    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                let mut bits = l.split("\t");
                let contig: String = bits.next().unwrap().to_string();
                let start: u32 = bits.next().unwrap().parse::<u32>().unwrap();
                let stop: u32 = bits.next().unwrap().parse::<u32>().unwrap();
                let base: String = bits.next().unwrap().to_string();
                let length: u32 = bits.next().unwrap().parse::<u32>().unwrap();
                homos.push(HomopolymerRecord{
                    contig,
                    start,
                    stop,
                    base,
                    length
                });
                
            }
        }
    }
    homos
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_fasta(filename: String) -> HashMap<String, String> {

    let mut fasta_map: HashMap<String, String> = HashMap::new();
    let mut header = String::new();
    let mut seq = String::new();

    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                if l.starts_with(">") {
                    if header.is_empty() {
                        header.push_str(&l[1..]);
                        seq.clear();
                    } else {
                        fasta_map.insert(header.to_string(), seq.to_string());
                        header.clear();
                        header.push_str(&l[1..]);
                        seq.clear();
                    }
                } else {
                    seq.push_str(l.as_str());
                }
            }
        }
        fasta_map.insert(header.to_string(), seq.to_string());
    }
    fasta_map
}
