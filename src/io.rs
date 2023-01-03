#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(unreachable_code)]

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::collections::HashMap;

use bam;

use crate::homopolymer::HomopolymerRecord;
use crate::read_alignment::ReadAlignment;

#[derive(Debug)]
pub struct FastaSequence {
    pub seq_idxs: HashMap<i32, String>,
    pub seq_map: HashMap<String, String>
}


pub fn read_homo_pol_file(filename: String) -> Vec<HomopolymerRecord> {

    let mut homos: Vec<HomopolymerRecord> = Vec::new();
    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                let mut bits = l.split("\t");
                let contig: String = bits.next().unwrap().to_string();
                let start: u32 = bits.next().unwrap().parse::<u32>().unwrap()-1;
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

pub fn read_fasta(filename: String) -> FastaSequence {

    let mut fasta_map: HashMap<String, String> = HashMap::new();
    let mut seq_idxs: HashMap<i32, String> = HashMap::new();
    let mut header = String::new();
    let mut seq = String::new();
    let mut seq_id: i32 = 0;

    if let Ok(lines) = read_lines(filename) {
        for line in lines {
            if let Ok(l) = line {
                if l.starts_with(">") {
                    if header.is_empty() {
                        header.push_str(&l[1..]);
                        seq.clear();
                    } else {
                        fasta_map.insert(header.to_string(), seq.to_string());
                        seq_idxs.insert(seq_id, header.to_string());
                        seq_id += 1;
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
        seq_idxs.insert(seq_id, header.to_string());
    }
    FastaSequence {
        seq_map: fasta_map,
        seq_idxs: seq_idxs
    }
}
