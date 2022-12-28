#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::env;

use clap::Parser;



/// hope (homopolymer performance). Identify portions of long reads that map to
/// specified homopolymers in an assembly. Report errors in the sequencing of
/// those homopolymers
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Alan Collins <Alan.Collins@IHRC.com>")]
struct Opts {
    /// file with homopolymer locations and bases
    #[clap(short, long)]
    input_homos: String,
    /// the input assembly file
    #[clap(short, long)]
    assembly: String,
    /// the input bam file
    #[clap(short, long)]
    bam: String,
    /// the outprefix
    #[clap(short, long)]
    outrefix: String,
}


#[derive(Debug)]
struct HomopolymerRecord {
    contig: String,
    start: u32,
    stop: u32,
    base: String,
    length: u32,
}

impl HomopolymerRecord {
    fn print(&self) {
        println!("{:?}", self);
    }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}


fn main() {
    let args = Opts::parse();
  
    let mut homos: Vec<HomopolymerRecord> = Vec::new();
    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(args.input_homos) {
        // Consumes the iterator, returns an (Optional) String
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
        println!("{:?}", homos);
    }
}