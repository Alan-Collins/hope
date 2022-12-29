#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

use clap::Parser;

mod homopolymer;
mod io;


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


fn main() {
    let args = Opts::parse();
    let homos = crate::io::read_homo_pol_file(args.input_homos);
    
    let fasta = crate::io::read_fasta(args.assembly);
    

}
