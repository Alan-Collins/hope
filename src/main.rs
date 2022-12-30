#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(unused_mut)]
#![allow(unreachable_code)]



use clap::Parser;

mod homopolymer;
mod io;
mod read_alignment;


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
    
    let fasta_seq = crate::io::read_fasta(args.assembly);

    let reads = io::read_bam(args.bam, fasta_seq);

    let mut results: Vec::<homopolymer::HomopolymerResult> = Vec::new();
    // let mut x = 0;
    for homo in homos {
        // x += 1;
        // if x > 20 {
        //     break
        // }
        for ra in reads.values() {
            // if read doesn't map to homopolymer, skip
            if ra.pos as u32 > homo.start {
                continue
            } 
            if homo.stop > ra.end as u32 {
                println!("{:?}", (ra.end, homo.stop));
                continue
            }

            let hr = homopolymer::HomopolymerResult::new(&homo, &ra);
            if hr.score == homopolymer::HomopolymerScore::Other("?".to_string()) {
                println!("{:?}", hr.score);
                hr.ra.print_alignment(hr.start-10, hr.stop+10);
                println!("{:?}", hr.read_upstream);
            } else {
                println!("{:?}", hr.score);
                hr.ra.print_alignment(hr.start, hr.stop);
            }

        }
        std::process::exit(1);
    }
    

}
