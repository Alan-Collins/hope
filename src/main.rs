#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(unused_mut)]
#![allow(unreachable_code)]



use clap::Parser;
use bam;

mod homopolymer;
mod io;
mod read_alignment;


/// hope (homopolymer performance). Identify portions of long reads that map to
/// specified homopolymers in an assembly. Report errors in the sequencing of
/// those homopolymers
#[derive(Parser)]
#[clap(version = "0.2.1", author = "Alan Collins <Alan.Collins@IHRC.com>")]
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
    outprefix: String,
    /// include sequence context in outfile?
    #[clap(short, long)]
    context: bool,
}


fn main() {
    std::env::set_var("RUST_BACKTRACE", "1");
    let args = Opts::parse();
    let homos = io::read_homo_pol_file(args.input_homos);
    
    let fasta_seq = io::read_fasta(args.assembly);

    let mut outlines: Vec<String> = Vec::new();
    if args.context {
        outlines.push("homopolymer_length\thomopolymer_base\tdifference\tread_context\tassembly_context\thomo_start\tread_ID\n".to_string());
    } else {
        outlines.push("homopolymer_length\thomopolymer_base\tdifference\thomo_start\tread_ID\n".to_string());
    }

    let reader = bam::BamReader::from_path(args.bam, 0).unwrap();
    for record in reader {
        let record = record.unwrap();
        // skip if map is secondary or supplementary
        if record.flag().is_secondary() | record.flag().is_supplementary() {
            continue
        }

        // extract read name
        let name: &str = std::str::from_utf8(&record.name()).unwrap();
        
        // extract cigar string as vector of tuples
        let mut cig: Vec<(String, u32)> = Vec::new();
        for (l, c) in record.cigar().iter() {
            cig.push((c.to_string(), l));
        }

        // extract basic details
        let contig_id: i32 = record.ref_id();

        // following line doesn't work directly for some reason. Need to split in two
        // let seq: &str = std::str::from_utf8(&record.sequence().to_vec()).unwrap();

        // This seems to be the same thing, but works. 
        let temp = &record.sequence().to_vec();
        let seq: &str = std::str::from_utf8(&temp).unwrap();

        let contig = fasta_seq.seq_idxs.get(&contig_id).unwrap();
        let ref_seq = fasta_seq.seq_map.get(&contig.to_string()).unwrap();
        let start = record.start();
        let end = record.calculate_end();
        let flag = record.flag().0;

        let mut ra = read_alignment::ReadAlignment { 
            cig: cig,
            contig: contig.to_string(),
            contig_id: contig_id,
            seq: seq.to_string(),
            pos: start,
            end: end,
            aligned_end: 0,
            name: name.to_string(),
            flag: flag,
        };
        ra.aligned_end = ra.get_aligned_index(ra.end as u32) as i32;
        // ra.generate_alignment(ref_seq);

        for homo in &homos {
            // if read doesn't map to homopolymer, skip
            if ra.contig != homo.contig{
                continue
            }
            if ra.pos as u32 > homo.start {
                continue
            } 
            if homo.stop > ra.end as u32 {
                continue
            }

            let mut hr = homopolymer::HomopolymerResult::new(&homo, &ra, &ref_seq);

            let score = match hr.score {
                homopolymer::HomopolymerScore::Other(score) => score,
                homopolymer::HomopolymerScore::Difference(score) => score.to_string(),
            };

            if args.context {
                outlines.push(format!("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n", hr.homo_length, hr.base, score, &hr.region_read_aln, &hr.region_ref_aln, hr.homo.start, hr.ra.name));
            } else {
                outlines.push(format!("{0}\t{1}\t{2}\t{3}\t{4}\n", hr.homo_length, hr.base, score, hr.homo.start, hr.ra.name));
            }
            
            
        }
    }
    let outcontents = outlines.join("");
    let outfile = format!("{}{}", args.outprefix, "out.txt");
    std::fs::write(outfile, outcontents).expect("Unable to write file");
}
