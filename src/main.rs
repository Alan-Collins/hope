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
    outprefix: String,
}


fn main() {
    let args = Opts::parse();
    let homos = io::read_homo_pol_file(args.input_homos);
    
    let fasta_seq = io::read_fasta(args.assembly);

    let mut outcontents = "homopolymer_length\thomopolymer_base\tdifference\tread_context\tassembly_context\thomo_start\tread_ID\n".to_string();

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

        let Some(contig) = fasta_seq.seq_idxs.get(&contig_id) else { continue };
        let Some(ref_seq) = fasta_seq.seq_map.get(&contig.to_string()) else { continue };
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
            name: name.to_string(),
            flag: flag,
            whole_read_alignment: "".to_string(),
            whole_ref_alignment: "".to_string(),
        };
        ra.generate_alignment(ref_seq);

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
            
            let start = ra.get_aligned_index(homo.start) as usize;
            let stop = ra.get_aligned_index(homo.stop) as usize;
            let upstart = std::cmp::max(start, 30) - 30;
            let downstop = std::cmp::min(stop, ra.whole_read_alignment.len()-30) + 30;
            let mut hr = homopolymer::HomopolymerResult {
                base: homo.base.to_string(),
                homo_length: homo.length,
                homo: homopolymer::HomopolymerRecord{
                    contig: homo.contig.clone(),
                    start: homo.start,
                    stop: homo.stop,
                    base: homo.base.clone(),
                    length: homo.length,
                },
                ra: &ra,
                start: start,
                stop: stop,
                read_alignment: ra.whole_read_alignment[start..stop].to_string(),
                ref_alignment: ra.whole_ref_alignment[start..stop].to_string(),
                read_upstream: ra.whole_read_alignment[upstart..start].to_string(),
                read_downstream: ra.whole_read_alignment[stop..downstop].to_string(),
                ref_upstream: ra.whole_ref_alignment[upstart..start].to_string(),
                ref_downstream: ra.whole_ref_alignment[stop..downstop].to_string(),
                length: (stop-start) as u32,
                score: homopolymer::HomopolymerScore::Difference(0), 
            };
            hr.score();

            let score = match hr.score {
                homopolymer::HomopolymerScore::Other(score) => score,
                homopolymer::HomopolymerScore::Difference(score) => score.to_string(),
            };
            let line = format!("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n", hr.homo_length, hr.base, score, format!("{}{}{}", &hr.read_upstream[std::cmp::max(hr.read_upstream.len() as isize -5, 0) as usize..], hr.read_alignment, &hr.read_downstream[..std::cmp::min(hr.read_downstream.len() as usize,5)]), format!("{}{}{}", &hr.ref_upstream[std::cmp::max(hr.ref_upstream.len() as isize -5, 0) as usize..], hr.ref_alignment, &hr.ref_downstream[..std::cmp::min(hr.ref_downstream.len() as usize,5)]), hr.homo.start, hr.ra.name);
            
            outcontents = format!("{}{}", &outcontents, &line);
        }
    }

    let outfile = format!("{}{}", args.outprefix, "out.txt");
    std::fs::write(outfile, outcontents).expect("Unable to write file");
}
