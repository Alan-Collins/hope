use std::cmp;


#[derive(Debug)]
pub struct ReadAlignment {
    pub cig: Vec<(String, u32)>,
    pub contig: String,
    pub contig_id: i32,
    pub pos: i32,
    pub name: String,
    pub seq: String,
    pub flag: u16,
    pub whole_ref_alignment: String,
    pub whole_read_alignment: String,
}

impl ReadAlignment {
    pub fn get_aligned_index(&self, pos: u32) -> u32 {
        let mut read_idx: u32 = 0;
        let mut ref_idx: u32 = self.pos.try_into().unwrap();
        for (c, l) in &self.cig {
            if vec!["H","S"].iter().any(|&i| i==c) {
                continue
            } else if c=="D" {
            	read_idx += l;
            	ref_idx += l;
            } else if c=="M" {
            	if ref_idx + l >= pos {
            		read_idx += pos - ref_idx;
                	break
            	} else {
                	read_idx += l;
                	ref_idx += l;
                }
            } else if c=="I" {
            	if ref_idx==pos {
            		break
            	} else {
            		read_idx += l;
            	}
            } else {
                println!("unrecognized cigar");
            }
        }
        read_idx
    }

    pub fn generate_alignment(&mut self, ref_seq: &String) {
        let mut read_idx: i32 = 0;
        let mut read_seq = String::new();
        let mut ref_idx: i32 = self.pos.try_into().unwrap();
        let mut aln_ref_seq = String::new();
        for (c, l) in &self.cig {
        	let l = *l as i32;
            if vec!["H","S"].iter().any(|&i| i==c) {
                read_idx += l;
            } else if c=="D" {
            	read_seq.push_str(&"-".repeat((l).try_into().unwrap()));
            	aln_ref_seq.push_str(&ref_seq[(ref_idx).try_into().unwrap()..(ref_idx+l).try_into().unwrap()]);
            	ref_idx += l;
            } else if c=="M" {
                read_seq.push_str(&self.seq[(read_idx).try_into().unwrap()..(read_idx+l).try_into().unwrap()]);
            	read_idx += l;
            	aln_ref_seq.push_str(&ref_seq[(ref_idx).try_into().unwrap()..(ref_idx+l).try_into().unwrap()]);
            	ref_idx += l;
            } else if c=="I" {
            	read_seq.push_str(&self.seq[(read_idx).try_into().unwrap()..(read_idx+l).try_into().unwrap()]);
         	    read_idx += l;
            	aln_ref_seq.push_str(&"-".repeat((l).try_into().unwrap()));
            } else {
                println!("unrecognized cigar");
            }
        }
	    self.whole_ref_alignment.push_str(&aln_ref_seq);
	    self.whole_read_alignment.push_str(&read_seq);
    }

    pub fn print_alignment(&self, start: usize, stop: usize) {
        for i in (start..stop).step_by(60) {
            let end = cmp::min(stop, i+60);
            let mut mism = String::new();
            let read_slice = &self.whole_read_alignment[i..end];
            let ref_slice = &self.whole_ref_alignment[i..end];
            for (a,b) in read_slice.chars().zip(ref_slice.chars()) {
                if a == b {
                    mism.push_str("*");
                } else {
                    mism.push_str(" ");
                }
            }
            println!("\tRead\t{read_slice}\t{end}");
            println!("\tRef\t{ref_slice}\t{end}");
            println!("\t\t{mism}");
        }
    }
}

// def print_alignment(self, start=0, stop=None):

//         for i in range(start, stop, 60):
//             print(f"read\t{self.whole_read_alignment[i:min(i+60, stop)]}  {min(i+60, stop)}")
//             print(f"ref \t{self.whole_ref_seq[i:min(i+60, stop)]}  {min(i+60, stop)}")
//             mism = ""
//             for a,b in zip(
//                 self.whole_read_alignment[i:i+60],
//                 self.whole_ref_seq[i:min(i+60, stop)]
//                 ):
//                 if a == b:
//                     mism += "*"
//                 else:
//                     mism += " "
//             print(f"\t{mism}\n")