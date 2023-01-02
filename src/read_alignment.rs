use std::cmp;


#[derive(Debug)]
pub struct ReadAlignment {
    pub cig: Vec<(String, u32)>,
    pub contig: String,
    pub contig_id: i32,
    pub pos: i32,
    pub end: i32,
    pub aligned_end: i32,
    pub name: String,
    pub seq: String,
    pub flag: u16,
    // pub whole_ref_alignment: String,
    // pub whole_read_alignment: String,
}

impl ReadAlignment {
    pub fn get_aligned_index(&self, pos: u32) -> u32 {
        let mut read_idx: u32 = 0;
        let mut ref_idx: u32 = self.pos.try_into().unwrap();
        for (c, l) in &self.cig {
            if vec!["H","S"].iter().any(|&i| i==c) {
                continue
            } else if c=="D" {
                if ref_idx + l >= pos {
                    read_idx += pos - ref_idx;
                    break
                } else {
                	read_idx += l;
                	ref_idx += l;
                }
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

    // pub fn generate_alignment(&mut self, ref_seq: &String) {
    //     let mut read_idx: i32 = 0;
    //     let mut read_seq = String::new();
    //     let mut ref_idx: i32 = self.pos.try_into().unwrap();
    //     let mut aln_ref_seq = String::new();
    //     for (c, l) in &self.cig {
    //     	let l = *l as i32;
    //         if vec!["H","S"].iter().any(|&i| i==c) {
    //             read_idx += l;
    //         } else if c=="D" {
    //         	read_seq.push_str(&"-".repeat((l).try_into().unwrap()));
    //         	aln_ref_seq.push_str(&ref_seq[(ref_idx).try_into().unwrap()..(ref_idx+l).try_into().unwrap()]);
    //         	ref_idx += l;
    //         } else if c=="M" {
    //             read_seq.push_str(&self.seq[(read_idx).try_into().unwrap()..(read_idx+l).try_into().unwrap()]);
    //         	read_idx += l;
    //         	aln_ref_seq.push_str(&ref_seq[(ref_idx).try_into().unwrap()..(ref_idx+l).try_into().unwrap()]);
    //         	ref_idx += l;
    //         } else if c=="I" {
    //         	read_seq.push_str(&self.seq[(read_idx).try_into().unwrap()..(read_idx+l).try_into().unwrap()]);
    //      	    read_idx += l;
    //         	aln_ref_seq.push_str(&"-".repeat((l).try_into().unwrap()));
    //         } else {
    //             println!("unrecognized cigar");
    //         }
    //     }
	//     self.whole_ref_alignment.push_str(&aln_ref_seq);
	//     self.whole_read_alignment.push_str(&read_seq);
    // }

    // pub fn print_alignment(&self, start: usize, stop: usize) {
    //     for i in (start..stop).step_by(60) {
    //         let end = cmp::min(stop, i+60);
    //         let mut mism = String::new();
    //         let read_slice = &self.whole_read_alignment[i..end];
    //         let ref_slice = &self.whole_ref_alignment[i..end];
    //         for (a,b) in read_slice.chars().zip(ref_slice.chars()) {
    //             if a == b {
    //                 mism.push_str("*");
    //             } else {
    //                 mism.push_str(" ");
    //             }
    //         }
    //         println!("\tRead\t{read_slice}\t{end}");
    //         println!("\tRef\t{ref_slice}\t{end}");
    //         println!("\t\t{mism}");
    //     }
    // }

    pub fn extract_alignment(&self, start: u32, stop: u32, ref_seq: &String) -> (String, String) {
        let mut read_idx: u32 = 0;
        let mut read_seq = String::new();
        let mut ref_idx: u32 = self.pos.try_into().unwrap();
        let mut aln_ref_seq = String::new();
        // keep track of whether we are finding start or stop
        let mut on_start = true;
        let mut on_stop = false;
        let mut pos = start;
        for (c, l) in &self.cig {
            println!("{:?}", (c, l));
            if vec!["H","S"].iter().any(|&i| i==c) {
                read_idx += l;
                continue
            } else if c=="D" {
                if ref_idx + l >= pos {
                    if on_start {
                        // Add first bit of sequence
                        read_seq.push_str(&"-".repeat((ref_idx + l - pos).try_into().unwrap()));
                        aln_ref_seq.push_str(&ref_seq[pos.try_into().unwrap()..(ref_idx + l).try_into().unwrap()]);
                        // now looking for end
                        on_start = false;
                        on_stop = true;
                        pos = stop;
                        // adjust ref index
                        ref_idx += l ;
                        continue
                    } else if on_stop {
                        // Add last bit of sequence
                        read_seq.push_str(&"-".repeat((pos - ref_idx).try_into().unwrap()));
                        aln_ref_seq.push_str(&ref_seq[ref_idx.try_into().unwrap()..pos.try_into().unwrap()]);
                        break
                    }
                } else if on_stop {
                    // Add sequence
                    read_seq.push_str(&"-".repeat((*l).try_into().unwrap()));
                    aln_ref_seq.push_str(&ref_seq[ref_idx.try_into().unwrap()..(ref_idx + l).try_into().unwrap()]);
                    ref_idx += l;
                } else {
                    ref_idx += l;
                }
            } else if c=="M" {
                if ref_idx + l >= pos {
                    if on_start {
                        // Add first bit of sequence
                        let intermediate_index = pos - ref_idx;
                        read_seq.push_str(&self.seq[(read_idx + intermediate_index).try_into().unwrap()..(read_idx + l).try_into().unwrap()]);
                        aln_ref_seq.push_str(&ref_seq[pos.try_into().unwrap()..(ref_idx + l).try_into().unwrap()]);
                        // now looking for end
                        on_start = false;
                        on_stop = true;
                        pos = stop;
                        // adjust index
                        read_idx += l;
                        ref_idx += l;
                        continue
                    } else if on_stop {
                        // Add last bit of sequence
                        let remaining = pos - ref_idx;
                        read_seq.push_str(&self.seq[read_idx.try_into().unwrap()..(read_idx + remaining).try_into().unwrap()]);
                        aln_ref_seq.push_str(&ref_seq[ref_idx.try_into().unwrap()..pos.try_into().unwrap()]);
                        break
                    }
                } else if on_stop {
                    // Add sequence
                    read_seq.push_str(&self.seq[read_idx.try_into().unwrap()..(read_idx + l).try_into().unwrap()]);
                    aln_ref_seq.push_str(&ref_seq[ref_idx.try_into().unwrap()..(ref_idx + l).try_into().unwrap()]);
                }
                read_idx += l;
                ref_idx += l;
            } else if c=="I" {
                if ref_idx == pos {
                    if on_stop {
                        // Add last bit of sequence
                        read_seq.push_str(&self.seq[read_idx.try_into().unwrap()..(read_idx + l).try_into().unwrap()]);
                        aln_ref_seq.push_str(&"-".repeat((*l).try_into().unwrap()));
                        break
                    }
                } else if on_stop {
                    // Add sequence 
                    read_seq.push_str(&self.seq[read_idx.try_into().unwrap()..(read_idx + l).try_into().unwrap()]);
                    aln_ref_seq.push_str(&"-".repeat((*l).try_into().unwrap()));
                    read_idx += l;
                } else {
                    read_idx += l;
                }
            } else {
                println!("unrecognized cigar");
            }
        }
        (read_seq, aln_ref_seq)
    }
}
