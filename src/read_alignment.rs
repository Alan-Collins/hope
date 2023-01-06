
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
}

impl ReadAlignment {
    pub fn get_aligned_index(&self, pos: u32) -> u32 {
        let mut ref_idx: u32 = self.pos.try_into().unwrap();
        for (c, l) in &self.cig {
            if vec!["H","S"].iter().any(|&i| i==c) {
                continue
            } else if c=="D" {
                if ref_idx + l >= pos {
                    break
                } else {
                	ref_idx += l;
                }
            } else if c=="M" {
            	if ref_idx + l >= pos {
            		ref_idx += pos - ref_idx;
                	break
            	} else {
                	ref_idx += l;
                }
            } else if c=="I" {
            	if ref_idx==pos {
            		break
            	}
            } else {
                println!("unrecognized cigar");
            }
        }
        ref_idx
    }


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
                        
                        // is stop also in this stretch?
                        if ref_idx + l >= pos {
                            // Add last bit of sequence
                            read_seq.push_str(&"-".repeat((pos - ref_idx).try_into().unwrap()));
                            aln_ref_seq.push_str(&ref_seq[ref_idx.try_into().unwrap()..pos.try_into().unwrap()]);
                            break
                        } else {
                            // adjust ref index
                            ref_idx += l;
                            continue
                        }
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

                        // is stop also in this stretch?
                        if ref_idx + l >= pos {
                            // Add last bit of sequence
                            let remaining = pos - ref_idx;
                            read_seq.push_str(&self.seq[read_idx.try_into().unwrap()..(read_idx + remaining).try_into().unwrap()]);
                            aln_ref_seq.push_str(&ref_seq[ref_idx.try_into().unwrap()..pos.try_into().unwrap()]);
                            break
                        } else {
                            // adjust index
                            read_idx += l;
                            ref_idx += l;
                            continue
                        }
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
