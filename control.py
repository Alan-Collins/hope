#!/usr/bin/env python3

# built-ins
import sys
import argparse

# external
import pysam

# local
from hope import (
	Homopolymer,
	ReadAlignment,
	fasta_to_dict,
	read_homo_locfile,
	read_bam
	)


class Control():
	def __init__(self):
		self.contig = ""
		self.start = 0
		self.stop = 0

class ControlResult():
	def __init__(self, ra, start, stop, contig):
		self.control_length = stop-start
		self.control_start = start
		self.control_stop = stop
		self.contig = contig
		self.ra = ra
		self.start = ra.get_aligned_index(start)
		self.stop = ra.get_aligned_index(stop)

		self.read_alignment = ra.whole_read_alignment[self.start: self.stop]
		self.ref_alignment = ra.whole_ref_seq[self.start: self.stop]
		self.length = len(self.read_alignment)



def cmdline_args():

	p = argparse.ArgumentParser(
		description=""
		)
	p.add_argument(
		"-a", "--assembly",
		required=True,
		help=""
		)
	p.add_argument(
		"-b", "--bam",
		required=True,
		help=""
		)
	p.add_argument(
		"-f", "--homopolymer_file",
		required=True,
		help=""
		)
	p.add_argument(
		"-o", "--outprefix",
		required=False,
		default="./",
		help=""
		)
	p.add_argument(
		"-l", "--length",
		required=False,
		default=10,
		type=int,
		help="length of non-homopolymer sequence to use"
		)
	p.add_argument(
		"-p", "--pad",
		required=False,
		default=10,
		type=int,
		help="minimum distance from nearest homopolymer"
		)
	# p.add_argument(
	# 	"-t", "--threads",
	# 	required=False,
	# 	default=1,
	# 	type=int,
	# 	help=""
	# 	)


	return p.parse_args()


def identify_non_homo_sites(args, all_homos, assembly_dict):
	length = args.length
	pad = args.pad
	space_needed = length + 2*pad

	regions = []

	for i, homo in enumerate(all_homos):
		if homo is not all_homos[-1]:
			region_end = all_homos[i+1].start
		else:
			region_end = len(assembly_dict[homo.contig])
		if region_end < homo.stop + space_needed:
			continue

		num_regions = (region_end - homo.stop) // space_needed
		for x in range(num_regions):
			c = Control()
			c.contig = homo.contig
			c.start = homo.stop + x*space_needed + pad
			c.stop = c.start + length
			regions.append(c)


	return regions


def process_homo(homo, reads):
	homo_results = []
	for read in reads:
		if homo.start < read.pos:
			continue
		if 256 & read.flag or 2048 & read.flag:
			continue
		
		h = HomoResult(read, homo)

		if score == "skip":
			continue
		homo_results.append(h)

	return homo_results


def main(args):
	assembly_dict = fasta_to_dict(args.assembly)
	all_homos = read_homo_locfile(args.homopolymer_file)

	all_homos.sort(reverse=False, key=lambda x: x.start)

	control_regions = identify_non_homo_sites(args, all_homos, assembly_dict)
	read_dict = {} # {rname: ReadAlignment}
	for header, seq in assembly_dict.items():
		d = read_bam(args.bam, header, 0, len(seq), assembly_dict)
		read_dict = read_dict | d

	control_results = []


	with pysam.AlignmentFile(args.bam, "rb") as samfile:
		for c in control_regions:
			for read in samfile.fetch(c.contig, c.start, c.stop):
				if 256 & read.flag or 2048 & read.flag:
					continue
				control_results.append(ControlResult(read_dict[read.qname], c.start, c.stop, c.contig))

	outcontents = "ref_length\tread_length\tdifference\tread_sequence\tref_sequence\tregion_contig\tregion_start\tread_id\n"
	for c in control_results:
		outcontents += "\t".join([str(i) for i in [
			c.length,
			args.length,
			c.length-args.length,
			c.read_alignment,
			c.ref_alignment,
			c.contig,
			c.control_start,
			c.ra.name
			]]) + "\n"

	with open(f"{args.outprefix}control_out.txt", "w") as fout:
		fout.write(outcontents)



if __name__ == '__main__':
	args = cmdline_args()
	main(args)	


# def read_bam(bam, contig, start, stop, assembly_dict):
# 	read_dict = {}
# 	# samfile = pysam.AlignmentFile(bam, "rb")
# 	with pysam.AlignmentFile(bam, "rb") as samfile:
# 		for read in samfile.fetch(contig, start, stop):
# 			# if read.qname != "8362b57e-621b-4106-9abc-ff3ea2257af7":
# 			# 	continue
# 			if 256 & read.flag or 2048 & read.flag:
# 				continue
# 			read_dict[read.qname] = ReadAlignment(read, assembly_dict)

# 	return read_dict