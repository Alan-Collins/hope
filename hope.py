#!/usr/bin/env python3

# /mnt/g/SynologyDrive/IHRC/work/ABiL/homopolymer
# ./hope.py -a data/Sterne_pilon.fasta -b maptest/sterne0_subseq_A2B10.bam -f data/sterne_5n.txt -o test_out

import sys
import os
import argparse
import multiprocessing
import pysam

import time
from datetime import timedelta


class Homopolymer():
	def __init__(self, line):
		bits = line.split()
		# contig	start	stop	base	length
		self.contig = bits[0]
		self.start = int(bits[1])-1
		self.stop = int(bits[2])-1
		self.base = bits[3]
		self.length = int(bits[4])

class HomoResult():
	def __init__(self, ra, homo):
		self.base = homo.base
		self.homo_length = homo.length
		self.homo = homo
		self.ra = ra
		self.start = ra.get_aligned_index(homo.start)
		self.stop = ra.get_aligned_index(homo.stop)+1
		# print(ra.print_alignment(self.start-5, self.stop+5))
		# print(ra.cigartuples[108:113])
		# sys.exit()
		self.read_alignment = ra.whole_read_alignment[self.start: self.stop]
		self.ref_alignment = ra.whole_ref_seq[self.start: self.stop]
		self.read_upstream = ra.whole_read_alignment[max(0, self.start-30): self.start]
		self.read_downstream = ra.whole_read_alignment[self.stop: min(self.stop+30, len(ra.whole_read_alignment))]
		self.ref_upstream = ra.whole_ref_seq[max(0, self.start-30): self.start]
		self.ref_downstream = ra.whole_ref_seq[self.stop: min(self.stop+30, len(ra.whole_ref_seq))]
		self.length = len(self.read_alignment)
		self.score = None

	def score_homo(self):
		# Takes a ReadAlignment that maps to homopolymer site
		base = self.base
		# Before scoring, check if we have upstream and downstreak sequence
		if self.start == 0 or self.stop >= len(self.ra.whole_read_alignment):
			# If not don't include this one as we can't confidently assess
			# homopolymer boundaries
			self.score = "skip"
			return self.score
		
		# First check for identical homopolymer with no flanking gaps
		if (all([i != "-" for i in self.ref_alignment])
			and all([i != "-" for i in self.read_alignment])
			and (self.ref_upstream[-1] not in [base, "-"]) 
			# and (self.ref_upstream[-1] == self.read_upstream[-1])
			and (self.ref_downstream[0] not in [base, "-"])
			# and (self.ref_downstream[0] == self.read_downstream[0])
			):
			if not all([i == base for i in self.read_alignment]):
				# non-homopolymer bases found
				self.score = "mm"
				return self.score
			
			else:
				self.score = 0
				return self.score

		# next handle identical homopolymer with flanking gaps in read
		if (all([i != "-" for i in self.read_alignment])
			and ((self.read_upstream[-1] == "-") 
				or (self.read_downstream[0] == "-"))):
			if self.read_upstream[-1] == "-":
				i = -1
				while self.read_upstream[i] == "-":
					if self.ref_upstream[i] == base:
						self.score = "?"
						return self.score

					i -= 1

				if self.read_upstream[i] == base:
					# indel flanked by homopolymer base. 
					# Call it homopolymer-associated error
					self.score = "?"
					return self.score

			if self.read_downstream[0] == "-":
				i = 0
				while self.read_downstream[i] == "-":
					if self.ref_downstream[i] == base:
						self.score = "?"
						return self.score

					i += 1

				if self.read_downstream[i] == base:
					# indel flanked by homopolymer base. 
					# Call it homopolymer-associated error
					self.score = "?"
					return self.score
			self.score = 0
			return self.score

		# next handle extension of homopolymer in read
		if any([i == "-" for i in self.ref_alignment]):
			# Check if non-homopolymer bases are the majority
			if len([i for i in self.read_alignment if i not in [base, "-"]]) > self.length / 2:
				self.score = "?"
			else:
				self.score = self.length - self.homo_length
			return self.score
		
		# next handle deletions in homopolymer
		if any([i == "-" for i in self.read_alignment]):
			# if any bases not the homopolymer base or gap, return "?"
			if any([i not in [base, "-"] for i in self.read_alignment]):
				self.score = "?"
				return self.score

			# If not flanked by gaps in read, simply truncated homopolymer
			if self.read_upstream[-1] != "-" and self.read_downstream[0] != "-":
				# Check if majority of non-gap sequence not homopolymer base
				if (len([i for i in self.read_alignment if i not in [base, "-"]])
					> len([i for i in self.read_alignment if i  == base])):
					self.score = "?"
				else:
					self.score = len([i for i in self.read_alignment if i == base]) - self.homo_length
				return self.score

			# else, flanking deletion includes non-homopolymer base, return ?
			self.score = "?"
			return self.score

		# Next handle insertions in read next to homopolymer
		if self.ref_upstream[-1] == "-" or self.ref_downstream[0] == "-":
			# Check if any inserted bases in the read are the homopolymer base
			# If so, return ?
			i = -1
			if self.ref_upstream[i] == "-":
				while self.ref_upstream[i] == "-":
					if self.read_upstream[i] == base:
						self.score = "?"
						return self.score

					i -= 1
					if -1*i == len(self.ref_downstream):
						self.score = "?"
						return self.score

			i = 0
			if self.ref_downstream[i] == "-":
				while self.ref_downstream[i] == "-":
					if self.read_downstream[i] == base:
						self.score = "?"
						return self.score

					i += 1
					if i == len(self.ref_downstream):
						self.score = "?"
						return self.score
			self.score = 0
			return self.score

class ReadAlignment():
	def __init__(self, read, ref_seq_dict):
		self.whole_ref_seq, self.whole_read_alignment = self.generate_alignment(
			read, ref_seq_dict[read.reference_name])
		self.cigartuples = read.cigartuples
		self.pos = read.pos
		self.name = read.qname
		self.flag = read.flag

	def generate_alignment(self, read, ref_seq):
		# cigar codes:
		# 0 = M
		# 1 = I
		# 2 = D
		# 3 = H???
		# 4 = S
		read_idx = 0
		read_seq = ''
		ref_idx = read.pos
		aln_ref_seq = ''
		for cig, ln in read.cigartuples:
			if cig == 1:
				read_seq += read.seq[read_idx: read_idx+ln]
				read_idx += ln
				aln_ref_seq += "-"*ln
			if cig in [3, 4]:
				read_idx += ln
			elif cig == 2:
				read_seq += '-'*ln
				aln_ref_seq += ref_seq[ref_idx: ref_idx+ln]
				ref_idx += ln
			elif cig == 0:
				read_seq += read.seq[read_idx: read_idx+ln]
				read_idx += ln
				aln_ref_seq += ref_seq[ref_idx: ref_idx+ln]
				ref_idx += ln
		return aln_ref_seq, read_seq


	def get_aligned_index(self, index):
		# cigar codes:
		# 0 = M
		# 1 = I
		# 2 = D
		# 3 = H???
		# 4 = S
		read_idx = 0
		ref_idx = self.pos

		for n, (cig, ln) in enumerate(self.cigartuples):
			if cig in [3, 4]:
				continue
			elif cig  == 1:
				if ref_idx == index:
					# If insertion at homopolymer start
					break
				else:
					read_idx += ln
			elif cig == 2:
				read_idx += ln
				ref_idx += ln
			elif cig == 0:
				if ref_idx + ln >= index:
					# if index in this match portion add remining
					# distance to read index
					read_idx += index - ref_idx
					break
				else:
					read_idx += ln
					ref_idx += ln
		# print(n)
		return read_idx


	def print_alignment(self, start=0, stop=None):
		if stop == None:
			stop = len(self.whole_read_alignment)

		for i in range(start, stop, 60):
			print(f"read\t{self.whole_read_alignment[i:min(i+60, stop)]}  {min(i+60, stop)}")
			print(f"ref \t{self.whole_ref_seq[i:min(i+60, stop)]}  {min(i+60, stop)}")
			mism = ""
			for a,b in zip(
				self.whole_read_alignment[i:i+60],
				self.whole_ref_seq[i:min(i+60, stop)]
				):
				if a == b:
					mism += "*"
				else:
					mism += " "
			print(f"\t{mism}\n")

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
		"-t", "--threads",
		required=False,
		default=1,
		type=int,
		help=""
		)
	p.add_argument(
		"fastqs", nargs="*", 
		help="reads to process"
		)

	return p.parse_args()


def fasta_to_dict(FASTA_file):
	"""Read a fasta file into a dict 

	Dict has headers (minus the > symbol) as keys and the associated 
	sequence as values.
	
	Args:
	  FASTA_file (str): 
		path to fasta format file

	Returns:
	  dict: 
		dict of format {fasta_header : sequence}

	Raises:
	  TypeError: If FASTA_file is not a str
	  OSError: If FASTA_file is not the path to an existing file
	"""
	
	if type(FASTA_file) is not str:
		raise TypeError(
			"FASTA_file must be str, not {}.".format(type(FASTA_file).__name__))

	if not os.path.exists(FASTA_file):
		raise OSError(
			"FASTA_file must be the path to an existing file.")


	fasta_dict = {}
	with open(FASTA_file, 'r') as f:
		multifasta = f.read()
	f.close()
	fastas = multifasta.split(">")
	trimmed_fastas = []
	for i in fastas:
		if len(i) != 0:
			trimmed_fastas.append(i)

	fastas = trimmed_fastas

	for i in fastas:
		header = i.split("\n")[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq

	return fasta_dict


def read_homo_locfile(file):
	homos = []
	with open(file) as f:
		for line in f:
			homos.append(Homopolymer(line))

	return homos


def process_homo(homo, reads):
	homo_results = []
	for read in reads:
		if homo.start < read.pos:
			continue
		if 256 & read.flag or 2048 & read.flag:
			continue
		
		h = HomoResult(read, homo)

		score = h.score_homo()

		if score == "skip":
			continue
		homo_results.append(h)

	return homo_results


def read_bam(bam, contig, start, stop, assembly_dict):
	read_dict = {}
	# samfile = pysam.AlignmentFile(bam, "rb")
	with pysam.AlignmentFile(bam, "rb") as samfile:
		for read in samfile.fetch(contig, start, stop):
			# if read.qname != "8362b57e-621b-4106-9abc-ff3ea2257af7":
			# 	continue
			if 256 & read.flag or 2048 & read.flag:
				continue
			read_dict[read.qname] = ReadAlignment(read, assembly_dict)

	return read_dict


def assemble_homo_list(homo, bam, read_dict):
	reads = []
	with pysam.AlignmentFile(bam, "rb") as samfile:
		# samfile = pysam.AlignmentFile(bam, "rb")
		for read in samfile.fetch(homo.contig, homo.start, homo.stop):
			if 256 & read.flag or 2048 & read.flag:
				continue
			# if read.qname != "8362b57e-621b-4106-9abc-ff3ea2257af7":
			# 	continue
			reads.append(read_dict[read.qname])

	return (homo, reads)


def main(args):
	start_time = time.time()
	assembly_dict = fasta_to_dict(args.assembly)
	all_homos = read_homo_locfile(args.homopolymer_file)
	pool = multiprocessing.Pool(processes=args.threads)

	bam_read_options = []
	for contig, seq in assembly_dict.items():
		interval = len(seq)//args.threads
		for start in range(0, len(seq), interval):
			stop = start + interval
			bam_read_options.append((args.bam, contig, start, stop, assembly_dict))

	chunksize = len(bam_read_options)//args.threads
	bam_read_results = pool.starmap(read_bam, bam_read_options, chunksize)
	pool.close()
	pool.join()
	pool.terminate()

	read_dict = {} # {rname: ReadAlignment}
	for d in bam_read_results:
		read_dict = read_dict | d

	pool = multiprocessing.Pool(processes=args.threads)
	chunksize = len(all_homos)//args.threads
	read_subset_options = [(homo, args.bam, read_dict) for homo in all_homos]
	chunksize = len(read_subset_options)//args.threads
	homo_process_options = pool.starmap(assemble_homo_list, read_subset_options, chunksize)
	pool.close()
	pool.join()
	pool.terminate()

	pool = multiprocessing.Pool(processes=args.threads)
	homo_results = []
	homo_results += pool.starmap(process_homo, homo_process_options, chunksize)
	pool.close()
	pool.join()
	pool.terminate()

	outcontents = "homopolymer_length\thomopolymer_base\tdifference\tread_context\tassembly_context\thomo_start\tread_ID\n"
	for h_list in homo_results:
		for h in h_list:
			outcontents += (f"{h.homo_length}\t{h.base}\t{h.score}\t"
				+ f"{h.read_upstream[-5:]+h.read_alignment+h.read_downstream[:5]}"
				+ f"\t{h.ref_upstream[-5:]+h.ref_alignment+h.ref_downstream[:5]}"
				+ f"\t{h.homo.start}\t{h.ra.name}\n")
	with open(f"{args.outprefix}out.txt", "w") as fout:
		fout.write(outcontents)


	end_time = time.time()
	time_taken = timedelta(seconds=end_time-start_time)
	sys.stderr.write("\nTotal run time: {}\n".format(time_taken))


if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
