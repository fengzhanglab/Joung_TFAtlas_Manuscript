#counts TF barcodes from fastq for TF screening with 10X scRNA-seq readout

from Bio import SeqIO
import csv
from collections import OrderedDict
import numpy as np
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import pairwise2

KEY_REGION_START = 25 #start index of key region
KEY_REGION_END = 50 #end index of key region
BARCODE_LENGTH = 24
MAX_BARCODE_DIFF = 1 #allows 1 mismatch for barcode mapping
KEY = "GAAAGGACGA" #identifies sequence before barcode to determine barcode position

def align_barcodes(barcode1, barcode2, max_barcode_diff):
	diff = 0
	if len(barcode1) != len(barcode2):
		return False
	for i, base in enumerate(barcode1):
		if barcode1[i] != barcode2[i]:
			diff += 1
		if diff > max_barcode_diff:
			return False
	return True

def count_barcodes(tf_input_file, cell_input_file, r1_fastq_file, r2_fastq_file, output_prefix): 
	"""
	reads TF barcodes from tf_input_file, reads cell barcodes from cell_input_file, creates a dictionary with cell barcodes and 
	associated TF counts from r1_fastq_file and r2_fastq_file, writes results to files with output_prefix
	tf_input_file: csv file of all possible TF barcodes with 1 barcode per line
	cell_input_file: csv file of all possible cell barcodes with 1 barcode per line
	r1_fastq_file: read 1 fastq file
	r2_fastq_file: read 2 fastq file
	output_prefix: prefix of files to write barcode count results to
	"""

	num_reads = 0 #total number of reads processed
	perfect_matches = 0 # barcodes with perfect match to library
	non_perfect_matches = 0 #number of barcodes without a perfect match to the library
	key_not_found = 0 #count of reads where key was not found
	pair_mismatch = 0

	try:
		# open TF barcodes and initiate dictionary of read counts for each barcode
		with open(tf_input_file, 'r') as infile:   
			reader = csv.reader(infile)
			dictionary = {rows[0]:0 for rows in reader}
		# open cell barcodes and initiate dictionary of read counts for each barcode
		with open(cell_input_file, 'r') as infile:   
			reader = csv.reader(infile)
			cell_dict = {rows[0]:{} for rows in reader}
	except:
		print('could not open input file')

	f_iter = FastqGeneralIterator(open(r1_fastq_file,"r"))
	r_iter = FastqGeneralIterator(open(r2_fastq_file,"r"))
	for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter,r_iter):
		try:
			assert f_id.split()[0] == r_id.split()[0]
		except:
			print(f_id)
			print(r_id)
			pair_mismatch += 1
			break
		num_reads += 1
		cell_barcode = str.upper(str(f_seq))
		if "N" in cell_barcode:
			continue
		read_sequence = str.upper(str(r_seq))
		key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
		key_index = key_region.find(KEY)
		if key_index >= 0:
			start_index = key_index + KEY_REGION_START + len(KEY)
			tf_barcode = read_sequence[start_index:(start_index + BARCODE_LENGTH)]
			barcode_found = False
			for b in dictionary.keys():
				if align_barcodes(b, tf_barcode, MAX_BARCODE_DIFF):
					dictionary[b] += 1
					perfect_matches += 1
					barcode_found = True
					if cell_barcode in cell_dict.keys():
						if b in cell_dict[cell_barcode].keys():
							cell_dict[cell_barcode][b] += 1
						else:
							cell_dict[cell_barcode][b] = 1
			if not barcode_found:
				non_perfect_matches += 1
		else:
			key_not_found += 1

	# create ordered dictionary with barcodes and respective counts and output as a csv file                      
	dict_sorted = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))
	with open(output_prefix + '_barcodes.csv', 'w') as csvfile:
		mywriter = csv.writer(csvfile, delimiter=',')
		for barcode in dict_sorted:
			count = dict_sorted[barcode]
			mywriter.writerow([barcode,count])
	with open(output_prefix + '_cells.csv', 'w') as csvfile:
		mywriter = csv.writer(csvfile, delimiter=',')
		for barcode in cell_dict.keys():
			tfs = cell_dict[barcode]
			tfs_keys = sorted(tfs, key=tfs.__getitem__, reverse=True)
			mywriter.writerow([barcode] + [(x + ':' + str(tfs[x])) for x in tfs_keys])

	# percentage of barcodes that matched perfectly
	if perfect_matches > 0:
		percent_matched = round(perfect_matches/float(perfect_matches + non_perfect_matches) * 100, 1)
	else:
		percent_matched = 0
	# percentage of undetected barcodes with no read counts
	all_barcode_counts = list(dictionary.values())
	barcodes_with_reads = np.count_nonzero(all_barcode_counts)
	barcodes_no_reads = len(all_barcode_counts) - barcodes_with_reads
	percent_no_reads = round(barcodes_no_reads/float(len(all_barcode_counts)) * 100, 1)
	# skew ratio of top 10% to bottom 10% of barcode counts
	top_10 = np.percentile(all_barcode_counts, 90)
	bottom_10 = np.percentile(all_barcode_counts, 10)
	if top_10 != 0 and bottom_10 != 0:
		skew_ratio = top_10/bottom_10
	else:
		skew_ratio = 'Not enough perfect matches to determine skew ratio'

	# write analysis statistics to stats.txt
	with open(output_prefix+'_stats.txt', 'w') as infile:
		infile.write('Number of perfect barcode matches: ' + str(perfect_matches) + '\n')
		infile.write('Number of nonperfect barcode matches: ' + str(non_perfect_matches) + '\n')
		infile.write('Number of reads where key was not found: ' + str(key_not_found) + '\n')
		infile.write('Number of reads processed: ' + str(num_reads) + '\n')
		infile.write('Percentage of barcodes that matched perfectly: ' + str(percent_matched) + '\n')
		infile.write('Percentage of undetected barcodes: ' + str(percent_no_reads) + '\n')
		infile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n')
		infile.write('Pairend mismatches: ' + str(pair_mismatch))
		infile.close()
	return 
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Analyze sequencing data for sgRNA library distribution')
	parser.add_argument('-f', '--fastq1', type=str, dest='r1_fastq_file',
						help='read 1 fastq file name', default='r1_NGS.fastq')
	parser.add_argument('-r', '--fastq2', type=str, dest='r2_fastq_file',
						help='read 2 fastq file name', default='r2_NGS.fastq')
	parser.add_argument('-o', '--output-prefix', type=str, dest='output_prefix',
						help='output file name prefix', default='counts')
	parser.add_argument('-i', '--input', type=str, dest='tf_input_file',
						help='input file name', default='library_sequences.csv')
	parser.add_argument('-ci', '--cell-input', type=str, dest='cell_input_file',
						help='cell input file name', default='cell_barcodes.csv')
	args = parser.parse_args()

	count_barcodes(args.tf_input_file, args.cell_input_file, args.r1_fastq_file, args.r2_fastq_file, args.output_prefix)
