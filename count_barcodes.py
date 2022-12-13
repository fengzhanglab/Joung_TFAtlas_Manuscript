#counts TF barcodes from fastq for bulk TF screening

from Bio import SeqIO
import csv
from collections import OrderedDict
import numpy as np
import sys
import argparse

KEY_REGION_START = 25 #start index of key region
KEY_REGION_END = 50 #end index of key region
BARCODE_LENGTH = 24
KEY = "GAAAGGACGA" #identifies sequence before barcode to determine barcode position

def count_barcodes(input_file, fastq_file, output_file): 
	"""
	reads barcodes from input_file, creates a dictionary with barcode counts from fastq_file, writes to output_file
	dictionary: csv file of all possible barcodes with 1 barcode per line
	fastq_file: forward read fastq file
	output_file: csv file to write barcode counts to
	"""

	num_reads = 0 #total number of reads processed
	perfect_matches = 0 # barcodes with perfect match to library
	non_perfect_matches = 0 #number of barcodes without a perfect match to the library
	key_not_found = 0 #count of reads where key was not found

	# open TF barcodes and initiate dictionary of read counts for each barcode
	try:
		with open(input_file, 'r') as infile:  
			reader = csv.reader(infile)
			dictionary = {rows[0]:0 for rows in reader}
	except:
		print('could not open ' + input_file)
	  
	# open fastq file
	try:
		handle = open(fastq_file, "r")
	except:
		print("could not find fastq file")
		return

	# process reads in fastq file
	readiter = SeqIO.parse(handle, "fastq")
	for record in readiter: #contains the seq and Qscore etc.
		num_reads += 1
		read_sequence = str.upper(str(record.seq))
		key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
		key_index = key_region.find(KEY)
		if key_index >= 0:
			start_index = key_index + KEY_REGION_START + len(KEY)
			barcode = read_sequence[start_index:(start_index + BARCODE_LENGTH)]
			if barcode in dictionary:
				dictionary[barcode] += 1
				perfect_matches += 1
			else:
				non_perfect_matches += 1
		else:
			key_not_found += 1

	# create ordered dictionary with barcodes and respective counts and output as a csv file                      
	dict_sorted = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))
	with open(output_file + '_counts.csv', 'w') as csvfile:
		mywriter = csv.writer(csvfile, delimiter=',')
		for barcode in dict_sorted:
			count = dict_sorted[barcode]
			mywriter.writerow([barcode,count])

	# percentage of barcodes that matched perfectly
	percent_matched = round(perfect_matches/float(perfect_matches + non_perfect_matches) * 100, 1)
	# percentage of undetected barcodes with no read counts
	barcode_counts = list(dictionary.values())
	barcodes_with_reads = np.count_nonzero(barcode_counts)
	barcodes_no_reads = len(barcode_counts) - barcodes_with_reads
	percent_no_reads = round(barcodes_no_reads/float(len(barcode_counts)) * 100, 1)
	# skew ratio of top 10% to bottom 10% of barcode counts
	top_10 = np.percentile(barcode_counts, 90)
	bottom_10 = np.percentile(barcode_counts, 10)
	if top_10 != 0 and bottom_10 != 0:
		skew_ratio = top_10/bottom_10
	else:
		skew_ratio = 'Not enough perfect matches to determine skew ratio'

	# write analysis statistics to statistics.txt
	with open(output_file + '_stats.txt', 'w') as infile:
		infile.write('Number of perfect barcode matches: ' + str(perfect_matches) + '\n')
		infile.write('Number of nonperfect barcode matches: ' + str(non_perfect_matches) + '\n')
		infile.write('Number of reads where key was not found: ' + str(key_not_found) + '\n')
		infile.write('Number of reads processed: ' + str(num_reads) + '\n')
		infile.write('Percentage of barcodes that matched perfectly: ' + str(percent_matched) + '\n')
		infile.write('Percentage of undetected barcodes: ' + str(percent_no_reads) + '\n')
		infile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio))
		infile.close()

	handle.close()           
	return 
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Analyze sequencing data for sgRNA library distribution')
	parser.add_argument('-f', '--fastq', type=str, dest='fastq_file',
						help='fastq file name', default='NGS.fastq')
	parser.add_argument('-o', '--output', type=str, dest='output_file',
						help='output file prefix', default='library')
	parser.add_argument('-i', '--input', type=str, dest='input_file',
						help='input file name', default='library_sequences.csv')
	args = parser.parse_args()

	count_barcodes(args.input_file, args.fastq_file, args.output_file)
