# Python 3.6
'''
# Author: Julia Joung
# Date: 12/13/2022
# Objective: extract TF barcodes from SHARE-seq TF mapping fastqs
'''
import argparse
import itertools
import gzip 
import os

KEY_REGION_START = 20 #start index of key region
KEY_REGION_END = 40 #end index of key region
BARCODE_LENGTH = 24
KEY = "AAGGACGA" #identifies sequence before barcode to determine barcode position

def extractFastq(r1, out, i1):
	R1 = gzip.open(r1, 'r')
	Index1 = gzip.open(i1, 'r')
	outR1 = gzip.open(out, 'w')
	print(("Added index to read name of " + os.path.basename(r1)))
	i = 0
	totalReads = 0
	filteredReads = 0
	valid = False
	for line in zip(R1, Index1):
		if (i == 0):
			Temp1 = line[0][0:-1].decode("utf-8")
			Temp1 = str.encode(Temp1.replace(" ", "_"))
			totalReads += 1
		elif (i == 1):
			read_sequence = str.upper(line[0].decode("utf-8"))
			key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
			key_index = key_region.find(KEY)
			if key_index >= 0:
				start_index = key_index + KEY_REGION_START + len(KEY)
				Temp2 = str.encode(read_sequence[start_index:(start_index + BARCODE_LENGTH)])
				id1 = line[1]
				valid = True
		elif (i == 2):
			Temp3 = line[0]
		elif (i == 3):
			if valid:
				Temp4 = line[0][start_index:(start_index + BARCODE_LENGTH)]
		i += 1
		if (i == 4):
			if valid and len(Temp2) == len(Temp4):
				outR1.write(Temp1 + id1)
				outR1.write(Temp2 + str.encode("\n"))
				outR1.write(Temp3)
				outR1.write(Temp4 + str.encode("\n"))
				filteredReads += 1
				valid = False
				id1 = str.encode("")
			i = 0
	print(str(totalReads))
	print(str(filteredReads))
	R1.close()
	outR1.close()

def main():
	parser = argparse.ArgumentParser(
		description="Splits fastqs based on splitATAC project",
		epilog="Intended for splitATAC data")
	parser.add_argument(
		"-R1",
		metavar="Read 1",
		required=True,
		help="Path to the Read 1 fastq file")
	parser.add_argument(
		"-Index1",
		metavar="Output",
		required=True,
		help="Path to the index1 files")
	parser.add_argument(
		"--out",
		metavar="Output",
		required=True,
		help="Path to the output fastq files")
	args = parser.parse_args()
	extractFastq(args.R1, args.out, args.Index1)


main()
