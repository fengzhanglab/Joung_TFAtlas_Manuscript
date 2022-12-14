# Python 3.6
'''
# Author: Julia Joung
# Date: 12/13/2022
# Objective: Map TF barcodes to cells using SHARE-seq TF mapping bam files. Outputs mapping as a csv file
'''

import pysam
import argparse
import subprocess
import csv

BARCODE_LENGTH = 24
READ_MIN = 5
MULTIPLE_MIN = 2
def mapTFs(bamfile, outfile):
	bam = pysam.AlignmentFile(bamfile, 'rb')
	num_reads = 0
	cellToTF = dict()
	for read in bam.fetch():
		num_reads += 1
		barcodeMatch = 0
		name = read.query_name
		nameIndex = name.find('R1.')
		cellBarcode = name[nameIndex:nameIndex+17]
		TFBarcode = read.get_reference_sequence().upper()
		if len(TFBarcode) != BARCODE_LENGTH:
			continue
		if cellBarcode in cellToTF.keys():
			TFdict = cellToTF[cellBarcode]
			if TFBarcode in TFdict.keys():
				TFdict[TFBarcode] += 1
			else:
				TFdict[TFBarcode] = 1
		else:
			TFdict = {TFBarcode:1}
			cellToTF[cellBarcode] = TFdict

	num_unmapped = 0
	num_mapped = 0
	num_amb = 0
	total = 0
	with open(outfile, 'w') as csvfile:
		csvwriter = csv.writer(csvfile)
		for cell, TFs in cellToTF.items():
			total += 1
			TFs_sorted = sorted(TFs.items(), key=lambda k:k[1], reverse=True)
			if TFs_sorted[0][1] >= READ_MIN:
				if len(TFs_sorted) > 1:
					if TFs_sorted[0][1] >= MULTIPLE_MIN*TFs_sorted[1][1]:
						csvwriter.writerow([cell, TFs_sorted[0][0], str(TFs_sorted[0][1])])
						num_mapped += 1
					else:
						num_amb += 1
				else:
					csvwriter.writerow([cell, TFs_sorted[0][0], str(TFs_sorted[0][1])])
					num_mapped += 1
			else:
				num_unmapped += 1

	print("Number of cells unmapped: " + str(num_unmapped))
	print("Number of cells mapped: " + str(num_mapped))
	print("Number of cells ambiguous: " + str(num_amb))
	print("Number of cells total: " + str(total))
	print("Number of reads total: " + str(num_reads))
	bam.close()

def main():
	parser = argparse.ArgumentParser(
		description="Maps TF barcodes to cells")
	parser.add_argument(
		"--bam",
		metavar="Input",
		required=True,
		help="Path to bam files annotated with RGID")
	parser.add_argument(
		"--out",
		metavar="Output",
		required=True,
		help="Path to csv file of cells mapped with TFs")
	args = parser.parse_args()
	pysam.index(args.bam)
	mapTFs(args.bam, args.out)

main()
