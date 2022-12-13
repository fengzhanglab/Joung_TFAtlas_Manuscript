#maps TF isoforms to cells using the counts_cells.csv file from the count_barcodes_10X.py output

import csv

#open csv file of TF ORF barcode mapping to TF isoforms. Each line contains a TF barcode followed by the corresponding TF isoform
with open("ORF_lib_barcode_map.csv", 'r') as infile:
	barcode_map = {row[0]:row[1] for row in csv.reader(infile.read().splitlines())}

#open csv file of cell barcode mapping to TF ORF barcodes. Each line contains a cell barcode followed by TF barcodes and associated read counts
with open("counts_cells.csv", 'r') as infile:
	data = [row for row in csv.reader(infile.read().splitlines())]

#map cells to TFs. NA indicates unmapped cells and AMB indicates ambiguous mapping
with open("Cell_ORFs.csv", 'w') as csvfile:
	csvwriter = csv.writer(csvfile)
	for row in data:
		if len(row) == 1:
			csvwriter.writerow(row + ['NA'])
		elif len(row) == 2:
			count = int(row[1].split(':')[1])
			if count > 1:
				orf = barcode_map[row[1].split(':')[0]]
				csvwriter.writerow([row[0], orf])
			else:
				csvwriter.writerow([row[0],'AMB'])
		else:
			top = int(row[1].split(':')[1])
			second = int(row[2].split(':')[1])
			if top < second*1.25:
				csvwriter.writerow([row[0],'AMB'])
			else:
				orf = barcode_map[row[1].split(':')[0]]
				csvwriter.writerow([row[0], orf])
		