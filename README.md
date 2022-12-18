# Joung_TFAtlas_Manuscript

Scripts from the Joung et al Cell 2023 manuscript on a Transcription Factor Atlas of Directed Differentiation

## Bulk TF screening
- count_barcodes.py counts TF barcodes from fastq
  - input_file: csv file of all possible barcodes with 1 barcode per line
  - fastq_file: read 1 fastq file
  - output_prefix: prefix of files to write barcode count results to

## TF screening with 10x Chromium scRNA-seq
- count_barcodes_10X.py counts TF barcodes from fastq
  - tf_input_file: csv file of all possible TF barcodes with 1 barcode per line
  - cell_input_file: csv file of all possible cell barcodes with 1 barcode per line
	- r1_fastq_file: read 1 fastq file
	- r2_fastq_file: read 2 fastq file
	- output_prefix: prefix of files to write barcode count results to
- map_cells_10X.py maps TF isoforms to cells using the counts_cells.csv file from the count_barcodes_10X.py output
  - ORF_lib_barcode_map.csv: csv file of TF ORF barcode mapping to TF isoforms. Each line contains a TF barcode followed by the corresponding TF isoform
  - counts_cells.csv: csv file of cell barcode mapping to TF ORF barcodes. Each line contains a cell barcode followed by TF barcodes and associated 
  read counts

## TF screening with SHARE-seq
- shareseq_mapTF.sh maps cells to TF ORFs for SHARE-seq using Google Cloud. Calls all other scripts. Based on code from the SHARE-seq manuscript 
DOI: 10.1016/j.cell.2020.09.056.
  - rawdirs: list of NGS run folders
  - writedir: output directory location
  - shareseqdir: directory with share-seq files
  - dir: output directory name
  - Project: project name and prefix for output files
  - Sequencer: Novaseq or Nextseq
  - lanes: number of sequencing lanes
  - bowtieGenome: build a bowtie2 reference for TFs using a fasta file with all TF barcodes
- shareseq_map_cellstoTF.py maps single TF barcodes to cells using SHARE-seq TF mapping bam files. Outputs mapping as a csv file. 
Called by shareseq_mapTF.sh
- shareseq_map_multi_cellstoTF.py maps multiple TF barcodes to cells using SHARE-seq TF mapping bam files. Outputs mapping as a csv file. 
Called by shareseq_mapTF.sh
- shareseq_primerTrim_mapTF.py extracts TF barcodes from SHARE-seq TF mapping fastqs. Called by shareseq_mapTF.sh
- shareseq_updateRGID_mapTF.py updates the RGID of SHARE-seq TF mapping bam files. Called by shareseq_mapTF.sh
