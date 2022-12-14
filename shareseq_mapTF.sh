#!/bin/bash
# Author: Julia Joung <julia@joung.science>
# Last modified date: 2022/11/28
# Designed for mapping cells with TF ORFs for SHARE-seq using Google Cloud. Based on code from the SHARE-seq manuscript DOI: 10.1016/j.cell.2020.09.056

rawdirs=(201218_SL-NVR_0049_AHMCTNDSXY) #list of NGS run folders
writedir=/mnt/disks/10tb #output directory location
shareseqdir=shareseq #directory with share-seq files
dir=$writedir/201218_SHAREseq_200k_dialout #output directory name
Project=(TFAtlas.200k) #project name and prefix for output files
Type=(RNA)
Species=TFs
Start=Bcl #specifies where to start the pipeline. Bcl or Fastq_Merge or Fastq_SplitLane
Runtype=full #QC or full, QC only analyze 12M reads
Sequencer=Novaseq #Novaseq or Nextseq	
cores=40 #number of computing cores
lanes=4 #number of sequencing lanes

myPATH=~/shareseq/
picardPATH=/opt/picard-tools/picard.jar
bowtieGenome=~/shareseq/refGenome/bowtie2/ #build a bowtie2 reference for TFs using a fasta file with all TF barcodes

echo "the number of projects is" ${#Project[@]}
echo "Running $Runtype pipeline"
echo "the name of the species is $Species"

for rawdir in ${rawdirs[@]}; do
	if [ ! -d $rawdir ]; then 
		gsutil cp gs://shareseq_bucket/$rawdir.tgz .
		tar xf $rawdir.tgz
		rm $rawdir.tgz
	fi
done

cd ~/
if [ ! -d ~/$shareseqdir ]; then 
	gsutil cp gs://shareseq_bucket/$shareseqdir.tgz 
	tar xf $shareseqdir.tgz
	rm $shareseqdir.tgz
fi

shareseqdir=~/$shareseqdir
if [ ! -d $dir ]; then mkdir $dir; fi
if [ ! -d $dir/fastqs ]; then mkdir $dir/fastqs ; fi
if [ ! -d $dir/temp ]; then mkdir $dir/temp ; fi
cp $writedir/shareseq_mapTF.sh $dir/
cd $dir 
if [ -f $dir/Run.log ]; then rm $dir/Run.log; fi

if (( $lanes > 1 )); then
	for i in $( seq 1 $lanes ); do
		Start=Bcl
		Name=${Project[0]}.L${i}
		echo "the name of the project is $Name"
		# if start with bcl file
		if [ "$Start" = Bcl ]; then
			echo "Bcl2fastq"
			for rawdir in ${rawdirs[@]}; do
				rawdir=$writedir/$rawdir
				if [ -f $rawdir/fastqs/Undetermined_S1_L1_R2_001.fastq.gz ]; then
					echo "Found Undetermined_S0_L001_I1_001.fastq.gz, skip Bcl2Fastq"
				else
					echo "Converting bcl to fastq"
					mkdir $rawdir/fastqs/
					bcl2fastq -p $cores -R $rawdir --mask-short-adapter-reads 0 -o $rawdir/fastqs/ --create-fastq-for-index-reads 2>>$dir/Run.log
				fi
				cd $rawdir/fastqs/
				mv Undetermined_S0_L00${i}_R1_001.fastq.gz Undetermined_S1_L${i}_R1_001.fastq.gz
				mv Undetermined_S0_L00${i}_I1_001.fastq.gz Undetermined_S1_L${i}_R2_001.fastq.gz
				rm filler*
			done
			rawdir=$writedir/${rawdirs[0]}/fastqs
		fi

		Start=Fastq_Merge
		if [ "$Start" = Fastq_SplitLane ] || [ "$Start" = Fastq_Merge ]; then
			echo "Skip bcltofastq"
			if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
			cd $rawdir
			if [ "$Start" = Fastq_SplitLane ]; then
			temp=$(ls 1_*.1.1.fastq.gz)
			Run=$(echo $temp | sed -e 's/^1_//' | sed -e 's/.1.1.fastq.gz//')
			echo "Run number is:" $Run
			echo "Spliting to 3 million reads files"
			mkdir $dir/smallfastqs/
			dosplit(){
				fastp -i $3/$1_$4.$1.barcode_1.fastq.gz -o $2/smallfastqs/$4.$1.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
				fastp -i $3/$1_$4.$1.barcode_2.fastq.gz -o $2/smallfastqs/$4.$1.I2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
				fastp -i $3/$1_$4.$1.1.fastq.gz -o $2/smallfastqs/$4.$1.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q
				fastp -i $3/$1_$4.$1.2.fastq.gz -o $2/smallfastqs/$4.$1.R2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
			}
			export -f dosplit
			dosplitQC(){
				zcat $3/$1_$4.$1.1.fastq.gz | head -n 12100000 | gzip > $2/temp1.$1.fastq.gz &
				zcat $3/$1_$4.$1.2.fastq.gz | head -n 12100000 | gzip > $2/temp2.$1.fastq.gz &
				zcat $3/$1_$4.$1.barcode_1.fastq.gz | head -n 12100000 | gzip > $2/temp3.$1.fastq.gz &
				zcat $3/$1_$4.$1.barcode_2.fastq.gz | head -n 12100000 | gzip > $2/temp4.$1.fastq.gz &
				wait
				
				fastp -i $2/temp1.$1.fastq.gz -o $2/smallfastqs/$4.$1.R1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
				fastp -i $2/temp2.$1.fastq.gz -o $2/smallfastqs/$4.$1.R2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q
				fastp -i $2/temp3.$1.fastq.gz -o $2/smallfastqs/$4.$1.I1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
				fastp -i $2/temp4.$1.fastq.gz -o $2/smallfastqs/$4.$1.I2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
			}
			export -f dosplitQC
			if [ -f $dir/smallfastqs/0001.$Run.1.I1.fastq.gz ]; then
					echo "Found 0001.Undetermined.1.I1.fastq, skip split fastqs"
			else
				if [ "$Runtype" = QC ]; then
					if [ "$Sequencer" = Novaseq ]; then
						parallel --delay 1 --jobs 2 dosplitQC {} $dir $rawdir $Run ::: {1..2}
					else
						parallel --delay 1 --jobs 4 dosplitQC {} $dir $rawdir $Run ::: {1..4}
					fi
			#		rm $dir/temp*.fastq.gz
				elif [ "$Runtype" = full ]; then
					if [ "$Sequencer" = Novaseq ]; then
						parallel --delay 1 --jobs 4 dosplit {} $dir $rawdir $Run ::: {1..4}
					else
						parallel --delay 1 --jobs 2 dosplit {} $dir $rawdir $Run ::: {1..2}
					fi
				else
					echo "Unknown sequencer type, exiting" && exit;
				fi
			fi
			elif [ "$Start" = Fastq_Merge ]; then
				temp=$(ls *S1_L1_R1_001.fastq.gz)
				Run=$(echo $temp | sed -e 's/\_\S1\_\L1\_\R1\_\001.fastq.gz//')
				echo "Run number is:" $Run
				echo "Spliting to 3 million reads files"
				mkdir $dir/smallfastqs/
				if [ -f $dir/smallfastqs/0001.$Run.L${i}.R1.fastq.gz ]; then
					echo "Found 0001.Undetermined.1.R1.fastq, skip split fastqs"
				else
					if [ "$Runtype" = QC ]; then
						echo "Runing QC pipeline"
						zcat $rawdir/"$Run"_S1_R1_001.fastq.gz | head -n 12100000 | sed 's/1:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp1.fastq.gz				
						fastp -i $dir/temp1.fastq.gz -o $dir/smallfastqs/$Run.1.R1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log
						rm $dir/temp*.fastq.gz
					elif [ "$Runtype" = full ]; then
						echo "Runing full pipeline"
						if [ ! -f $dir/temp3.L${i}.fastq.gz ]; then
							echo "Modify fastqs"
							zcat $rawdir/"$Run"_S1_L${i}_R1_001.fastq.gz | sed 's/1:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp1.L${i}.fastq.gz &
							zcat $rawdir/"$Run"_S1_L${i}_R2_001.fastq.gz | sed 's/2:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp3.L${i}.fastq.gz &
							wait
						fi
						echo "Split fastqs to small files"
						fastp -i $dir/temp1.L${i}.fastq.gz -o $dir/smallfastqs/$Run.L${i}.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
						fastp -i $dir/temp3.L${i}.fastq.gz -o $dir/smallfastqs/$Run.L${i}.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
						wait
						rm $dir/temp*.fastq.gz
					else
						echo "Unknown sequencer type, exiting" && exit;
					fi
				fi
			fi
			
			if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
			ls $dir/smallfastqs | grep L${i}.R1 > $dir/filesr1.xls
			ls $dir/smallfastqs | grep L${i}.I1 > $dir/filesi1.xls
			
			cd $dir/
			mkdir $dir/Indexed/
			if [ -f $dir/Indexed/Sub.0001.$Run.L${i}.R1.fastq.gz ]; then
				echo "Found Sub.0001.Undetermined.1.R1.fastq.gz, skip adding index"
			else
				echo "Adding index to fastqs"
				paste filesr1.xls filesi1.xls | awk -v OFS='\t' '{print $1, $2}'> Filelist2.xls
				parallel --jobs 4 --colsep '\t' ' if [ -f '$dir'/Indexed/Sub.{1} ]; then echo "found '$dir'/Indexed/Sub.{1}"; \
						  else  python3 '$myPATH'/shareseq_primerTrim_mapTF.py -R1 '$dir'/smallfastqs/{1} \
						  --out '$dir'/Indexed/Sub.{1} -Index1 '$dir'/smallfastqs/{2}; fi 2>>'$dir'/Run.log' :::: Filelist2.xls
				cd $dir/Indexed
				cat Sub*L${i}.R1.fastq.gz > $Name.L${i}.R1.trim.fastq.gz
				gunzip $Name.L${i}.R1.trim.fastq.gz
				mv $Name.L${i}.R1.trim.fastq $dir/fastqs/
				rawdir=$dir/Indexed/
			fi
		fi 

		if [ ! -d $dir/$Name ]; then mkdir $dir/$Name; fi
		cd $dir/$Name
		if [ -f $Name.$Species.L${i}.bam ]; then
			echo "Found $Name.$Species.bam, skip alignment"
		else
			(bowtie2 -p $cores --rg-id $Name -x $bowtieGenome/$Species/$Species -U $dir/fastqs/$Name.L${i}.R1.trim.fastq -N 1 -L 24 | \
			 samtools view -bS - -o $Name.$Species.L${i}.bam) 2>$Name.$Species.L${i}.align.log
		fi
		if [ -f $Name.$Species.L${i}.st.bam ]; then
			echo "Found $Name.$Species.st.bam, skip sorting"
		else
			java -Xmx24g -Djava.io.tmpdir=$dir/temp/ -jar $picardPATH SortSam SO=coordinate I=$Name.$Species.L${i}.bam O=$Name.$Species.L${i}.st.bam VALIDATION_STRINGENCY=LENIENT TMP_DIR=$dir/temp/ 2>>$dir/Run.log
			samtools index -@ $cores $Name.$Species.L${i}.st.bam
		fi	
		# Update RGID's and Headers
		if [ -f $Name.$Species.rigid.reheader.L${i}.st.bam.bai ]; then
			echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RGID"
		else
			echo "Update RGID for $Name.$Species.st.bam"
			if [ -f $Name.$Species.rigid.L${i}.st.bam ]; then
			echo "Found $Name.$Species.rigid.st.bam, skip update RG tag"
			else
				python3 $myPATH/shareseq_updateRGID_mapTF.py --bam $Name.$Species.L${i}.st.bam --out $Name.$Species.rigid.L${i}.st.bam --err $Name.$Species.discard.L${i}.st.bam --libtype ${Type[$index]}
			fi
			samtools view -H $Name.$Species.L${i}.st.bam > $Name.$Species.L${i}.st.header.sam
			samtools view -@ $cores $Name.$Species.rigid.L${i}.st.bam | cut -f1 | sed 's/_/\t/g' | cut -f2 | sort | uniq | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.L${i}.sam
			sed -e '/\@RG/r./header.temp.L${i}.sam' $Name.$Species.L${i}.st.header.sam > $Name.$Species.rigid.L${i}.st.header.sam
			samtools reheader -P $Name.$Species.rigid.L${i}.st.header.sam $Name.$Species.rigid.L${i}.st.bam > $Name.$Species.rigid.reheader.L${i}.st.bam
			samtools index -@ $cores $Name.$Species.rigid.reheader.L${i}.st.bam
			rm *rigid.L${i}.st*
			rm $Name.$Species.L${i}.st.header.sam header.temp.L${i}.sam
		fi
		#map cells to TFs
		if [ -f $Name.$Species.TFmap.L${i}.csv ]; then
			echo "Found $Name.$Species.TFmap.L${i}.csv, skip TFmap"
		else
			echo "Map cells to TFs"
			python3 $myPATH/shareseq_map_cellstoTF.py --bam $Name.$Species.rigid.reheader.L${i}.st.bam --out $Name.$Species.TFmap.L${i}.csv
		fi
	done
else
    Start=Bcl
    Name=${Project[0]}
    echo "the name of the project is $Name"
	if [ "$Start" = Bcl ]; then
	    echo "Bcl2fastq"
	    rawdir=$writedir/$rawdir
		if [ -f $rawdir/fastqs/Undetermined_S0_R1_001.fastq.gz ] || [ -f $rawdir/fastqs/Undetermined_S1_R1_001.fastq.gz ]; then
		echo "Found Undetermined_S0_L001_I1_001.fastq.gz, skip Bcl2Fastq"
		else
		echo "Converting bcl to fastq"
		mkdir $rawdir/fastqs/
		bcl2fastq -p $cores -R $rawdir --no-lane-splitting --mask-short-adapter-reads 0 -o $rawdir/fastqs/ --create-fastq-for-index-reads 2>>$dir/Run.log
		fi
		cd $rawdir/fastqs/	
		mv Undetermined_S0_R1_001.fastq.gz Undetermined_S1_R1_001.fastq.gz
		mv Undetermined_S0_I1_001.fastq.gz Undetermined_S1_R2_001.fastq.gz
		rm filler*
		if (( ${#rawdirs[@]} > 1 )); then
			if [ ! -f $writedir/${rawdirs[0]}/original_fastqs/Original_S1_R2_001.fastq.gz ]; then
				mkdir $writedir/${rawdirs[0]}/original_fastqs/
				for i in 1 2 3 4; do
					cat $writedir/${rawdirs[0]}/fastqs/Undetermined_S1_R${i}_001.fastq.gz $writedir/${rawdirs[1]}/fastqs/Undetermined_S1_R${i}_001.fastq.gz > ${writedir}/${rawdirs[0]}/fastqs/temp${i}.fastq.gz
					mv $writedir/${rawdirs[0]}/fastqs/Undetermined_S1_R${i}_001.fastq.gz $writedir/${rawdirs[0]}/original_fastqs/Original_S1_R${i}_001.fastq.gz
					mv $writedir/${rawdirs[0]}/fastqs/temp${i}.fastq.gz $writedir/${rawdirs[0]}/fastqs/Undetermined_S1_R${i}_001.fastq.gz
				done
			fi
		fi
		rawdir=$writedir/${rawdirs[0]}/fastqs
	fi

	Start=Fastq_Merge
	if [ "$Start" = Fastq_SplitLane ] || [ "$Start" = Fastq_Merge ]; then
		echo "Skip bcltofastq"
		if [ -f $dir/fastqs/filesi1.xls ]; then rm $dir/fastqs/filler*; fi
		cd $rawdir
		if [ "$Start" = Fastq_SplitLane ]; then
		temp=$(ls 1_*.1.1.fastq.gz)
		Run=$(echo $temp | sed -e 's/^1_//' | sed -e 's/.1.1.fastq.gz//')
		echo "Run number is:" $Run
		echo "Spliting to 3 million reads files"
		mkdir $dir/smallfastqs/
		dosplit(){
			fastp -i $3/$1_$4.$1.barcode_1.fastq.gz -o $2/smallfastqs/$4.$1.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
			fastp -i $3/$1_$4.$1.barcode_2.fastq.gz -o $2/smallfastqs/$4.$1.I2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
			fastp -i $3/$1_$4.$1.1.fastq.gz -o $2/smallfastqs/$4.$1.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q
			fastp -i $3/$1_$4.$1.2.fastq.gz -o $2/smallfastqs/$4.$1.R2.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 
		}
		export -f dosplit
		dosplitQC(){
			zcat $3/$1_$4.$1.1.fastq.gz | head -n 12100000 | gzip > $2/temp1.$1.fastq.gz &
			zcat $3/$1_$4.$1.2.fastq.gz | head -n 12100000 | gzip > $2/temp2.$1.fastq.gz &
			zcat $3/$1_$4.$1.barcode_1.fastq.gz | head -n 12100000 | gzip > $2/temp3.$1.fastq.gz &
			zcat $3/$1_$4.$1.barcode_2.fastq.gz | head -n 12100000 | gzip > $2/temp4.$1.fastq.gz &
			wait
			
			fastp -i $2/temp1.$1.fastq.gz -o $2/smallfastqs/$4.$1.R1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
			fastp -i $2/temp2.$1.fastq.gz -o $2/smallfastqs/$4.$1.R2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q
			fastp -i $2/temp3.$1.fastq.gz -o $2/smallfastqs/$4.$1.I1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
			fastp -i $2/temp4.$1.fastq.gz -o $2/smallfastqs/$4.$1.I2.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 
		}
		export -f dosplitQC
		if [ -f $dir/smallfastqs/0001.$Run.1.I1.fastq.gz ]; then
				echo "Found 0001.Undetermined.1.I1.fastq, skip split fastqs"
		else
			if [ "$Runtype" = QC ]; then
				if [ "$Sequencer" = Novaseq ]; then
					parallel --delay 1 --jobs 2 dosplitQC {} $dir $rawdir $Run ::: {1..2}
				else
					parallel --delay 1 --jobs 4 dosplitQC {} $dir $rawdir $Run ::: {1..4}
				fi
		#		rm $dir/temp*.fastq.gz
			elif [ "$Runtype" = full ]; then
				if [ "$Sequencer" = Novaseq ]; then
					parallel --delay 1 --jobs 4 dosplit {} $dir $rawdir $Run ::: {1..4}
				else
					parallel --delay 1 --jobs 2 dosplit {} $dir $rawdir $Run ::: {1..2}
				fi
			else
				echo "Unknown sequencer type, exiting" && exit;
			fi
		fi
		elif [ "$Start" = Fastq_Merge ]; then
			temp=$(ls *S1_R1_001.fastq.gz)
			Run=$(echo $temp | sed -e 's/\_\S1\_\R1\_\001.fastq.gz//')
			echo "Run number is:" $Run
			echo "Spliting to 3 million reads files"
			mkdir $dir/smallfastqs/
			if [ -f $dir/smallfastqs/0001.$Run.1.R1.fastq.gz ]; then
				echo "Found 0001.Undetermined.1.R1.fastq, skip split fastqs"
			else
				if [ "$Runtype" = QC ]; then
					echo "Runing QC pipeline"
					zcat $rawdir/"$Run"_S1_R1_001.fastq.gz | head -n 12100000 | sed 's/1:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp1.fastq.gz				
					fastp -i $dir/temp1.fastq.gz -o $dir/smallfastqs/$Run.1.R1.fastq.gz -S 1000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log
					rm $dir/temp*.fastq.gz
				elif [ "$Runtype" = full ]; then
					echo "Runing full pipeline"
					if [ ! -f $dir/temp2.fastq.gz ]; then
						echo "Modify fastqs"
						zcat $rawdir/"$Run"_S1_R1_001.fastq.gz | sed 's/1:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp1.fastq.gz &
						zcat $rawdir/"$Run"_S1_R2_001.fastq.gz | sed 's/2:N:0:1/1:N:0:/g' | sed 's/1:N:0:0/1:N:0:/g' | sed 's/^$/N/' | gzip > $dir/temp3.fastq.gz &
						wait
					fi
					echo "Split fastqs to small files"
					fastp -i $dir/temp1.fastq.gz -o $dir/smallfastqs/$Run.1.R1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
					fastp -i $dir/temp3.fastq.gz -o $dir/smallfastqs/$Run.1.I1.fastq.gz -S 12000000 --thread 1 -d 4 -A -G -L -Q 2>>$dir/split.log &
					wait
					rm $dir/temp*.fastq.gz
				else
					echo "Unknown sequencer type, exiting" && exit;
				fi
			fi
		fi
		
		if [ -f $dir/fastp.json ]; then rm $dir/fastp.json $dir/fastp.html; fi
		ls $dir/smallfastqs | grep R1 > $dir/filesr1.xls
		ls $dir/smallfastqs | grep I1 > $dir/filesi1.xls
		
		cd $dir/
		mkdir $dir/Indexed/
		if [ -f $dir/Indexed/Sub.0001.$Run.1.R1.fastq.gz ]; then
			echo "Found Sub.0001.Undetermined.1.R1.fastq.gz, skip adding index"
		else
			echo "Adding index to fastqs"
			paste filesr1.xls filesi1.xls | awk -v OFS='\t' '{print $1, $2}'> Filelist2.xls
			parallel --jobs 4 --colsep '\t' ' if [ -f '$dir'/Indexed/Sub.{1} ]; then echo "found '$dir'/Indexed/Sub.{1}"; \
					  else  python3 '$myPATH'/shareseq_primerTrim_mapTF.py -R1 '$dir'/smallfastqs/{1} \
					  --out '$dir'/Indexed/Sub.{1} -Index1 '$dir'/smallfastqs/{2}; fi 2>>'$dir'/Run.log' :::: Filelist2.xls
			cd $dir/Indexed
			cat Sub*R1.fastq.gz > $Name.R1.trim.fastq.gz
			gunzip $Name.R1.trim.fastq.gz
			mv $Name.R1.trim.fastq $dir/fastqs/
			rawdir=$dir/Indexed/
		fi
	fi 

	if [ ! -d $dir/$Name ]; then mkdir $dir/$Name; fi
	cd $dir/$Name
	if [ -f $Name.$Species.bam ]; then
		echo "Found $Name.$Species.bam, skip alignment"
	else
		(bowtie2 -p $cores --rg-id $Name \
			 -x $bowtieGenome/$Species/$Species \
			 -U $dir/fastqs/$Name.R1.trim.fastq | \
			 samtools view -bS - -o $Name.$Species.bam) 2>$Name.$Species.align.log
	fi
	if [ -f $Name.$Species.st.bam ]; then
		echo "Found $Name.$Species.st.bam, skip sorting"
	else
		java -Xmx24g -Djava.io.tmpdir=$dir/temp/ -jar $picardPATH SortSam SO=coordinate I=$Name.$Species.bam O=$Name.$Species.st.bam VALIDATION_STRINGENCY=LENIENT TMP_DIR=$dir/temp/ 2>>$dir/Run.log
		samtools index -@ $cores $Name.$Species.st.bam
	fi	
	# Update RGID's and Headers
	if [ -f $dir/$Name/$Name.$Species.rigid.reheader.st.bam.bai ]; then
		echo "Found $Name.$Species.rigid.reheader.st.bam, skip update RGID"
	else
		echo "Update RGID for $Name.$Species.st.bam"
		if [ -f $Name.$Species.rigid.st.bam ]; then
		echo "Found $Name.$Species.rigid.st.bam, skip update RG tag"
		else
			python3 $myPATH/shareseq_updateRGID_mapTF.py --bam $Name.$Species.st.bam --out $Name.$Species.rigid.st.bam --err $Name.$Species.discard.st.bam --libtype ${Type[$index]}
		fi
		samtools view -H $Name.$Species.st.bam > $Name.$Species.st.header.sam
		samtools view -@ $cores $Name.$Species.rigid.st.bam | cut -f1 | sed 's/_/\t/g' | cut -f2 | sort | uniq | awk -v OFS='\t' '{print "@RG", "ID:"$1, "SM:Barcode"NR}' > header.temp.sam
		sed -e '/\@RG/r./header.temp.sam' $Name.$Species.st.header.sam > $Name.$Species.rigid.st.header.sam
		samtools reheader -P $Name.$Species.rigid.st.header.sam $Name.$Species.rigid.st.bam > $Name.$Species.rigid.reheader.st.bam
		samtools index -@ $cores $Name.$Species.rigid.reheader.st.bam
		rm *rigid.st*
		rm $Name.$Species.st.header.sam header.temp.sam
	fi
	#map cells to TFs
	echo "Map cells to TFs"
	python3 $myPATH/shareseq_map_cellstoTF.py --bam $Name.$Species.rigid.reheader.st.bam --out $Name.$Species.TFmap.csv
fi

exit;
