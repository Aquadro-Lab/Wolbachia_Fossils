#!/bin/bash

#Use as bash readfile SRAfilename.txt Wolbachia_accession
#where SRAfilename.txt is a file with the SRA run IDs with an extra "end" line after all the IDs

timestamp() {
	date +"%T"
}

timestamp


echo "$1"

while IFS= read -r line; do

    echo "------------ SRA FILE: $line ------------"
    timestamp

#comment out this section if you already have the fastq reads because it takes a long time to download them.

echo "------------ Downloading reads with fastq-dump ------------" 

fasterq-dump --split-3 --progress $line
mkdir $line.dir 
mv *.fastq $line.dir

pwd

#resume pipeline here if you commented out the above section. this section needs prior "bwa index" of the reference fasta in order to work.

    echo "------------ Aligning reads to reference with BWA mem------------" 

    i3=$line.dir
    i2="${i3%.dir}"
    i4="${i2}_1"
    i5="${i2}_2"
    i6="${i2}.bam"

    echo $i3
	echo $i2
	echo $i4
	echo $i5
	echo $i6

	pwd

	bwa mem -t 24 GCF_004382195.2_ASM438219v2_genomic-Wol_pan_genome95.fasta $i3/$i4* $i3/$i5* | samtools sort -o $i6 #have to change Dsim_wNo.fasta if using a different reference sequence to map the reads to.


	pwd #added this 4.10.23 b/c got an error that said can't find the bam file. maybe we are in the wrong directory? not sure.

	echo "------------ Sorting and Indexing with Samtools ------------"

#	cd $line.dir

	i7="${i2}.sort.bam"
	echo $i7

	samtools sort $i6 -o $i7
	
	samtools index $i7


	echo "------------ Subset Wolbachia hits from sorted bam file ------------"

	i14="${i2}.sort2.bam"
	echo $i14

	samtools view $i7 | grep 'GCF' > $i14


echo "------------ Finding BND structural variants with Delly ------------"

	export PATH=/programs/delly-0.8.7:$PATH

	i8="${i2}.bcf"

#have to change Dsim_wNo.fasta if using a different reference sequence to map the reads to 
	
	/programs/delly-0.8.7/delly call -t BND -o $i8 -g GCF_004382195.2_ASM438219v2_genomic-Wol_pan_genome95.fasta $i7 
#	delly call -t BND -o $i8 -g /../../workdir/mjw349/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic-Wol_pan_genome95.fasta $i7 
	i9="${i2}.vcf"

	bcftools view $i8 > $i9

	echo "------------ Extract lines with Wolbachia hits ------------"

	i13="${i2}_GCF_004382195.2_ASM438219v2_genomic-Wol_pan_genome95_HITS.txt"


	grep 'GCF' $i9 > $i13


	echo "------------ Moving everything to the $line folder ------------"

	mv $i6 $line.dir
	mv $i7 $line.dir
	mv $i8 $line.dir
	mv $i9 $line.dir
##	mv $i12 $line.dir

	i10="${i8}.csi"
	i11="${i7}.bai"

	mv $i10 $line.dir
	mv $i11 $line.dir
	mv $i13 $line.dir
	mv $i14 $line.dir

done < $1 