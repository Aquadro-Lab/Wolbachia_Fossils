#!/bin/bash

# usage: bash nucmer_blast.sh SRAfilename.txt
# make sure we have the taxdb.btd and taxdb.bti in the working directory folder
# also make sure you have the BLAST nt database in the working directory folder. (Neither of these were working in the pipeline)

timestamp() {
	date +"%T"
}

timestamp

echo "$1"
echo "$2"

echo "------------ STEP 0.0: Copy BLAST database to local directory ------------"

#cp /shared_data/genome_db/BLAST_NCBI/nt* ./

# This didn't work for some reason on 10/29/23. idk why. If you copy it directly from the biohpc page, though, it works. 
# https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=16#c

timestamp

##echo "------------ STEP 0.1: Download NCBI taxdb to have taxids later ------------"
##
##
## Following this: https://www.biostars.org/p/76551/
## None of this worked on 10/29/23 on the cluster, so I just copied the files needed to the working directory on the cluster.
## This includes taxdb.bti and taxdb.btd files
##
##wget -r ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz #this is not working on the server on 10/29/23. Connection timed out error and Failed to connect to network error
##
##tar -xf ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz #thus commenitng out this
##
##rm -r ftp.ncbi.nlm.nih.gov #and this
##
## The instructions say: "As per the documentation make sure taxdb database is in the path defined by BLASTDB environment variable. Easiest is if you just make sure the files are in the same folder as the working directory"

timestamp

while IFS= read -r line; do

    echo "------------ SRA FILE: $line ------------"
    timestamp

echo "------------ STEP 1: Downloading reads with fasterq-dump ------------" 

##comment out this section if you already have the fastq reads because it takes a long time to download them.

#fasterq-dump --split-3 --progress -e 10 $line # UNCOMMENT THIS
#mkdir $line.dir # UNCOMMMENT THIS

pwd

##resume pipeline here if you commented out the above section. 

	echo "------------ STEP 2: setting variable names and testing them ------------" 

	SRRdir=$line.dir
	SRRname="${SRRdir%.dir}"

	fastq1="${SRRname}_1.fastq"
	fastq2="${SRRname}_2.fastq"

	fasta1="${SRRname}_1.fasta"
	fasta2="${SRRname}_2.fasta"

	fastq1stem="${SRRname}_1"
	fastq2stem="${SRRname}_2"

	delta1="${fastq1stem}.delta"
	delta2="${fastq2stem}.delta"

	coords1="${fastq1stem}.coords"
	coords2="${fastq2stem}.coords"

	SRRnumlistwithheader1="${SRRname}_1_withheaderlist.txt"
	SRRnumlistnoheader1="${SRRname}_1_noheaderlist.txt"
	SRRnumlist1="${SRRname}_1_list.txt" #this file has all the WGS reads that mapped to the Wolbachia pangenome

	SRRnumlistwithheader2="${SRRname}_2_withheaderlist.txt"
	SRRnumlistnoheader2="${SRRname}_2_noheaderlist.txt"
	SRRnumlist2="${SRRname}_2_list.txt" #this file has all the WGS reads that mapped to the Wolbachia pangenome

	SRRwolbachiafastq1="${SRRname}_1_wolbachia.fastq"
	SRRwolbachiafasta1="${SRRname}_1_wolbachia.fasta"
	SRRwolbachiafasta_noduplicates1="${SRRname}_1_wolbachia_noduplicates.fasta"

	SRRwolbachiafastq2="${SRRname}_2_wolbachia.fastq"
	SRRwolbachiafasta2="${SRRname}_2_wolbachia.fasta"
	SRRwolbachiafasta_noduplicates2="${SRRname}_2_wolbachia_noduplicates.fasta"

	blastresultstable1="${SRRname}_1_blastn_results.out"
	blastresultstable2="${SRRname}_2_blastn_results.out"

	results_noWol1="${SRRname}_1_blastn_results_noWol.txt"
	results_noWol2="${SRRname}_2_blastn_results_noWol.txt"

	results_onlydros1="${SRRname}_1_blastn_results_noWolonlyDros.txt"
	results_onlydros2="${SRRname}_2_blastn_results_noWolonlyDros.txt"

	results_onlyWol1="${SRRname}_1_blastn_results_onlyWol.txt"
	results_onlyWol2="${SRRname}_2_blastn_results_onlyWol.txt"

#	echo $SRRdir
#	echo $SRRname
#	echo $fastq1
#	echo $fastq2
#
#	echo $fasta1
#	echo $fasta2
#
#	echo $fastq1stem
#	echo $fastq2stem
#
#	echo $delta1
#	echo $delta2
#
#	echo $coords1
#	echo $coords2
#
#	echo $SRRnumlistwithheader1
#	echo $SRRnumlistnoheader1
#	echo $SRRnumlist1
#
#	echo $SRRnumlistwithheader2
#	echo $SRRnumlistnoheader2
#	echo $SRRnumlist2
#
#	echo $SRRwolbachiafastq1
#	echo $SRRwolbachiafasta1
#	echo $SRRwolbachiafasta_noduplicates1
#
#	echo $SRRwolbachiafastq2
#	echo $SRRwolbachiafasta2
#	echo $SRRwolbachiafasta_noduplicates2
#
#	echo $blastresultstable1
#	echo $blastresultstable2

	timestamp

	echo "------------ STEP 3: convert to fasta ------------" 

	/programs/seqtk/seqtk seq -a $fastq1 > $fasta1 #UNCOMMENT
	/programs/seqtk/seqtk seq -a $fastq2 > $fasta2 #UNCOMMENT

	timestamp

	echo "------------ STEP 4: nucmer ------------" 

	nucmer Wol_pan_genome95.fasta $fasta1 -p $fastq1stem
	nucmer Wol_pan_genome95.fasta $fasta2 -p $fastq2stem

	show-coords -r -c -l $delta1 > $coords1
	show-coords -r -c -l $delta2 > $coords2

	timestamp

	echo "------------ STEP 5: extract names of reads that map to Wolbachia pangenome in three simple steps ------------" 

	sed 's/.*SRR//p' $coords1 > $SRRnumlistwithheader1
	tail -n +7 $SRRnumlistwithheader1 > $SRRnumlistnoheader1
	sed -e 's/^/SRR/' $SRRnumlistnoheader1 > $SRRnumlist1

	sed 's/.*SRR//p' $coords2 > $SRRnumlistwithheader2
	tail -n +7 $SRRnumlistwithheader2 > $SRRnumlistnoheader2
	sed -e 's/^/SRR/' $SRRnumlistnoheader2 > $SRRnumlist2

	timestamp

	echo "------------ STEP 6: extract fasta sequences of reads identified in step 5 ------------"

	/programs/seqtk/seqtk subseq $fastq1 $SRRnumlist1 > $SRRwolbachiafastq1 #extract the fastq reads that we want
	/programs/seqtk/seqtk seq -a $SRRwolbachiafastq1 > $SRRwolbachiafasta1
	/programs/seqkit-0.15.0/seqkit rmdup -s < $SRRwolbachiafasta1 > $SRRwolbachiafasta_noduplicates1 #get rid of the duplicate fastas

	/programs/seqtk/seqtk subseq $fastq1 $SRRnumlist2 > $SRRwolbachiafastq2
	/programs/seqtk/seqtk seq -a $SRRwolbachiafastq2 > $SRRwolbachiafasta2
	/programs/seqkit-0.15.0/seqkit rmdup -s < $SRRwolbachiafasta2 > $SRRwolbachiafasta_noduplicates2 #get rid of the duplicate fastas

	timestamp

	echo "------------ STEP 7: BLAST the Wolbachia identified reads to the nr/nt database ------------"

	blastn -db nt -query $SRRwolbachiafasta_noduplicates1 -num_threads 4 -out $blastresultstable1 -outfmt '7 qseqid sseqid evalue bitscore sgi length pident mismatch gaps sacc staxids ssciname stitle'

	echo "completed for $SRRwolbachiafasta_noduplicates1"

	blastn -db nt -query $SRRwolbachiafasta_noduplicates2 -num_threads 4 -out $blastresultstable2 -outfmt '7 qseqid sseqid evalue bitscore sgi length pident mismatch gaps sacc staxids ssciname stitle'

	echo "completed for $SRRwolbachiafasta_noduplicates2"

	timestamp

	echo "------------ STEP 8: Parse BLAST results ------------"

	python parseBLASTtable_noWol.py $blastresultstable1 > $results_noWol1
	python parseBLASTtable_noWol.py $blastresultstable2 > $results_noWol2

	python parseBLASTtable_onlyDros.py $results_noWol1 > $results_onlydros1
	python parseBLASTtable_onlyDros.py $results_noWol2 > $results_onlydros2

	python parseBLASTtable_onlyWol.py $blastresultstable1 > $results_onlyWol1
	python parseBLASTtable_onlyWol.py $blastresultstable2 > $results_onlyWol2

	echo "------------ STEP 9: Move everything to the $line folder and delete unneccessary and/or large files ------------"

	rm $fastq1
	rm $fastq2

	rm $fasta1
	rm $fasta2

	mv $delta1 $line.dir
	mv $coords1 $line.dir

	mv $delta2 $line.dir
	mv $coords2 $line.dir

	rm $SRRnumlistwithheader1
	rm $SRRnumlistnoheader1

	rm $SRRnumlistwithheader2
	rm $SRRnumlistnoheader2

	mv $SRRnumlist1 $line.dir
	mv $SRRnumlist2 $line.dir

	mv $SRRwolbachiafasta1 $line.dir
	mv $SRRwolbachiafastq1 $line.dir
	mv $SRRwolbachiafasta_noduplicates1 $line.dir

	mv $SRRwolbachiafasta2 $line.dir
	mv $SRRwolbachiafastq2 $line.dir
	mv $SRRwolbachiafasta_noduplicates2 $line.dir

	mv $blastresultstable1 $line.dir
	mv $blastresultstable2 $line.dir

	mv $results_noWol1 $line.dir
	mv $results_noWol2 $line.dir

	mv $results_onlydros1 $line.dir
	mv $results_onlydros2 $line.dir

	mv $results_onlyWol1 $line.dir
	mv $results_onlyWol2 $line.dir

	timestamp
	
done < $1
