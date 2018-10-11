BatVI user manual:

1. INSTALLATION
===============

Prior to the installation, please ensure blastn, Picard tools, samtools,
bedtools and bwa programs are installed. BLAST can be downloaded from
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/, Picard tools from
broadinstitute.github.io/picard/, samtools from samtools.sourceforge.net/,
bedtools from code.google.com/p/bedtools/ and BWA can be downloaded from
http://bio-bwa.sourceforge.net/.

To install, unpack the BatVI package and change to the directory containing
the extracted files.  To compile the package run the script build.sh.
The following commands give the full procedure:
tar -zxvf batvi1.00.tar.gz
cd batvi1.00
./build.sh

2. CONFIGURATION FILE
=====================

Within the BatVI directory or the directory containing the files to be
processed, you need to create a file named batviconfig.txt. This file states
the location of several indexes and files required by BatVI. In the next
section, we will describe how to build these indexes.

	INDEX=<path_to_pathogen_batmis_index>
	PATHOGEN_BLAST_DB=<path_to_pathogen_blast_index>
	HG_BLAST_DB=<path_to_human_blast_index>
	HG_GENOME=<compressed_human_genome>
	HG_BWA=<path_to_human_BWA_index>
	PATHOGEN_BWA=<path_to_pathogen_BWA_index>
	BLAST_PATH=<path_to_blast_binary>
	BWA_PATH=<path_to_BWA_binary>
	PICARD_PATH=<path_to_picard_jarfiles>
	SAMTOOLS_PATH=<path_to_samtools_binary>
	BEDTOOLS_PATH=<path_to_bedtools_binary>

The content of these directories are as follows.
	INDEX is the batmis index of the virus database.
	HG_BLAST_DB is the blast index of the virus database.
	PATHOGEN_BLAST_DB is the blast index of the human genome.
	HG_GENOME contains a compressed version of the human genome.
	HG_BWA is the BWA index of human.
	PATHOGEN_BWA is the BWA index of the virus database.
	BLAST_PATH is the path to the blast binary
	BWA_PATH is the path to the BWA binary
	PICARD_PATH is the path to the picard jarfiles
	SAMTOOLS_PATH is the path to the samtools binary
	BEDTOOLS_PATH is the path to the bedtools binary

Below is an example file for batviconfig.txt

#This is an example batviconfig.txt file.
#Comments can be written preceded by a hash..
INDEX=/mnt/projects/HBVall/batmis/HBVall.fa
PATHOGEN_BLAST_DB=/mnt/projects/blast/HBVall/HBVall.fa
HG_BLAST_DB=/mnt/projects/blast/hg19/hg19.fa
HG_GENOME=/mnt/projects/hg19/hg19.fa
HG_BWA=/mnt/projects/bwa/hg19/hg19.fa
PATHOGEN_BWA=/mnt/projects/bwa/HBVall/HBVall.fa
BLAST_PATH=/usr/bin
BWA_PATH=/usr/bin
PICARD_PATH=/usr/bin/picard/
SAMTOOLS_PATH=/usr/bin
BEDTOOLS_PATH=/usr/bin

You can let BatVI automatically try to generate the entries for paths of the
binaries by running the script gen_hits.sh 

3. PREPARING THE INDEXES
========================
Before you prepare the indexes, you need to have two fasta files: (1) the
fasta file for the human genome (say hg19.fa) and (2) the fasta file for the
database of viruses (say HBVall.fa).

1.To build the batmis index of the virus database:
Copy the virus database fasta file to $INDEX, and run the command
BatMis-3.00/bin/build_index $INDEX
This will produce two files with extentions .pac and .ann.locations

2.To build the blast index of the virus database:
Copy the virus database fasta file to $PATHOGEN_BLAST_DB, and run the command
makeblastdb -in $PATHOGEN_BLAST_DB -dbtype nucl

3.To build the blast index of the human:
Copy the human fasta file to $HG_BLAST_DB, and run the command makeblastdb -in
$HG_BLAST_DB -dbtype nucl

4.To build the batmis index of the human:
Copy the human genome fasta files to $HG_GENOME, and run the command
BatMis-3.00/bin/bwtformatdb -p $HG_GENOME
This will produce two files with extentions .pac and .ann.locations

5.To build the BWA index of the human:
Copy the human genome fasta files to $HG_BWA, and run the command bwa index -a
bwtsw $HG_BWA

6.To build the BWA index of the virus database:
Copy the virus database fasta file to $PATHOGEN_BWA, and run the command bwa
index -a bwtsw $PATHOGEN_BWA

4. RUNNING THE PROGRAM
======================

Create a directory. In the directory, create a file named "filelist.txt" that
contains the forward read filename, reverse read filename and the insert size
separated by semicolons. (The forward read file and reverse read file can
either be in fastq format or in fastq.gz format.) Below is a sample
filelist.txt
A_1.fq;A_2.fq;800
B_1.fq;B_2.fq;800

Integrations are called with the command
    call_integratons.sh <processing_directory> [ options ]

The options are:
-l|--log <log file name> 	- Name of the log files to write processing
information to
-t|--threads <thread number> 	- Number of threads to use
-f|--filterdup			- Filter out duplicate reads

5.OUTPUT FILES AND FORMAT
=========================

The list of best predictions are reported in the file final_hits.txt. The
columns are explained below.
	LIB : Library name    
	Chr : Chromosome of the human integration    
	Human Pos : Location of the human integration        
	Sign : Orientation of the human integration    
	Viral Sign : Orientation of the viral integration      
	Viral Pos : Location of the viral integration        
	Read Count : Number of reads used in the prediction of integration. If
the entry is marked MSA, the integration has been found using assembly. 
	Split Reads : Number of split reads involved in the predicion of the
integration. A higher number indicreases the confidence of a prediction.     
	Uniquely Mapped Reads : Number of unique mappings to the human genome
involved in the prediction. A higher number increases the confidence of a
prediction. 
	Multiply Mapped Reads : Number of multiple mappings used in predicting
the integration. 
	Rank1 Hits : Number of reads which have the tophit near the
prediction. A higher number increases the confidence of a prediction.
	The last two columns of this file are to be ignored for this version
of BatVI.

The list of all possible breakpoints found by BatVI are listed in the file
predictions.opt.subopt.txt
The format of this file is as follows:
	Sign : Orientation of the human integration    
	Chr : Chromosome of the human integration    
	Human St : Location of the human integration        
	Human Ed : Location of the end of the human integration found by reads 
	Viral Sign : Orientation of the viral integration      
	Viral St : Location of the viral integration        
	Viral Ed : Location of the end of the viral integration found by reads
	Median_Rank : Median rank of the prediction
	Read Count : Number of reads used in the prediction of integration.
	Split Reads : Number of split reads involved in the predicion of the
integration. 
	Uniquely Mapped Reads : Number of unique mappings to the human genome
involved in the prediction. 
	Multiply Mapped Reads : Number of multiple mappings used in predicting
the integration. 
	Rank1 Hits : Number of reads which have the tophit near the
prediction. 
	The other columns of this file are to be ignored for this version of
BatVI.

The predictions using the assembly method are in the file
tmp.batvi/XXXX.predictionsx.msa.txt, where XXXX is the name of the directory
containing the fastq files.
	Chr : Chromosome of the human integration    
	Human Pos : Location of the human integration        
	Sign : Orientation of the human integration    
	Viral Pos : Location of the viral integration        
	Viral Sign : Orientation of the viral integration      
	Integration type : If the entry is marked TOPHIT, the assembly maps to
this location as the best hit. If the entry is marked REPEAT, the assembly can
map to at least one other location with similar confidence, and is therefore
ambiguous.

A list of clusters formed by BatVI and the reads belonging to them can be
found in the file \verb|t.opt.subopt.cluster|. The description of the entries
are given below.

	ID : Cluster ID 
	Chr : Chromosome read given in column 4 maps to. 
	Pos : Position read given in column 4 maps to. 
	Read ID : Read ID of a read belonging to the cluster. 

6. TEST EXAMPLE
===============

There is a set of test fastq files in the directory test. The filelist.txt
file is already present. To run a test, a HBV database fasta file is given in
test/HBVall.fa. You can download the hg19 fasta file from UCSC genome browser.
After you can run the program on this data set, the expected output is given
in test/expected.txt.

7. SPECIAL NOTES
================
 1) Please ensure  that the white spaces in the fasta files and genome names are removed or replaces with a character like an underscore.
 2) The index builder bwtformatdb might crash with a message like 
   "Building cached SA index...
    xxxx/build_index: 47 line: 168087 Segmentation fault      (core dump) ....". 
    This is OK as the required indexes would have already been built at this stage.
 3) The installation can be cleaned using the command 'build.sh clean '.
 4) A directory can be prepared for a fresh run of BatVI with the command 'clean_run.sh '.
