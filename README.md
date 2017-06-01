#sequencing-alignment-sanger

#Scripts to handle Sanger sequencing data

Sanger sequencing is used very often in molecular biology. Yet, unless your desired sequence is short, if you try to align your forward and reverse runs against one-another, you do not have 100% of your sequence with good confidence. Sometimes, you will have inconsistencies/conflicting results (indels, SNVs) between the 2 sequencing reactions - what should you believe? 

Biologists tend to handle this manually in this case: visualize the sequences with software and look at the phred scores. And compare with the reference sequence to see what they should expect. But this is not good enough for high-throughput approaches: manually curating hundreds of sequencing runs is not only mind-numbing, but also more and more error-prone as time goes.

There are free common tools that perform alignments very well, but more complex options usually require $$$ software, and are not designed around high-throughput (isoform alignment, variant alignment, quality scores handling for discriminating. 

This is a simple approach to take FASTQ files and do global alignment using MAFFT command-line tools. When supplied with reference sequence for different isoforms, the script performs multiple alignment and returns the alignment in CLUSTAL format, asks for the user input to select the appropriate isoform, then generates a clustal output and fasta file. (sequencingalignment.py)


If necessary, the user can also supply fasta sequences with containing variants, which is useful to detect if a construct has acquired a desirable mutation after site-directed mutagenesis.

Since many sequencing platforms supply their results only in ABI and FASTA format, I also included a small script to transform ABI files to a single FASTQ file using a simple textfile legend. (abitofastq.py)


## Requirements:
 * Python 2.7+
 * Have MAFFT command line tools installed . Can be obtained at http://mafft.cbrc.jp/alignment/software/

 * Python Libraries :

	Pandas (+ Numpy)
	BioPython
	

 
 * Sequencing data. Either ABI files, or FASTQ files. FASTA can of course also be aligned, but has a loss in precision due to the absence of PHRED scores. 
 * Reference Sequence




## USER GUIDE:

This script was designed with high-throughput in mind: meaning if you only need one or 2 sequences aligned, you may find you are jumping through unnecessary hoops. However, it remains possible to use with this in mind. 

The example files included in the docs directory show how data is meant to be organized for the script to function properly. You do NOT need to use all files if you only intend to do general de novo alignment with your forward and reverse reactions, or do not need to determine if your mutagenesis experiment got you the intended mutation. Those are simply options. 

The flags have default values that use the example files included in the docs directory, so the script can be used without any flags. 
However, you will need to use the flags if you wish to use other files.


	python sequencingalignment.py -fq path/to/my/fastqfile.fq -f path/to/my/fastqinfofile.txt -r path/to/my/reference/sequences.fa



If you have more than one reference sequence (VERY LARGE FASTA file, with hundreds of entries and likely irrelevant sequences), the script will find the entries with the same Ensembl ID as your sequencing reaction form the filename file you supplied . 

	
However, if one wants to determine which isoform his amplified sequence corresponds to, thus more than 1 reference sequence, the script will show the alignment with all of them, letting the user choose which isoform he wants to align to. 


The script also includes a variant mode, which lets the user visualize his alignments with variant(SNV, mutations, indels) fasta sequences of the reference sequence. This is useful after cloning a potnetially mutated gene, or inducing directly mutagenesis. Variant FASTA sequences can be obtained form sitedirmutagen.py script, from the same author. 

	python sequencingalignment.py  -vm -v path/to/my/variantFASTA.fa 


Finally, the script has optional FASTA clipping functionality to clean out output and facilitate visualization, by removing unnecessary bases outside of the reference sequence. 

	python sequencingalignment.py  -c



## KNOWN BUGS:
 IPYTHON has problems with the argparse module, other Python distributions do not. I suggest using Official Python in the meantime. 
 Script has not been tested in Python 3.


## IN DEVELOPMENT:
 * Add functionality for multiple FORWARD and REVERSE contigs for a given reference sequence at once. Useful for long genes that require multiple reactions  
   to span the full length. 
 * Add support for regional phred scoring comparison, instead of per base. Will judge whether upon encounter of a discrepancy, whether the base call is a in 
   less reliable region of the sequencing reaction, such as often found in the extremities. 
 * Add optional automatic isoform detection, instead of user choice. Biologists tend to prefer the latter, but high-throughput is usually better suited with the former. 



 Use the -h flag for further options and documentation