  
###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### sequencingalignment.py
### Copyright 2017 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 25-06-16
### Version 2.01


from Bio import SeqIO
import os
import pandas as pd
import numpy as np
import subprocess
import sys
import argparse
import shutil
from adamP_BioTools.dna_tools import reversecomp

########################################################################################################################################################
########################################################## USER INPUT  AND OUTPUT ######################################################################
cwd = os.path.dirname(os.path.realpath(__file__))

pd.options.mode.chained_assignment = None  # default='warn

parser = argparse.ArgumentParser(description="""Aligns FWD and REV contigs from Sanger sequencing run and aligns them to 
		their reference sequences using the MAFFT algorithm. Can also integrate alignment to variant forms of the reference genes if FASTA sequences are supplied.
		Outputs MSA in both the FASTA format and CLUSTAL output, for easy visualization and validation""" )
parser.add_argument("-fq","--fastqfile",
					help="FASTQ file containing sequencing results. Defaults to 'docs/ex_fastq_seq.fq'", default= str(cwd)+"/docs/ex_fastq_seq.fq")

parser.add_argument("-r", "--reference", default=str(cwd)+"/docs/ex_ccds_seq.fa",
					help="File Containing reference FASTA sequences. Defaults to 'docs/ex_ccds_seq.fa'")

parser.add_argument("-vm", "--variantmode", action="store_true",
					help="Activates variant mode, to make alignments based on the variant fasta flag.")

parser.add_argument("-c", "--clip", action="store_true",
					help="Clip terminal bases before and after the REFERENCE sequence in the final output.")

parser.add_argument("-v", "--variantfasta", default=str(cwd)+"/docs/ex_variantfasta.fa",
					help="File Containing Variant FASTA sequences. Can be generated automatically using the sitedirmutagen.py script. Defaults to 'docs/ex_variantfasta.fa'")

parser.add_argument("-f","--filenames",
					help="File containing the info for each individual FastQ Entry. Defaults to 'docs/ex_filenames.txt'", default= str(cwd)+"/docs/ex_filenames.txt")

parser.add_argument("-o", "--outputdir", default=str(cwd)+"/docs/Alignments/",
                    help="Directory for Clustal and FASTA alignments files. Will create directory if it does not exist. Defaults to docs/Alignments/")

parser.add_argument("-mp","--minphred",
					help="Minimal threshold of quality to trust a phred score comparison between 2 reactions. Defaults to 15. 20 is a 1/100 miscall, and 10 a 1/10. ", default= 15)

args = parser.parse_args()


fastqfile = args.fastqfile
ccdsfile = args.reference
variantfile = args.variantfasta
filenames = args.filenames
outputdir = args.outputdir
variantmode = args.variantmode
minphred = args.minphred
clip = args.clip

print "\nsequencingalignment.py script for automating alignment of SANGER sequencing reactions  \n"
print " ** use the -h flag for detailed help on how to use this script **\n"

print "Using %s for FASTQ ..." % fastqfile
print "Using %s for Reference FASTA ..." % ccdsfile
print "Using %s for FastQ Entries List ..." % filenames
print "Using %s directory to for outfiles ..." % outputdir

if variantmode == True:
	print "\n\n"
	print "Running MAFFT alignment with variants..."
	print "Using %s for Variant Fasta sequences..." % variantfile
else: 
	print "Running MAFFT alignment without variants..."

print "\n\n"

########################################################################################################################################################
########################################################################################################################################################




def dicttofasta(dictio,filename):
	""" Takes a dictionary with keys as ID WITHOUT the > and fastasequences as values to generate a fasta file"""
	with open(filename, "w") as tempoutput:
			tempoutput.write("")
	for i in dictio:
		
		with open(filename, "a") as tempoutput:
			tempoutput.write(">" + str(i) + "\n" + dictio[i] + "\n")

def clipfasta(fastafile):
	changedict = {}

	for record in SeqIO.parse(fastafile, "fasta"):
		if "REFERENCE" in str(record.id):
			reference = str(record.seq)
		
		changedict[record.id] = list(record.seq) #add all aligned sequences to a dictionary

	for x in changedict:
		upstream = True ## clip 5' first: everytime reference has a - , remove that base from consensus and other sequences, until a base is found. 
		while upstream == True:
			for j in xrange(len(reference)):
				if reference[j] != "-":
					upstream = False
					break  # if no break, the loop will continue iterating in the actual reference sequence, and if there is a gap (like an insertion, it would delete them from consensus. not good)
				else:
					changedict[x][j] = ""
	
	for x in changedict: # invert all sequences to perform the same operation on the 3' portion. 
		testint = changedict[x]
		changedict[x] = testint[::-1]			
	
	referencer = reference[::-1]
	for x in changedict:
		downstream = True
		while downstream == True: 
			for j in xrange(len(reference)):
				if referencer[j] != "-":
					downstream = False
					break
					
				else:
					changedict[x][j] = ""

	for x in changedict:
		testint = changedict[x]
		changedict[x] = testint[::-1]	

	for j in changedict:
		concat = "".join(changedict[j])
		changedict[j] = concat

	return changedict



def mafftalignment(idindex,filenamedf,ccdsfile):
	""" Starts a MAFFT alignment for the supplied ID. Will look for additional info for that id in the filenames.
	 Returns files in CLUSTAL and FASTA format. """
	

	## Load the fastq id, sequence and phred score corresponding to the Sequenced reactions ID, as in the filenames file. 
	fwd_id = filenamedf.iloc[idindex]["Fwd_file"]
	rev_id = filenamedf.iloc[idindex]["Rev_file"]
	idname = filenamedf.iloc[idindex]["ID"]
	genename = filenamedf.iloc[idindex]["Genename"]
	ensemblid = filenamedf.iloc[idindex]["EnsID"]
	outputfile = filenamedf.iloc[idindex]["outputfile"]
	ccdslist = []
	ccds_seqf = []

	seqdict = {} #Dictionary that will contain FWD, REV and Reference info to make a dataframe. 
	
	seqdfhead = ["Number","ID", "Sequence", "REFERENCE"]
	for record in SeqIO.parse(fastqfile, "fastq"):
		if record.id == fwd_id:
			seqdict[0] = ["","FWD", str(record.seq).upper(), "N"]
			fphred = (record.letter_annotations['phred_quality']) #phred scores for Forward reaction
			
		elif record.id == rev_id:
			seqdict[1] = ["","REV", reversecomp(str(record.seq).upper()), "N"]
			rphred = (record.letter_annotations['phred_quality']) #phred scores for Reverse reaction
			rphred.reverse()

	ccdsindex = 5 #this index serves in generating the options for choosing isoforms. Choices 0 to 4 are taken by other functions.
	ccdsnumber = 2 # actual index number in the df
	for recordc in SeqIO.parse(ccdsfile, "fasta"):
		if filenamedf.iloc[idindex]["EnsID"] in str(recordc.id):
			seqdict[ccdsnumber] = [str(ccdsindex),str(recordc.id)+"|REFERENCE", str(recordc.seq).upper(), "Y"]
			ccdsindex +=1 
			ccdsnumber +=1

	seqdf = pd.DataFrame.from_dict(seqdict).T
	seqdf.columns = seqdfhead
	
	seqdf["concatID"] = seqdf["Number"] + ":" + seqdf["ID"] ## Generate a concatenated version of IDs for menu lists during filtering. Visual purposes
	seqdf["removeID"] = seqdf["Number"] + ": REMOVE " + seqdf["ID"]
	
	seqdfc = seqdf.copy() # Generate a copy of original df so we can delete rows as we filter isoforms, but can restore the original df if necessary

	
	
	mafftdict = {0:0, 1:1}  ## this dictionary serves to store variables needed for filtering and orienting contigs

	while mafftdict[0] == 0: #as long as this sequence has not been accepted as correct...
		ccdsidlist = list(seqdfc[seqdfc["REFERENCE"] == "Y"]["ID"])#list of reference sequence IDs
		if mafftdict[1] == 1:  #if 1, align ALL sequences from the df except the REV. 
			outdict = seqdfc[seqdfc["ID"] != "REV"].set_index("concatID")["Sequence"].to_dict()
		elif mafftdict[1] == -1: #if 1, align ALL sequences from the df except the FWD. 
			outdict = seqdfc[seqdfc["ID"] != "FWD"].set_index("concatID")["Sequence"].to_dict()
		## This portion is necessary for longer alignments,as errors in the extremities of either reactions tend to make them align improperly to
		## the reference sequence, as they also try to align unto one another, without success. This creates gaps etc, making it harder to align to isoforms. 
		## Aligning one, better allowing to vizualise the correct isoform, then the other if needed, makes it much easier. 
		## They will be aligned together after isoform selection. 

		
		dicttofasta(outdict,"internal_files/fastaoutputf.fa") 
		# MAFFT is made to align sequences from a file. BioPython can normally handle it directly, but it does not have all options and is a little awkward. 
		# So I prefer to output to file, then use MAFFT command line tools to align this file in both CLUSTAL output for visualization, and FASTA for further processing in parallel. 
		

		subprocess.call("mafft --auto --clustalout --preservecase --quiet internal_files/fastaoutputf.fa > internal_files/fastaoutputftemp.afa", shell=True)
		subprocess.call("mafft --auto --preservecase --quiet internal_files/fastaoutputf.fa > internal_files/fastaoutputf.afa", shell=True)


		with open("internal_files/fastaoutputftemp.afa", "r") as f:  #Let user see the result of this alignment in clustal format. 
			print f.read()

		mafftlist = ["0:ACCEPT ALIGNMENT", "1:VERIFY REVERSE CONTIG", "2: RESTORE DELETED ISOFORMS", "3: REVERSE COMPLEMENT CONTIGS", "4: REPORT ALIGNMENT TO THE END"]
		mafftlist = mafftlist + list(seqdfc[seqdfc["REFERENCE"] == "Y"]["removeID"]) #visual tool for user prompt
			
		print "\n" 
		print mafftlist
		print "\n" 
		
		optionchoice= raw_input("CHOOSE OPTION:  ")

		choice = False

		while choice == False: # choice validation loop. Repeat until acceptable answer is given, then execute corresponding code. 
			if int(optionchoice) not in range(ccdsindex):
				print str(optionchoice) + " is not a valid choice. Try again"
				print "\n" 
				print mafftlist
				print "\n" 
				optionchoice= raw_input("CHOOSE OPTION:  ")
			

			elif int(optionchoice) == 0: # if answer is 0, set mafftdict[0] to 1, which allows to escape the first while loop. Accepts alingment.
				mafftdict[0] = 1
				choice = True

				
			elif int(optionchoice) == 1: # sets mafftdict[1] value to its opposite, to see the reverse complement of the current sequencing reaction shown. 
				val = mafftdict[1]
				mafftdict[1] = (val *-1)
				choice = True
					
			elif int(optionchoice) >= 5 and len(ccdsidlist)>1:  #Delete one of the REFERENCE isoforms from alignment, as long as it<s not the last. 
				seqdfc= seqdfc[seqdfc["Number"] != str(optionchoice)]
				choice = True

			elif int(optionchoice) >= 5 and len(ccdsidlist)==1: #block deletion of last isoform
				print "Cannot delete last reference sequence! Accept sequence or Restore isoforms..."
				print "\n" 
				print mafftlist
				print "\n" 
				optionchoice= raw_input("CHOOSE OPTION:  ")
			
			elif int(optionchoice) == 2: #restore deleted isoforms from df. Try again, better luck next time!
				seqdfc = seqdf.copy()
				choice = True

			elif int(optionchoice) == 3: #Reverse complement your reactions, in case they were inverted, to see if they alig better. Useful when alignment is messed up, to determine if the problem is due to bad data or orientation
				fwd = [seqdfc.loc[seqdfc["ID"] == "FWD"]["Sequence"].index.values[0],reversecomp(seqdfc.loc[seqdfc["ID"] == "FWD"]["Sequence"].values[0])]
				seqdfc.set_value(fwd[0],"Sequence", fwd[1])
				rev = [seqdfc.loc[seqdfc["ID"] == "REV"]["Sequence"].index.values[0],reversecomp(seqdfc.loc[seqdfc["ID"] == "REV"]["Sequence"].values[0])]
				seqdfc.set_value(rev[0],"Sequence", rev[1])
				choice = True

			elif int(optionchoice) == 4: #When not sure, save this id for later. 
				savelist.append(idindex)
				choice = True

		outdict = seqdfc.set_index("ID")["Sequence"].to_dict()  #after validation, save the contents of the df to a fasta file and align it. 

	dicttofasta(outdict,"internal_files/fastaoutputf.fa")
	subprocess.call("mafft --auto --preservecase --quiet internal_files/fastaoutputf.fa > internal_files/both.afa", shell=True)
	subprocess.call("mafft --auto --preservecase --quiet --clustalout internal_files/fastaoutputf.fa > internal_files/test.afa", shell=True)	
	



### Score comparison. 
	### Verify, on aligned sequences, whether forward and reverse reactions have the same base. If they do, carry on. 
	### IF they do not, use PHRED scores to see which to choose with the higest confidence. 
	for record in SeqIO.parse("internal_files/both.afa", "fasta"):
		
		if "FWD" in str(record.id):
			fseq = str(record.seq)
	
		elif "REV" in str(record.id):
			rseq = str(record.seq)
			
		else:
			refseq = str(record.seq)
			refid = record.id
	
	 
	phred_fwd = []
	phred_rev = []
	
	# test = fseq
	### Phred scores lists from FASTQ files were not made to be aligned. thus, the gaps created by alignment need to be compensated by adding equivalent in the
	## phred scores lists. 
	for j in xrange(len(fseq)):
		if fseq[j] == "-":
			phred_fwd.append(0)
		else:
			phred_fwd.append(fphred.pop(0))  # if this position has a gap in fwd strand, add a 0. OTherwise, transfer the first item from the phred list and 

	for j in xrange(len(rseq)):
		if rseq[j] == "-":
			phred_rev.append(0)
		else:
			phred_rev.append(rphred.pop(0))
	
	


	consensus = ""  #build a consensus sequence after phred score comparison. 
	for j in xrange(len(fseq)):
		if j < 20: ## before the first 20 bases, base calls are unreliable in a  sanger reaction. add whichever letter, it won't matter. 
			consensus += fseq[j]
		elif refseq[j] == fseq[j] or refseq[j] == rseq[j] : 
			consensus += refseq[j]
			#this is a method used by biologists: if one of the discrepant bases is identical to the reference, 
			#it's safe to assume it to be correct. Can be commented out.
		elif fseq[j] == rseq[j]:
			consensus += fseq[j]
		elif phred_fwd[j] and phred_rev[j] < minphred: #if phred score is low for both, take average of the flanking 10 bases of that basecall to get an estimate. Normally, bad bases calls are found in regions, not alone.  
			fphredreg = phred_fwd[j-10:j+10]
			rphredreg = phred_rev[j-10:j+10]
			fphredavg = sum(fphredreg) / len(fphredreg)
			rphredavg = sum(rphredreg) / len(rphredreg)
			if fphredavg > rphredavg:
				consensus += fseq[j]
				
			elif fphredavg < rphredavg:
				consensus += rseq[j]
			
		elif phred_fwd[j] > phred_rev[j]: #if phred score of F is higher than R, assume F is correct. 
			consensus += fseq[j]
		
		elif phred_fwd[j] < phred_rev[j]:
			consensus += rseq[j]
		else:
			consensus += "N"
	consensus2 = consensus.replace("-","")   # remove the gaps from previous alignment, a new alignment with the consensus will be done. 
	

	outdict = {(str(idname) + "_"+ str(genename) + "_Consensus"):str(consensus2)}
	outdict.update(seqdfc[seqdfc["REFERENCE"]=="Y"].set_index("ID")["Sequence"].to_dict()) ## make a dict from consensus and reference sequences, export ot fasta. 



	if variantmode == True:   #if variant mode is activated, also add variant sequences to this dictionary. 
		variantdict = {}
		for record in SeqIO.parse(variantfile, "fasta"):
			if ensemblid in str(record.id):
				variantdict[str(record.id)] = str(record.seq)
	
		outdict.update(variantdict)


	dicttofasta(outdict,"internal_files/outputfile3.txt")	
	subprocess.call("mafft --auto --preservecase --quiet internal_files/outputfile3.txt > internal_files/outputfile4.afa" , shell=True) #alignment



	##### Removal of trailing ends section : everything before and afetr the reference sequence is of no potential biological interest
	
	if clip == True:
		changedict = clipfasta("internal_files/outputfile4.afa")
	else:
		changedict = {}
		for record in SeqIO.parse("internal_files/outputfile4.afa", "fasta"):
			changedict[str(record.id)] = str(record.seq)

	dicttofasta(changedict,"internal_files/outputfile4.txt")  #align all sequences to final output. 
	fastacall = "mafft --auto --preservecase --quiet internal_files/outputfile4.txt > %sFASTA/%s.afa" % (outputdir, outputfile)
	clustalcall = "mafft --auto --preservecase --quiet --clustalout internal_files/outputfile4.txt > %sCLUSTAL/%s.afa" % (outputdir, outputfile)
	subprocess.call(fastacall , shell=True)
	subprocess.call(clustalcall , shell=True)



	cfile = "%sCLUSTAL/%s.afa" % (outputdir,outputfile)
	with open(cfile, "r") as f:
			print f.read()

	print "\n"
	print "\n"

	
	pauseinput = False
	while pauseinput == False:
		pause = raw_input("Is this alignment correct?  Keep and move on (Y) or start over at the end(N)?  (Y/N):  ")
		if pause.upper() == "Y":
			pauseinput = True
			

			pass
		elif pause.upper() == "N":
			savelist.append(i)
			pauseinput = True

		else:
			print "%s is an incorrect answer. Try again" % pause


	print "Saved FASTA and CLUSTAL outputs to %s" % outputdir




### Start execution here,
## Determine if the user wants to align his sequencing data with associated supplied variants 


	
if variantmode == True: #load variant dictionary only if variant mode is on, to save memory. 
	variantdict = {}
	for record in SeqIO.parse(variantfile, "fasta"):
		variantdict[str(record.id)] = str(record.seq)
	
 # make necessary directories 
try:
    os.makedirs("internal_files/")
except OSError:
    if not os.path.isdir("internal_files/"):
        raise 
try:
    os.makedirs(outputdir)
except OSError:
    if not os.path.isdir(outputdir):
        raise 
try:
    os.makedirs(outputdir+"FASTA/")
except OSError:
    if not os.path.isdir(outputdir+"FASTA/"):
        raise 
try:
    os.makedirs(outputdir+"CLUSTAL/")
except OSError:
    if not os.path.isdir(outputdir+"CLUSTAL/"):
        raise 


filenamedf = pd.read_csv(filenames, sep = "\t")
ids = filenamedf["ID"].tolist()
savelist = []
for i in xrange(len(ids)):
	mafftalignment(i, filenamedf,ccdsfile)
	


while len(savelist) > 0:
	save = savelist
	savelist = []
	for i in save:
		mafftalignment(i,filenamedf,ccdsfile)

print "Alignment Over!" 


#remove the temporary files and directory for cleaning up. 
shutil.rmtree("internal_files")