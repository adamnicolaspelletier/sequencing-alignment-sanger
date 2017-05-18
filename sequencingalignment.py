  
###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### scriptname.py
### Copyright 2016 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 29-03-16



from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
import glob  #module for file wildcards. can use genename as reference to look for several transcript fasta sequences at once. CHANGE PATH.
import os
import pandas as pd
import numpy as np
import fileinput
import operator
import subprocess
import pprint
import sys
import argparse


########################################################################################################################################################
########################################################## INSTRUCTIONS ################################################################################

# """ 
# 1. Have the FASTQ sequence of your sequencing run in the current directory
# 2. Also include the CCDS sequence in FASTA format of your gene of interest . Have ther EnsGeneID, transcript ID and geneSymbol in the header for each FASTA format
# 3. Have the filename information/primer/orientation in a textfile, refer to the template.txt files for formatting.
# 4. Edit the USER input and output section in this file for your files. 
# 5. Run the script.
# 6. To get a clustal readable output, run the muscle.sh script after this. In terminal run ./muscle.sh , which reads all fasta files in 
# 	clustal_alignment_genes dircetory and makes a clustalw alignemnt in the final
# 7. When done, remove all FASTA and CLUSTAL alignment files from the ABI_files, clustal_alignments_genes and final files.
# """
########################################################################################################################################################
########################################################## USER INPUT  AND OUTPUT ######################################################################



parser = argparse.ArgumentParser(description="""Aligns  FWD and REV contigs from Sanger sequencing run and aligns them to 
		their reference sequences using the MAFFT alogrithm. Can also integrate alignment to variant forms of the reference genes if FASTA sequences are supplied.
		Outputs MSA in both the FASTA format and CLUSTAL output, for easy visualization and validation""" )
parser.add_argument("-fq","--fastqfile",
					help="FASTQ file containing sequencing results. Defaults to 'FastQ/fastaSeq7776.txt'", default= "FastQ/fastaSeq7776.txt")
parser.add_argument("-r", "--reference", default="ccds_seq2.fasta",
					help="File Containing reference FASTA sequences. Defaults to 'ccds_seq2.fasta'")
parser.add_argument("-i", "--isoform", default="isoform_list2.txt",
					help="List of possible isoforms file. Defaults to 'isoform_list2.txt'")
parser.add_argument("-snv", "--snvfasta", default="SNVfasta.txt",
					help="File Containing Variant FASTA sequences. Defaults to 'SNVfasta.txt'")
parser.add_argument("-si", "--snvinfo", default="SDM_primer_output.txt",
					help="File Containing Additional Info on the Variants. Defaults to 'SDM_primer_output.txt'")
parser.add_argument("-f","--filenames",
					help="File containing the info for each individual FastQ Entry. Defaults to 'ABI_files/filenames3.txt'", default= "ABI_files/filenames3.txt")
parser.add_argument("-go","--gapopen", type=int, 
					help="""Value for the Gap Open Penalty for the Sequence Aligner. Higher values tends to prioritize 
					large blocks and facilitate SNV detection in highly similar sequences. Defaults to -100""", default= -100.0)
args = parser.parse_args()


fastqfile = args.fastqfile
ccdsfile = args.reference
isoformfile = args.isoform
snvfile = args.snvfasta
snpinfo = args.snvinfo
filenames = args.filenames
gapopenvalue = args.gapopen

print "Using %s for FASTQ ..." % fastqfile
print "Using %s for Reference FASTA ..." % ccdsfile
print "Using %s for Possible Isoforms File ..." % isoformfile
print "Using %s for Vairant Fasta sequences..." % snvfile
print "Using %s for Additional SNV Information ..." % snpinfo
print "Using %s for FastQ Entries List ..." % filenames
print "Using %i for FastQ Entries List ...\n\n" % gapopenvalue


########################################################################################################################################################
########################################################################################################################################################

def reversecomp(rprimsequence): ## make a complement version of the sequence, and reverse it so it has the proper orientation
	a = ""
	tempzrev = rprimsequence
	tempzrev = tempzrev.replace("T","X")
	tempzrev = tempzrev.replace("A","T")
	tempzrev = tempzrev.replace("X","A")
	tempzrev = tempzrev.replace("C","Y")
	tempzrev = tempzrev.replace("G","C")
	tempzrev = tempzrev.replace("Y","G")
	templist = list(tempzrev)
	templist.reverse()
	for i in templist:
		a += i
	return a

snv = raw_input("Include Variant in alignment? (y/n)")


if snv == "y":
	snvdict = {}
	snpinfodf = pd.read_csv(snpinfo, sep = "\t")
	snpinfodf2= snpinfodf[["Variation Name", "GENENAME", "Ensembl Gene ID"]]
	snplist = snpinfodf2["Variation Name"].tolist()
	mutseqlist = []

	print "Running Muscle alignment with variants..."
	for recordz in SeqIO.parse(snvfile, "fasta"):
		snvdict[str(recordz.id)] = str(recordz.seq)
	
	for i in snplist:
		try:
			mutseqlist.append(snvdict[i])
		except KeyError:
			mutseqlist.append("---")

	snpinfodf2["Sequence"] = mutseqlist

		# snpinfodf["Sequence"] = recordz.
		# fwd_id = snpinfodf.iloc[i]["Fwd_file"]
	# print snpinfodf2
elif snv == "n":
	print "Running Muscle alignment without variants..."

else: 
	print str(snv) + " is an invalid choice. Aborting."
	sys.exit()


# fasta_sequences = SeqIO.parse(open(input_file),'fasta')
# with open(output_file) as out_file:
#     for fasta in fasta_sequences:
#         name, sequence = fasta.id, fasta.seq.tostring()
#         new_sequence = some_function(sequence)
#         write_fasta(out_file)

# Generate the list of PWMs to use
# filelist = []
# for i in os.listdir("Homosapiens_CISBP/pwmsUnif"):
#   if i.endswith(".txt"):
#       filelist.append("Homosapiens_CISBP/pwmsUnif"+"/")



# promoterlist = []
# for i in os.listdir("promoter_files/seeded/"+str(n)+"bp/"):
#     if i.endswith(".txt"):
#         promoterlist.append("promoter_files/seeded/"+str(n)+"bp/"+i)

isoforms = open(isoformfile).read().splitlines()
outputfile = "ABI_files/fasta/tempoutput.txt"

phred_dict = {}

handle = open(fastqfile, "rU")

for record in SeqIO.parse(handle, "fastq"):
	phred_dict[record.id] = record.letter_annotations['phred_quality']



def fastaconvert(fastalist):  #convert a conventional fasta file into a list of IDs and whole sequences (merges the 50 characters per line)
	a = ""
	x = ""
	z = []
	for i in fastalist:
		if ">" in i:
			a = a.replace("\n","")
			z.append(a)
			x += i
			x = x.replace("\n","")
			z.append(x)
			a = ""
			x = ""
		else:
			a += i
	del z[0]
	return z

def idsplitter(ID):
	ID.pop()
	realID = []
	for i in ID:
		name = str(i)
		idsplit = name.split("|")
		realID.append(idsplit[2])
	return realID

def matrixconverter(seqmatrix, letters): 
	"""generates a DNA sequence string based on a NumPy matrix"""
	a = np.transpose(np.nonzero(np.transpose(seqmatrix))).tolist()
	seqstring = ""
	for i in a:
		seqstring += letters[i[1]]
	return seqstring

def matrixmaker(dnastring):   
	""" Generates a numpy array from a DNA sequence string, made from ones and zeros. 2D representation of a DNA sequence. Complements the matrixconverter function""" 
	matrix = np.zeros( (5, len(dnastring)) )
	index = []
	for i in range(len(dnastring)):
		index.append([lettersinv[list(dnastring)[i]], i]) 
	a = np.array(index)  
	matrix[a[:,0], a[:,1]] = 1
	return matrix

def pwmwalk(pwm, sequence, pos): 
	""" Tests each position of a SHORT DNA sequence as a numpy array against a pwm array. Returns the best alignment score. 
	The Pos argument lets you define over which portion of the pwm you want to align. Add 0 for the whole thing. 
	Useful when you add buffer random sequences on the site to align earlier than the 0 position for both array and sequence, or further than the length of the array"""
	
	alphapos = 0
	alphascore = -1000              
	for i in xrange(pos,len(pwm.transpose())-pos):
		
		try:
			betapos = i 
			betascore = np.sum(sequence * pwm[:,betapos:(betapos + len(sequence.transpose()))])
			#print betapos, betascore, alphapos, alphascore
			if betascore > alphascore:
				alphascore = betascore
				alphapos = betapos
		except ValueError:
			pass
	return [alphascore,alphapos]

def dicttofasta(dictio,filename):
	""" Takes a dictinary with keys as ID WITHOUT the > and fastasequences as values to generate a fasta file"""
	with open(filename, "w") as tempoutput:
			tempoutput.write("")
	for i in dictio:
		print i
		with open(filename, "a") as tempoutput:
			tempoutput.write(">" + str(i) + "\n" + dictio[i] + "\n")




outputfile2 = "internal_files/muscle2.txt"
outputfile3 = "internal_files/muscle3.txt"
outputfile4 = "internal_files/muscle4.txt"
#finaloutput = "output_clustal/%s.txt"

isoformfile = "isoform_list.txt"

letters = {0:"A",1:"C",2:"G",3:"T", 4:"N"}  # dictionary of indexes of each nucleotide for matrices
lettersinv = {"A":0,"C":1,"G":2,"T":3, "N":4}


filenamedf = pd.read_csv(filenames, sep = "\t")
wells = filenamedf["Well"].tolist()

def phredscores(fastqfile, sequenceid):
	"""Retrieves a list of phred scores for a given contig"""
	for record in SeqIO.parse(fastqfile, "fastq"):
		phred = (record.letter_annotations['phred_quality'])

""" Get phred score from last value and use it on -. Otherwise, a badly read letter will always 
be prefered to a -, even if that - is in a region of the contig with a higher confidence score. """
# for record in SeqIO.parse(fastqfile, "fastq"):


def mafftalignment(listofids):
	""" Starts a MUSCLE alignment for all IDs in the list. Returns files in CLUSTAL and FASTA format, along with a list of IDs specified by the user to save for later use"""
	for i in range(len(listofids)):
		
		contigdict = {}
		contigdictr = {}
		ccdsdict = {}
		

		fwd_id = filenamedf.iloc[i]["Fwd_file"]
		rev_id = filenamedf.iloc[i]["Rev_file"]
		ccdslist = []
		ccds_seqf = []

		genename = filenamedf.iloc[i]["Genename"]
		outputfile = filenamedf.iloc[i]["outputfile"]
		for recordc in SeqIO.parse(ccdsfile, "fasta"):
			if filenamedf.iloc[i]["EnsID"] in str(recordc.id):
				ccdsdict[str(recordc.id)] = str(recordc.seq.upper())
		
		
		for record in SeqIO.parse(fastqfile, "fastq"):
			if record.id == fwd_id:
				contigdict["FWD"]= str(record.seq.upper())
				
				fphred = (record.letter_annotations['phred_quality'])
				
			elif record.id == rev_id:
				contigdictr["REV"]= reversecomp(str(record.seq.upper()))
				rphred = (record.letter_annotations['phred_quality'])
				rphred.reverse()
				
		print ccdsdict
		nextseq = 0
		orient = 0
		muscledict = {0:0, 1:1, 2:2 , 3:3} 
		ccdsidlist = []
		canceldict = ccdsdict.copy()


		while muscledict[0] == 0:
			idlist = []
			ccdsindex = 4
			ccdsidlist = []
			ccdsdict2 = {}
			
			for j in ccdsdict:
				idlist.append(j)
				ccdsidlist.append(j)
			for j in ccdsidlist:
				muscledict[ccdsindex] = j
				ccdsdict2[str(ccdsindex) + ": " + j] = ccdsdict[j]
				ccdsindex += 1
			
			
			mergedict = ccdsdict2.copy()
			if muscledict[1] == 1:
				mergedict.update(contigdict)
			elif muscledict[1] == -1:
				mergedict.update(contigdictr)	
			
			dicttofasta(mergedict,"internal_files/fastaoutputf.fa")
			

			
			subprocess.call("mafft --auto --clustalout --preservecase --quiet internal_files/fastaoutputf.fa > internal_files/fastaoutputftemp.afa", shell=True)

			subprocess.call("mafft --auto --preservecase --quiet internal_files/fastaoutputf.fa > internal_files/fastaoutputf.afa", shell=True)


			with open("internal_files/fastaoutputftemp.afa", "r") as f:
				print f.read()

			musclelist = ["0:ACCEPT ALIGNMENT", "1:VERIFY REVERSE CONTIG", "2: RESTORE DELETED ISOFORMS", "3: REVERSE COMPLEMENT CONTIGS", "REPORT ALIGNMENT TO THE END"]

			ccdsindex = 2
			for i in ccdsidlist:
				musclelist.append(str(ccdsindex)+": REMOVE "+ str(i))
				ccdsindex += 1
				
			print "\n" 
			print musclelist
			print "\n" 

			choice = True
			canceldict.update(muscledict)
			optionchoice= raw_input("CHOOSE OPTION:  ")

			while int(optionchoice) > ccdsindex-1:
				print str(optionchoice) + " is not a valid choice. Try again"
				print "\n" 
				print musclelist
				print "\n" 
				optionchoice= raw_input("CHOOSE OPTION:  ")
				

			if optionchoice == str(0):
				muscledict[0] = 1
				choice == False
			elif int(optionchoice) == 1:
				val = muscledict[1]
				muscledict[1] = (val *-1)
				choice = False
					
			elif int(optionchoice) <= ccdsindex-1:
				
				del ccdsdict[muscledict[int(optionchoice)]]
				choice = False
			
			elif int(optionchoice) == 2:
				ccdsdict.update(canceldict)
				choice = False
			elif int(optionchoice) == 3:
				ccdsdict.update()
				choice = False
			elif int(optionchoice) == 4:
				ccdsdict.update()
				choice = False
		
		mergedictr = ccdsdict.copy()
		mergedict = ccdsdict.copy()
		mergedict.update(contigdict)
		mergedictr.update(contigdictr)

		dicttofasta(mergedict,"internal_files/fastaoutputf.fa")
		dicttofasta(mergedictr,"internal_files/fastaoutputr.fa")

		for record in SeqIO.parse("internal_files/fastaoutputr.fa", "fasta"):
			if "REV" in str(record.id):
				
				rseq = str(record.seq)
				rseqid = str(record.id)
				
		with open("internal_files/fastaoutputf.fa", "a") as f:
			f.write("\n>"+rseqid+"\n"+ rseq+"\n")

		subprocess.call("mafft --auto --preservecase --quiet internal_files/fastaoutputf.fa > internal_files/both.afa", shell=True)
		subprocess.call("mafft --auto --preservecase --quiet --clustalout internal_files/fastaoutputf.fa > internal_files/test.afa", shell=True)	
		# subprocess.call("mafft --auto --preservecase --quiet internal_files/fastaoutputr.fa > internal_files/fastaoutputr.afa", shell=True)

		# subMSA1 = 0
		# subMSA2 = 0

		# with open("internal_files/fastaoutputf.afa", "r") as f:
		# 	for line in f.readlines():
		# 		if line.count(">") > 0 :
		# 			subMSA1 += 1
		
		# with open("internal_files/fastaoutputr.afa", "r") as f:
		# 	for line in f.readlines():				
		# 		if line.count(">") > 0 :
		# 			subMSA2 += 1
		# sub1 = range(subMSA1)
		# sub2 = range(subMSA2)
		

		# sub1 = [x+1 for x in sub1]
		# sub2 = [x+(1+subMSA1) for x in sub2]

		# subMSAlist = []

		# subMSAlist.append('\t'.join(map(str, sub1)))
		# subMSAlist.append('\t'.join(map(str, sub2)))
		

		# with open("internal_files/subMSA.txt", "w") as f:
		# 	f.write("")

		# with open("internal_files/subMSA.txt", "a") as f:
		# 	for i in subMSAlist:
		# 		f.write(i + "\n")

		# filenamesin = ["internal_files/fastaoutputf.afa", "internal_files/fastaoutputr.afa"]
		# with open('internal_files/merged.afa', 'w') as outfile:
		# 	for fname in filenamesin:
		# 		with open(fname) as infile:
		# 			outfile.write(infile.read())


		# subprocess.call("mafft --merge internal_files/subMSA.txt internal_files/merged.afa  > internal_files/both.afa", shell=True)

		# alignment = AlignIO.convert("internal_files/both.afa", "fasta", "internal_files/test.afa", "clustal")


		# subprocess.call("mafft --merge internal_files/fastaoutputf.afa internal_files/fastaoutputr.afa  > internal_files/test.afa", shell=True)



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
		
		test = fseq
		
		for j in xrange(len(fseq)):
			if fseq[j] == "-":
				phred_fwd.append(0)
			else:
				phred_fwd.append(fphred.pop(0))

		for j in xrange(len(rseq)):
			if rseq[j] == "-":
				phred_rev.append(0)
			else:
				phred_rev.append(rphred.pop(0))
		
		
		consensus = ""
		for j in xrange(len(fseq)):
			
			
			if j < 20:
				consensus += fseq[j]
			elif refseq[j] == fseq[j] or refseq[j] == rseq[j] :
				consensus += refseq[j]
			elif fseq[j] == rseq[j]:
				consensus += fseq[j]
			elif phred_fwd[j] and phred_rev[j] < 15: #subjective threshold.  20 is a 1/100 chance of a wrong call, 10 is a 1 in 10. 
				fphredreg = phred_fwd[j-10:j+10]
				rphredreg = phred_rev[j-10:j+10]
				fphredavg = sum(fphredreg) / len(fphredreg)
				rphredavg = sum(rphredreg) / len(rphredreg)
				if fphredavg > rphredavg:
					consensus += fseq[j]
					
				elif fphredavg < rphredavg:
					consensus += rseq[j]
				
			elif phred_fwd[j] > phred_rev[j]:
				consensus += fseq[j]
			
			elif phred_fwd[j] < phred_rev[j]:
				consensus += rseq[j]
			else:
				consensus += "N"
		consensus2 = consensus.replace("-","")
		

		with open("internal_files/outputfile3.txt", "w") as foutput:	
			foutput.write("")

		with open("internal_files/outputfile3.txt", "a") as foutput:	
			foutput.write('>' + str(genename) + "_Consensus" + '\n' + str(consensus2)+"\n")

		# print snpinfodf2
		ensemblid = []
		for x in idlist:
			with open("internal_files/outputfile3.txt", "a") as foutput:	
				foutput.write('>' + x + '\n' + str(ccdsdict[x])+"\n")
		
			ensemblid.append(x.split("|")[0])

		
		print ensemblid
		if snv == "y":
			snvdf = snpinfodf2[snpinfodf2["Ensembl Gene ID"].isin(ensemblid)]
			for x in snvdf["Variation Name"]:
				with open("internal_files/outputfile3.txt", "a") as foutput:	
					foutput.write('>' + x + '\n' + str(snvdict[x])+"\n")


		fastacall = "mafft --auto --preservecase --quiet internal_files/outputfile3.txt > internal_files/outputfile4.afa"
			
		
		# fastacall = "muscle -in internal_files/outputfile3.txt  -out internal_files/outputfile4.txt -gapopen -200 -gapextend 20" 
		subprocess.call(fastacall , shell=True)

		referencedict = {}
		changedict = {}

		for record in SeqIO.parse("internal_files/outputfile4.afa", "fasta"):
			if ensemblid[0] in str(record.id):
				reference = str(record.seq)
			
			changedict[record.id] = list(record.seq)
		
		for x in changedict:
			upstream = True
			while upstream == True:
				for j in xrange(len(reference)):
					if reference[j] != "-":
						upstream = False
						break
						# print changedict
					else:
						changedict[x][j] = ""
		
		for x in changedict:
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
						# print changedict
					else:
						changedict[x][j] = ""

		for x in changedict:
			testint = changedict[x]
			changedict[x] = testint[::-1]	

		for j in changedict:
			concat = "".join(changedict[j])
			changedict[j] = concat

		dicttofasta(changedict,"internal_files/outputfile4.txt")

		
		fastacall = "mafft --auto --preservecase --quiet internal_files/outputfile4.txt > final/fasta/%s.afa" % outputfile
		clustalcall = "mafft --auto --preservecase --quiet --clustalout internal_files/outputfile4.txt > final/clustal/%s.afa" % outputfile
		

		subprocess.call(fastacall , shell=True)
		subprocess.call(clustalcall , shell=True)

		cfile = "final/clustal/%s.afa" % outputfile
		with open(cfile, "r") as f:
				print f.read()

		print "\n"
		print "\n"
		irrelevant = raw_input("Ready for next gene?  Press any key to continue...")			

savelist = []

mafftalignment(wells)

while len(savelist) > 0:
	save = savelist
	savelist = []
	musclealignment(save)