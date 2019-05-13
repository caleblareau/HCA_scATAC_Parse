#!/usr/bin/env python

# Caleb Lareau, Broad Institute

# This program will split a fastq file into several corresponding parts
# And extract reads with valid barcodes after automatically detecting
# Constant sequences unl
# By default, 5 million reads

##### IMPORT MODULES #####
import os
import re
import sys
import gzip
from optparse import OptionParser
from multiprocessing import Pool, freeze_support
from itertools import repeat

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes"
opts = OptionParser(usage=usage)
opts.add_option("-a", "--fastq1", help="<Read1> .fastq.gz of the left read on the fragment")
opts.add_option("-b", "--fastq2", help="<Read2> .fastq.gz of the right read on the fragment")

opts.add_option("-i", "--index1", help="<Read1> first index fastq.gz file returned")
opts.add_option("-j", "--index2", help="<Read2> second index fastq.gz file returned")

opts.add_option("-n", "--nreads", default = 5000000, help="Number of reads in each split output file / processed in a batch.")
opts.add_option("-c", "--ncores", default = 4, help="Number of cores for parallel processing.")
opts.add_option("--constant1", default = "", help="Constant sequence between Adapter 1 and Primer 1; if provided nothing, this is learned from the first 10k reads; default is compatible with Next-seq")
opts.add_option("--constant2", default = "CGACCGTTAGCAAAGCTCCG", help="Constant sequence between Adapter 2 and Primer 2; if provided nothing, this is learned from the first 10k reads; default is compatible with Next-seq")

opts.add_option("-o", "--out", default = "out", help="Output sample convention")

options, arguments = opts.parse_args()

a = options.fastq1
b = options.fastq2

filei = options.index1
filej = options.index2

n = int(options.nreads)
cpu = int(options.ncores)

user_constant1 = options.constant1
user_constant2 = options.constant2

o = options.out

# Parse input files
extension = a.split('.')[-1]
if extension == "fastq" or extension == "fq":
	sys.exist("Quitting... GZIP your .fastq files!")
elif extension == "gz":
	print("Found supplied .fastq.gz files")
else:
	sys.exit("ERROR! The input files (-a , -b) a *.fastq.gz")


# Define global variables
dumb = "N"*4 + "_" + "N"*4 + "_" + "N"*4 + "_" + "N"*4
ADAPTER1 = {'GCGATCTA': 'A101',
            'ATAGAGAG': 'A102',
            'AGAGGATA': 'A103',
            'TCTACTCT': 'A104',
            'CTCCTTAC': 'A105',
            'TATGCAGT': 'A106',
            'TACTCCTT': 'A107',
            'AGGCTTAG': 'A108',
            'GATTTCCA': 'A109',
            'ATCATGTT': 'A110',
            'TTTCATCA': 'A111',
            'AGTCCGAC': 'A112'}
ADAPTER2 = {'TCGCCTTA': 'A201',
            'CTAGTACG': 'A202',
            'TTCTGCCT': 'A203',
            'GCTCAGGA': 'A204',
            'AGGAGTCC': 'A205',
            'CATGCCTA': 'A206',
            'GTAGAGAG': 'A207',
            'CCTCTCTG': 'A208',
            'AGCGTAGC': 'A209',
            'CAGCCTCG': 'A210',
            'TGCCTCTT': 'A211',
            'TCCTCTAC': 'A212',
            'CAGATCCA': 'A213',
            'ACAAACGG': 'A214',
            'ACCCAGCA': 'A215',
            'CACCACAC': 'A217'}
PRIMER1 = {'GCGATCTA': 'P101',
           'ATAGAGAG': 'P102',
           'AGAGGATA': 'P103',
           'TCTACTCT': 'P104',
           'CTCCTTAC': 'P105',
           'TATGCAGT': 'P106',
           'TACTCCTT': 'P107',
           'AGGCTTAG': 'P108',
           'GATTTCCA': 'P109',
           'ATCATGTT': 'P110',
           'TTTCATCA': 'P111',
           'AGTCCGAC': 'P112',
           'GCTAGAAA': 'P113',
           'CTTGGTTA': 'P114',
           'CGATACAC': 'P115',
           'TTGATGGA': 'P116',
           'TGCACGAA': 'P117',
           'GGCAACCT': 'P118',
           'ACATAAGG': 'P119',
           'CGTTGCTG': 'P120',
           'ATTGAACC': 'P121',
           'ACGAATGT': 'P122',
           'TGGGAATC': 'P123',
           'GCAGTCCG': 'P124',
           'GAACGGCT': 'P125',
           'GACCCAAT': 'P126',
           'AGTATGCA': 'P127',
           'CCAAGCCC': 'P128',
           'GCCACGTC': 'P129',
           'AAATTTGC': 'P130',
           'GAGGCTGC': 'P131',
           'AACTCGGA': 'P132',
           'CTTAATGC': 'P133',
           'GTTATCGT': 'P134',
           'CCCGCAGG': 'P135',
           'AACAATCA': 'P136'}
PRIMER2 = {'TCGCCTTA': 'P201',
           'CTAGTACG': 'P202',
           'TTCTGCCT': 'P203',
           'GCTCAGGA': 'P204',
           'AGGAGTCC': 'P205',
           'CATGCCTA': 'P206',
           'GTAGAGAG': 'P207',
           'CCTCTCTG': 'P208',
           'AGCGTAGC': 'P209',
           'CAGCCTCG': 'P210',
           'TGCCTCTT': 'P211',
           'TCCTCTAC': 'P212',
           'CAGATCCA': 'P213',
           'ACAAACGG': 'P214',
           'ACCCAGCA': 'P215',
           'CCCAACCT': 'P216',
           'CACCACAC': 'P217',
           'GAAACCCA': 'P218',
           'TGTGACCA': 'P219',
           'AGGGTCAA': 'P220',
           'TGTCCGCG': 'P221',
           'ATATGGAA': 'P222',
           'AACGAATT': 'P223',
           'TCGACGCC': 'P224',
           'CACTTTGT': 'P225',
           'TTCAAGTA': 'P226',
           'GCTATCAC': 'P227',
           'AATCTACT': 'P228'}
           
#------------------------------

def prove_barcodeDictionary(bc, dictionary):
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''
	
	values = list(dictionary.keys())
	
	if(bc in values):
		return(dictionary[bc])
	else:
		eo = process.extractOne(bc, values)
		if(eo[1] >= 74): # 74 comes from 6/8... the score is the score homology
			return(dictionary[eo[0]])
		else:
			return("NNNN")

def batch_iterator(iterator, batch_size):
	"""
	Returns lists of tuples of length batch_size.
	"""
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch
			
def determineConstantSequence(fastqFile):
	"""
	Method to determine constant sequence based on prevalence
	in the first 10k sequences. These are different depending
	on which chemistry is used. For simplicity, call these mi-seq 
	and next-seq
	"""
	constants = dict()
	nPeek = 10000
	
	ita = batch_iterator(FastqGeneralIterator(gzip.open(fastqFile, "rt")), nPeek)
	for i, fastqdat in enumerate(ita):
		for (header, sequence, quality) in fastqdat:
			barcode_long = sequence
			
			# Determine sequence constant which is surrounded by 8bp barcodes
			const = (barcode_long)[8:(len(barcode_long)-8)]
			constants[const] = constants.get(const, 0) + 1
		break

	const_pick = sorted(constants, key=constants.get, reverse=True)[0]

	print("Parsing constant sequences in : " + fastqFile)
	print("Determined " + const_pick + " at " + str(constants.get(const)/nPeek *100) + "% in reads")
	return(const_pick)

# Determine constant sequences from index files
const1 = determineConstantSequence(filei)
const2 = determineConstantSequence(filej)

# Determine whether or not to use mi-seq or next-seq parsing
miseq = False
if(const2 =="CGGAGCTTTGCTAACGGTCG"):
	miseq = True
elif(const2 =="CGACCGTTAGCAAAGCTCCG"):
	miseq = False
else:
	print("WARNING: supplied sequence not clearly miseq or nextseq; assuming nextseq by default")

# Check for user-supplied constant sequence
if(user_constant1 != ""):
	print("User over-ride for constant sequence; now const1 is " + user_constant1)
if(user_constant2 != ""):
	print("User over-ride for constant sequence; now const2 is " + user_constant2)

def formatRead(title, sequence, quality):
	"""
	Takes three components of fastq file and stiches them together in a string
	"""
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))

def extractbarcodeseq_sciatac(barcodeI1, barcodeI2):
	'''
	Function to extract barcodes
	'''
	# Parse out sequence features and split based on constant sequences
	left = barcodeI1
	right = barcodeI2
	
	idxl = find_near_matches(const1, left, max_l_dist=3)
	idxr = find_near_matches(const2, right, max_l_dist=3)
	
	# Parse out barcodes if we can ID the constants
	try:
		if(idxl[0] and idxr[0]):
		
			# Deal with the left barcodes
			endA2 = idxl[0][0]
			startP2 = idxl[0][1]
			if(endA2 == 7):
				bc_a2 = "N" + left[0:7]
			elif(endA2 == 8):
				bc_a2 = left[0:8]
			elif(endA2 == 9):
				bc_a2 = left[1:9]
			else:
				bc_a2 = "NNNN"
			
			bc_p2 = left[startP2:(startP2+8)]
			
			# Determine barcode from coordinate matches
			bc_a2 = prove_barcodeDictionary(bc_a2, ADAPTER2)
			bc_p2 = prove_barcodeDictionary(bc_p2, PRIMER2)
			
			# Deal with the right barcodes which varies from technology to technology
			# Basically, assume that we are next-seq, but if it's determined that we
			# are actually mi-seq, then switch the bc_a1 and bc_p1
			
			endA1 = idxr[0][0]
			startP1 = idxr[0][1]
			if(endA1 == 7):
				bc_a1 = "N" + right[0:7]
			elif(endA1 == 8):
				bc_a1 = right[0:8]
			elif(endA1 == 9):
				bc_a1 = right[1:9]
			else:
				bc_a1 = "NNNN"
			
			bc_p1 = right[startP1:(startP1+8)]
			
			if(miseq): #switch before validation
				temp = bc_p1
				bc_p1 = bc_a1
				bc_a1 = temp
			
			# Determine barcode from designated known values
			bc_a1 = prove_barcodeDictionary(bc_a1, ADAPTER1)
			bc_p1 = prove_barcodeDictionary(bc_p1, PRIMER1)
			
			return(bc_a2 + "_" + bc_p2 + "_" + bc_a1 + "_" + bc_p1)
		else:
			return(dumb)
	except:
		return(dumb)
	

def debarcode_sciatac(quad):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = quad[0];  listRead2 = quad[1]
	indexRead1 = quad[2]; indexRead2 = quad[3]

	
	# parameters to return
	fq1 = ""
	fq2 = ""
	na2 = 0
	np2 = 0
	na1 = 0
	np1 = 0
	npass = 0
	nfail = 0
	
	# Grab attributes for biological reads
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	
	# Grab attributes for technical reads
	index_sequence1 = indexRead1[1]
	index_sequence2 = indexRead2[1]
	
	# Modify read title just in case for other processing
	title1 = title1.replace("_", "-")
	title2 = title2.replace("_", "-")
	
	barcode = extractbarcodeseq_sciatac(index_sequence1, index_sequence2)
	four = barcode.split("_")
	if("NNNN" in four):
		# Something went wrong
		nfail = nfail + 1
		if(barcode != dumb):
			if("NNNN" == four[0]):
				na2 = 1
			if("NNNN" == four[1]):
				np2 = 1
			if("NNNN" == four[2]):
				na1 = 1
			if("NNNN" == four[3]):
				np1 = 1
	else:
		npass = 1
		fq1 = formatRead("".join(four) + "_" + title1, sequence1, quality1)
		fq2 = formatRead("".join(four) + "_" + title2, sequence2, quality2)
	return([[fq1, fq2], [na1, na2, np1, np2, npass, nfail]])


def chunkWriterGzip(filename, what):
	'''
	Basic function to write a chunk of a fastq file
	to a gzipped file
	'''
	with gzip.open(filename, 'wt') as out_write:
				out_write.writelines(what)
	return(filename)

# Define variables to keep track of things that fail
na1 = 0
na2 = 0
np1 = 0
np2 = 0
npass = 0
nfail = 0

with gzip.open(a, "rt") as f1:
	with gzip.open(b, "rt") as f2:
		with gzip.open(filei, "rt") as i1:
			with gzip.open(filej, "rt") as i2:
		
				# Establish iterators
				it1 = batch_iterator(FastqGeneralIterator(f1), n)
				it2 = batch_iterator(FastqGeneralIterator(f2), n)
				
				it_i1 = batch_iterator(FastqGeneralIterator(i1), n)
				it_i2 = batch_iterator(FastqGeneralIterator(i2), n)
		
				# iterate over batches of length n
				for i, batch1 in enumerate(it1):
					batch2 = it2.__next__()
					batch_i1 = it_i1.__next__()
					batch_i2 = it_i2.__next__()
					
					output = o +  "-split" + str(i+1).zfill(3)
			
					# parallel process the barcode processing and accounting of failures.
					pool = Pool(processes=cpu)
					pm = pool.map(debarcode_sciatac, zip(batch1, batch2, batch_i1, batch_i2))
					pool.close()
			
					# Aggregate output
					fqs = list(map(''.join, zip(*[item.pop(0) for item in pm])))
					counts = list(map(sum, zip(*[item.pop(0) for item in pm])))
					na1 = na1 + counts[0]
					na2 = na2 + counts[1]
					np1 = np1 + counts[2]
					np2 = np2 + counts[3]
					npass = npass + counts[4]
					nfail = nfail + counts[5]
			
					# Export one chunk in parallel
					filename1 = output +'_1.fastq.gz'
					filename2 = output +'_2.fastq.gz'
			
					pool = Pool(processes=2)
					toke = pool.starmap(chunkWriterGzip, [(filename1, fqs[0]), (filename2, fqs[1])])
					pool.close()
			
with open(o + "-debarcode" + '.log', 'w') as logfile:
	# give summary statistics
	logfile.write("\nParsing read pairs:\n" + a + "\n" + b + "\n")
	logfile.write("\n"+str(npass)+" reads parsed with barcodes ("+str(round(npass/(npass+nfail)*100, 2))+"% success)\n")
	logfile.write("Total reads that failed: "+str(nfail)+"\n\n")
	logfile.write("\nOf reads that could not be parsed that had valid, detectable constants:\n")
	logfile.write(str(na1) + " had a bad A1 barcode sequence\n")
	logfile.write(str(na2) + " had a bad A2 barcode sequence\n")
	logfile.write(str(np1) + " had a bad P1 barcode sequence\n")
	logfile.write(str(np2) + " had a bad P2 barcode sequence\n")