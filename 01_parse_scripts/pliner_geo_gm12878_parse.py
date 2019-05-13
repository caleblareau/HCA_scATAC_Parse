#!/usr/bin/env python

# Caleb Lareau, Broad Institute

# The following program will assign barcodes to biological reads
# from scTHS-seq data

# append established barcodes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4836442/#SD1

##### IMPORT MODULES #####
import os
import re
import sys
import gzip
import string
from optparse import OptionParser
import itertools
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from optparse import OptionParser
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes"
opts = OptionParser(usage=usage)
opts.add_option("-a", "--fastq1", help="<Read1> .fastq.gz of the biological read (to be preserved)")
opts.add_option("-b", "--fastq2", help="<Read2> .fastq.gz of the i7 read (to be assembled in barcode)")
opts.add_option("-c", "--index", help="Pain text file of the read IDs to be appended")
opts.add_option("-o", "--out", help="Output sample convention")
options, arguments = opts.parse_args()

print(options)
f1_in = options.fastq1
f2_in = options.fastq2
f3_in = options.index
outname = options.out

# Parse input files
extension = f1_in.split('.')[-1]
if extension == "fastq" or extension == "fq":
	sys.exist("Quitting... GZIP your .fastq files!")
elif extension == "gz":
	f1 = gzip.open(f1_in, 'rt')
	f2 = gzip.open(f2_in, 'rt')
	idx = gzip.open(f3_in, 'rt')
else:
	sys.exit("ERROR! The input files (-a , -b, -c) a *.fastq.gz")

# Define output file handles
out1 = gzip.open(outname + '_1.bc.fastq.gz', 'wt')
out2 = gzip.open(outname + '_2.bc.fastq.gz', 'wt')

# Parse input files
with open(f3_in) as f:
    barcodes = f.readlines()

barcodes = [x.strip() for x in barcodes] 

def matchKnownBarcode(bc):
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''
	if(bc in barcodes):
		return(bc)
	else:
		eo = process.extractOne(bc, barcodes)
		if(eo[1] >= 94): # 94 comes from 34/36... the score is the score homology
			return(eo[0])
		else:
			return("N")


nfail = 0
npass = 0	
##### Loop
while 1:

	# process the first file
	seqhead1 = f1.readline().rstrip()
	if not seqhead1: break
	seq1 = f1.readline().rstrip()
	qualhead1 = f1.readline().rstrip()
	qual1 = f1.readline().rstrip()

	# process the second file
	seqhead2 = f2.readline().rstrip()
	seq2 = f2.readline().rstrip()
	qualhead2 = f2.readline().rstrip()
	qual2 = f2.readline().rstrip()

	# Parse out barcodes if we can ID the constant
	barcode = seqhead1.split("_")[1].replace("+","").split(" ")[0]
	readname = seqhead1.split(" ")[0].strip("@")
	tru_barcode = matchKnownBarcode(barcode)
	
	if "N" == tru_barcode:
		nfail = nfail + 1
	else:
		npass = npass + 1
		seqhead1 = "@" + tru_barcode + "_" + readname
		qualhead1 = "+" + tru_barcode + "_" + readname
		seqhead2 = "@" + tru_barcode + "_" + readname
		qualhead2 = "+" + tru_barcode + "_" + readname
		
		# Write out to file of sequences that pass
		out1.write(seqhead1+"\n");out1.write(seq1+"\n")
		out1.write(qualhead1+"\n");out1.write(qual1+"\n")
		out2.write(seqhead2+"\n");out2.write(seq2+"\n")
		out2.write(qualhead2+"\n");out2.write(qual2+"\n")


out1.close(); out2.close()
f1.close();f2.close();idx.close()

with open(outname+ '.log', 'w') as logfile:
	# give summary statistics
	logfile.write("\nParsing read files:\n" + f1_in + "\n" + f2_in + "\n" + f3_in + "\n")
	logfile.write("\n"+str(npass)+" reads parsed with barcodes ("+str(round(npass/(npass+nfail)*100, 2))+"% success)\n")
	logfile.write(str(nfail)+" reads that could not be parsed\n")
	logfile.write("Total reads that failed: "+str(nfail)+"\n\n")

