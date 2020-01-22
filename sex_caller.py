#!usr/bin/env python3

import argparse
import sys
import os
import pysam
from datetime import datetime
from multiprocessing import Pool
import csv

extended_help = """
Input the path to the directory containing the bam files.
"""

# Input command line arguments
parser = argparse.ArgumentParser(
	description='sex caller',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)

parser.add_argument(
	'--samples',
	required=True,
	type=str,
	metavar='<path>',
	help='path to a directory of bam files')

parset.add_argument(
	'--cores',
	required=False,
	type=int,
	default = 1,
	metavar='<int>',
	help='number of cores to use for parallel processing')
	
arg = parser.parse_args()

# Write print statement outputs to file
# sys.stdout = open(datetime.now().strftime('%I:%M%p_%b%d') + '_sex_caller.out', 'w')

sexChrm = ["chrX", "chrY"]

def getBamfiles(path):
	bamfilePaths = []
	for bamfile in os.listdir(path):
		if bamfile.endswith('.bam'):
			bamfilePaths.append(path + bamfile)
	return (bamfilePaths)

def sexCaller(bamfile):
	print(bamfile)
	with open("cffDNA_sex_caller_output.txt", 'a', newline='') as outfile:
		outfile = csv.writer(outfile, delimiter='\t')
		stats = {}
		print('bamfile = %s' % bamfile)
		tstart = datetime.now()
		
		samfile = pysam.AlignmentFile(bamfile, 'rb') 
		stats[bamfile] = {"chrX": 0, "chrY":0}		
		for read in samfile:
			if read.reference_name in sexChrm and read.is_duplicate == False:
				stats[bamfile][read.reference_name] += 1
		
		filename = bamfile.split("/")[-1]
		samplename = filename.split("_")[0]	
		outfile.writerow([
			samplename,
			stats[bamfile]["chrX"],
			stats[bamfile]["chrY"],
			stats[bamfile]["chrY"]/stats[bamfile]["chrX"],
			round((stats[bamfile]["chrY"]/stats[bamfile]["chrX"]) * 100, 5)
			#,stats[bamfile]["sex"]
		])	
		
		print('time for bamfile = %s' % (datetime.now() - tstart))

# sort output file by Name?
# def sortOutput(outputFile):

	
if __name__ == '__main__':
	timestart = datetime.now()	
	print('timestart = %s' % timestart)

	# "../../cffDNA/bam/"
	bamfilePaths = getBamfiles(arg.samples)
	end_1 = datetime.now()
	print(*bamfilePaths, sep='\n')
	print(len(bamfilePaths))
	print('getBamfiles() %s' % (end_1 - timestart))

	with open("cffDNA_sex_caller_output.txt", 'w', newline='') as outfile:
		outfile = csv.writer(outfile, delimiter='\t')
		outfile.writerow(['Sample_name', 'ChrX_reads', 'ChrY_reads', 'ChrY:ChrX_ratio', 'ChrY:ChrX_%'])		
		
	# arg.cores = 56
	pool = Pool(processes = arg.cores) 
	pool.map(sexCaller, bamfilePaths)
	pool.close()
	pool.join()
	
	print('timeend = %s' % (datetime.now() - timestart))
