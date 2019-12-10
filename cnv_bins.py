#!usr/bin/env python3

import argparse
import sys
import os
import csv
import pysam
from datetime import datetime
import math

extended_help = """
This is the extended help.

To run:
python3 window_cnv.py --controls ../samples/down_syndrome/controls/subset/ --experimental ../samples/down_syndrome/ds/subset/

To subset bam file:
[https://bioinformatics.stackexchange.com/questions/3565/subset-smaller-bam-to-contain-several-thousand-rows-from-multiple-chromosomes]
samtools view -bo 7015.subset.bam -s 123.01 7015.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chrM chrX chrY

samtools view -b input.bam chrM chr21 > output.bam

For final:
'../samples/down_syndrome/final_controls/'
'../samples/down_syndrome/final_experimentals/'

For testing: 
'../samples/down_syndrome/final_controls/subsetChr21/'
'../samples/down_syndrome/final_experimentals/subsetChr21/'
"""

# Input command line arguments
parser = argparse.ArgumentParser(
	description='CNV caller.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument(
	'--controls', 
	required=False, 
	type=str, 
	default='../hh_cnv_caller/samples/down_syndrome/final_controls/',
	metavar='<path>', 
	help='path to a directory of bam control files')
parser.add_argument(
	'--experimental', 
	required=False, 
	type=str, 
	default='../hh_cnv_caller/samples/down_syndrome/final_experimentals/',
	metavar='<path>', 
	help='path to a directory of bam files')
parser.add_argument(
	'--window', 
	required=False, 
	type=int, 
	default=5000,
	metavar='<int>', 
	help='size of the window')
arg = parser.parse_args()

# Write print statement outputs to file
sys.stdout = open(datetime.now().strftime('%I:%M%p_%m%d') + '_cnv_bins_full_no_XY.out', 'w')

# Chromosomes of interest (only autosomal for now)
# Need to work on integrating sex chromosomes
chromosomes = []
for i in range (1, 23): # 1..22
        chromosomes.append('chr' + str(i))



def getBamfiles(controlPath, expPath):
	"""
	Function:
	get BAM files from specified paths.
	
	Parameters:
	controlPath (str): path to directory containing control BAM files.
	expPath (str): path to directory containing experimental BAM files.
	
	Returns:
	bamfilePaths (list): list containing all BAM file paths.
	numControl (int): number of control BAM files.
	numExp (int): number of experimental BAM files.   
	"""
	bamfilePaths = []
	numControl = 0
	numExp = 0
	for bamfile in os.listdir(controlPath):
		if bamfile.endswith('.bam'):
			bamfilePaths.append(controlPath + bamfile)
			numControl += 1
	for bamfile in os.listdir(expPath):
		if bamfile.endswith('.bam'):
			bamfilePaths.append(expPath + bamfile)
			numExp += 1	
	return (bamfilePaths, numControl, numExp)



def initCoverage(bamfilePaths):
	"""
	Function: 
	read BAM files and store coverage for each chrm:bin:sample.
	
	Parameters:
	bamfilePaths (list): list containing all BAM file paths.
	
	Returns:
	coverage: dictionary containing coverage of each chrm:bin:sample.
	bamfiles: list containing BAM files. 
	readStats: dictionary containing the total reads and total read lengths for each sample.
	"""
	coverage = {chrm: {} for chrm in chromosomes}
	readStats = {}
	print(chromosomes)
	print(readStatsKeys)	
	bamfiles = []
	samfileStats = []
	
	for bamfile in bamfilePaths:	
		bamfiles.append(bamfile)
		print('bamfile = %s' % bamfile)
		tstart = datetime.now()
		samfile = pysam.AlignmentFile(bamfile, 'rb') 
		# print(samfile.mapped)
		# print(samfile.lengths)
		# print(samfile.references)
		# print(samfile.get_reference_length('chr1'))
					
		readStats[bamfile] = {key: 0 for key in [*chromosomes, "totalReads", "totalLen"]}	
		for read in samfile:
			# Skip duplicate reads
			if read.is_duplicate == False:
				chrm = read.reference_name
				readLength = read.reference_length
				pos = read.reference_start
				if chrm in coverage:
					# Increment readStats components
					readStats[bamfile][chrm] += 1
					readStats[bamfile]['totalReads'] += 1				
					readStats[bamfile]['totalLen'] += readLength		
					
					# Set binkey	
					binkey = int(pos/arg.window)
					if binkey in coverage[chrm]:
						if bamfile in coverage[chrm][binkey]:
						# binkey and sample exist
							coverage[chrm][binkey][bamfile]['rawReadsNum'] += 1
							coverage[chrm][binkey][bamfile]['rawReadsLen'] += (1*readLength)
						else:
						# binkey exists but add new sample
							coverage[chrm][binkey].update({bamfile: {\
								'rawReadsNum': 1,\
								'rawReadsLen': 1*readLength,\
								'cov': 0
							}})
					else: 
						# add new binkey and new sample with:
						#	rawReadsNum
						#	rawReadsLen
						#	cov	
						coverage[chrm].update({binkey: {bamfile: {\
							'rawReadsNum': 1,\
							'rawReadsLen': 1*readLength,\
							'cov': 0
						}}})
				
		print('time for bamfile = %s' % (datetime.now() - tstart))
	return (coverage, bamfiles, readStats)



def finalizeCoverage(coverage, readStats, numControl):
	"""
	Function: 
	Update "coverage" dictionary with normalized copy number estimations.
	
	Parameters:
	coverage: dictionary containing raw reads and lengths in each chrm:bin:sample.
	numControl: number of control BAM files (samples).
	
	Returns:
	coverage: updated dictionary containing copy number estimations in each chrm:bin:sample.
	"""
	effGenomeSize = 3088269832
	for chrm in coverage:
		for binkey in coverage[chrm]:
			normfactorCov = 0
			for bamfile in coverage[chrm][binkey]:
				if "control" in bamfile:
					if bamfile in coverage[chrm][binkey].keys():
						# Coverage 	= sum of (reads * read length) in each bin / bin size
						# To normalize: divide by total coverage of sample = sum of (reads * read length) in sample
						normfactorCov += ( coverage[chrm][binkey][bamfile]['rawReadsLen'] / arg.window ) / \
							 readStats[bamfile]['totalLen'] 
							
			# center to the median instead of the mean
			# instead of each bin, center to median of ALL autosomes and use that value as the reference copy number of ~2?
			# or better to do it for each bin? or better yet, for each sample. 
			# Find reference copy number of each sample by finding mean and median of all autosomes in all control samples 	
	
			# Normalize coverage values to a copy number of ~2 (for normal diploid)	
			normfactorCov = (0.5 * normfactorCov) / numControl
			for bamfile in coverage[chrm][binkey]:
				
				if bamfile in coverage[chrm][binkey].keys():
					if normfactorCov != 0:
						coverage[chrm][binkey][bamfile]['cov'] = \
							(coverage[chrm][binkey][bamfile]['rawReadsLen'] / arg.window ) / \
							readStats[bamfile]['totalLen'] / normfactorCov
				
	return coverage



def printCoverage(coverage, samples):
	"""
	Function: 
	output to file.
	
	Parameters:
	coverage: dictionary containing CNV estimations
	samples: list contanining BAM files paths

	Returns:
	None
	"""
	with open(datetime.now().strftime('%I:%M%p_%m%d') + '_cnv_bins_full_no_XY.txt', 'w', newline='') as fullC:
		fullC = csv.writer(fullC, delimiter='\t')
		samples.sort()
		fullC.writerow(['chr', 'start', 'end', *samples])
        
		for chrm in coverage:
			for binkey in coverage[chrm]:
				readsList = []
				covList = []
				for sample in sorted (coverage[chrm][binkey]):
					covList.append(coverage[chrm][binkey][sample]['cov'])
				covList = [round(elem, 3) for elem in covList]

				if len(coverage[chrm][binkey]) == len(samples):
					fullC.writerow([
						chrm,
						binkey*arg.window,
						(binkey+1)*arg.window,
						*covList
					])
	return None



if __name__ == '__main__':
	timestart = datetime.now()	
	print('timestart = %s' % timestart)
	print(arg.controls)
	print(arg.experimental)
	print(arg.window)

	s_getBamfiles = datetime.now()	
	bamfilePaths, numControl, numExp = getBamfiles(arg.controls, arg.experimental)
	e_getBamfiles_s_initCoverage = datetime.now()
	print('getBamfiles() %s' % (e_getBamfiles_s_initCoverage - s_getBamfiles))

	coverage, bamfiles, readStats = initCoverage(bamfilePaths)
	e_initCoverage_s_finalizeCoverage = datetime.now()
	print('initCoverage() %s' % (e_initCoverage_s_finalizeCoverage - e_getBamfiles_s_initCoverage))
	print(readStats)	

	coverage = finalizeCoverage(coverage, readStats, numControl)
	e_finalizeCoverage_s_printCoverage = datetime.now()
	print('finalizeCoverage() %s' % (e_finalizeCoverage_s_printCoverage -  e_initCoverage_s_finalizeCoverage))
		
	printCoverage(coverage, bamfiles)
	e_printCoverage = datetime.now()
	print('printCoverage() %s' % (e_printCoverage - e_finalizeCoverage_s_printCoverage))

	print('timeend = %s' % (datetime.now() - timestart))
