#!usr/bin/env python3

import argparse
import sys
import os
import csv
import pysam
from datetime import datetime
import math

extended_help = """
To run:
python3 window_cnv.py --controls ../samples/down_syndrome/controls/subset/ --experimental ../samples/down_syndrome/ds/subset/

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
sys.stdout = open(datetime.now().strftime('%I:%M%p_%b%d') + '_cnv_bins_full.out', 'w')

# Chromosomes of interest
chromosomes = []
for i in range (1, 23): # 1..22
        chromosomes.append('chr' + str(i))
chromosomes.extend(["chrX", "chrY"])


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



def initCnv(bamfilePaths):
	"""
	Function: 
	read BAM files and store coverage for each chrm:bin:sample.
	
	Parameters:
	bamfilePaths (list): list containing all BAM file paths.
	
	Returns:
	cnv: dictionary containing cnv info of each chrm:bin:sample.
	bamfiles: list containing BAM files. 
	stats: dictionary containing the total reads and total read lengths for each sample.
	"""
	cnv = {chrm: {} for chrm in chromosomes}
	stats = {}
	print(chromosomes)
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
					
		stats[bamfile] = {key: 0 for key in [*chromosomes, "totalReads", "totalLen"]}	
		for read in samfile:
			# Skip PCR and optical duplicate reads
			if read.is_duplicate == False:
				chrm = read.reference_name
				readLength = read.reference_length
				pos = read.reference_start
				if chrm in cnv:
					# Increment stats components
					stats[bamfile][chrm] += 1
					stats[bamfile]['totalReads'] += 1				
					stats[bamfile]['totalLen'] += readLength		
					
					# Set binkey	
					binkey = int(pos/arg.window)
					if binkey in cnv[chrm]:
						if bamfile in cnv[chrm][binkey]:
						# binkey and sample exist
							cnv[chrm][binkey][bamfile]['sumReadCounts'] += 1
							cnv[chrm][binkey][bamfile]['sumReadLengths'] += (1*readLength)
						else:
						# binkey exists but add new sample
							cnv[chrm][binkey].update({bamfile: {\
								'sumReadCounts': 1,\
								'sumReadLengths': 1*readLength,\
								copynum: 0
							}})
					else: 
						# add new binkey and new sample
						cnv[chrm].update({binkey: {bamfile: {\
							'sumReadCounts': 1,\
							'sumReadLengths': 1*readLength,\
							copynum: 0
						}}})
				
		print('time for bamfile = %s' % (datetime.now() - tstart))
	return (cnv, bamfiles, stats)

def findSampleSex(bamfiles, stats):
	"""
	Function:
	Find the sex of input sample.
	
	Parameters:
	bamfiles: list containing the input sample files.
	stats: dictionary containing sample info.

	Return:
	stats: updated dictionary with sex info as a string ("male" or "female").
	numMaleSamples: number of male samples.
	numFemaleSamples: number of female samples.
	numMaleControl: number of male control samples.
	numFemaleControl: number of female control samples.
	"""

	numMaleSamples = 0
	numFemaleSamples = 0
	numMaleControl = 0
	numFemaleControl = 0

	for bamfile in bamfiles:
		# chrY / chrX reads ratio = ~11-12% for males, ~0.4-0.5% for females
		# set male cutoff to 6%
		# if (reads in chrY / reads in chrX) chrY is > 0.06, sample is male
		if stats[bamfile]["chrY"]/stat[bamfile]["chrX"] > 0.06: 
			stats[bamfile]["sex"] == "male"
			numMaleSamples += 1
			if "control" in bamfile:
				numMaleControl += 1
		else:
			stats[bamfile]["sex"] == "female"
			numFemaleSamples += 1
			if "control" in bamfile:
				numFemaleControl += 1

	return (stats, numMaleSamples, numFemaleSamples, numMaleControl, numFemaleControl)


def finalizeCnv(cnv, stats, numMaleSamples, numFemaleSamples, numMaleControl, numFemaleControl):
	"""
	Function: 
	Update "cnv" dictionary with normalized copy number estimations.
	
	Parameters:
	cnv: dictionary containing raw reads and lengths in each chrm:bin:sample.
	numControl: number of control BAM files (samples).
	
	Returns:
	cnv: updated dictionary containing copy number estimations in each chrm:bin:sample.
	"""
	for chrm in cnv:
		for binkey in cnv[chrm]:
			normfactorCov = 0
			normfactorMaleX = 0
			normfactorFemaleX = 0
			normfactorMaleY = 0
			normfactorFemaleY = 0

			for bamfile in cnv[chrm][binkey]:
				if "control" in bamfile:
					if bamfile in cnv[chrm][binkey].keys():
						# Coverage 	= sum of (reads * read length) in each bin / bin size
						# To normalize: divide by total coverage of sample = sum of (reads * read length) in sample
						binCoverage = cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window
						totalReadsAndLengths = stats[bamfile]['totalLen']		
						normalizedBinCoverage = binCoverage / totalReadsAndLengths
							
						# if bamfile is from female, stays the same, diploid X, no Y
						# if bamfile is from male, use normfactorSexChrm and normalize to 1x each
						
						sampleSex = stats[bamfile]["sex"]
						if chrm == "chrX":
							if sampleSex == "male":
								normfactorMaleX += normalizedBinCoverage
							if sampleSex == "female":
								normfactorFemaleX += normalizedBinCoverage
						elif chrm == "chrY":
							if sampleSex == "male":
								normfactorMaleY += normalizedBinCoverage
							if sampleSex == "female":
								normfactorFemaleY += normalizedBinCoverage
						else:
							normfactorCov += normalizedBinCoverage
						 	
												
	
			# Normalize coverage values to a copy number of ~2 for autosomal chromosomes (for normal diploid)	
			normfactorCov = (0.5 * normfactorCov) / (numMaleControl + numFemaleControl)
			# Normalize coverage value to a copy number of ~1 for sex chromosomes
			if numMaleControl != 0:
				normfactorMaleX = normfactorMaleX / numMaleControl	
				normfactorMaleY = normfactorMaleY / numMaleControl
			if numFemaleControl != 0:
				normfactorFemaleX = normfactorFemaleX / numFemaleControl
				normfactorFemaleY = normfactorFemaleY / numFemaleControl
			
			for bamfile in cnv[chrm][binkey]:
				if bamfile in cnv[chrm][binkey].keys():
					sampleSex = stats[bamfile]["sex"]
					if chrm == "chrX":
						if sampleSex == "male" and normfactorMaleX != 0:
							cnv[chrm][binkey][bamfile][copynum] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorMaleX
						if sampleSex == "female" and normfactorFemaleX != 0:
							cnv[chrm][binkey][bamfile][copynum] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorFemaleX
					elif chrm == "chrY":
						if sampleSex == "male" and normfactorMaleY != 0:
							cnv[chrm][binkey][bamfile][copynum] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorMaleY
						if sampleSex == "female" and normfactorFemaleY != 0:
							cnv[chrm][binkey][bamfile][copynum] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorFemaleY
					else:
						if normfactorCov != 0:
							cnv[chrm][binkey][bamfile][copynum] = \
						 		(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorCov
	return cnv



def printCoverage(cnv, samples):
	"""
	Function: 
	output to file.
	
	Parameters:
	cnv: dictionary containing CNV estimations
	samples: list contanining BAM files paths

	Returns:
	None
	"""
	with open(datetime.now().strftime('%I:%M%p_%b%d') + '_cnv_bins_full.txt', 'w', newline='') as outfile:
		outfile = csv.writer(outfile, delimiter='\t')
		samples.sort()
		outfile.writerow(['chr', 'start', 'end', *samples])
        
		for chrm in cnv:
			for binkey in cnv[chrm]:
				readsList = []
				copynumList = []
				for sample in sorted (cnv[chrm][binkey]):
					copynumList.append(cnv[chrm][binkey][sample][copynum])
				copynumList = [round(elem, 3) for elem in copynumList]

				if len(cnv[chrm][binkey]) == len(samples):
					outfile.writerow([
						chrm,
						binkey*arg.window,
						(binkey+1)*arg.window,
						*copynumList
					])
	return None



if __name__ == '__main__':
	timestart = datetime.now()	
	print('timestart = %s' % timestart)
	print(arg.controls)
	print(arg.experimental)
	print(arg.window)

	bamfilePaths, numControl, numExp = getBamfiles(arg.controls, arg.experimental)
	end_getBamfiles = datetime.now()
	print('getBamfiles() %s' % (end_getBamfiles - timestart))

	cnv, bamfiles, stats = initCnv(bamfilePaths)
	end_initCnv = datetime.now()
	print('initCnv() %s' % (end_initCnv - end_getBamfiles))

	stats, numMaleSamples, numFemaleSamples, numMaleControl, numFemaleControl = findSampleSex(bamfiles, stats)
	end_findSampleSex = datetime.now()
	print('findSampleSex() %s' % (end_findSampleSex - end_initCnv)
	print(stats)	

	cnv = finalizeCnv(cnv, stats, numMaleSamples, numFemaleSamples, numMaleControl, numFemaleControl)
	end_finalizeCnv = datetime.now()
	print('finalizeCnv() %s' % (end_finalizeCnv - end_findSampleSex))
		
	printCoverage(cnv, bamfiles)
	end_printCoverage = datetime.now()
	print('printCoverage() %s' % (end_printCoverage - end_finalizeCnv))

	print('timeend = %s' % (datetime.now() - timestart))
