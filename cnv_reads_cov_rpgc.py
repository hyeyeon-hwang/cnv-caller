#!usr/bin/env python3

# control samples: JLKD, 7015, 7016, 6978, 6980
# exp samples: 7005, 7006, 7007, 7008

# https://medium.com/ngs-sh/coverage-analysis-from-the-command-line-542ef3545e2c
# vim tabs: https://stackoverflow.com/questions/1878974/redefine-tab-as-4-spaces
# set tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
# hi Comment ctermfg=LightGray

# scaling, standardizing
# https://www.theanalysisfactor.com/rescaling-variables-to-be-same/
# normalizing, scaling to range ** used this
# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
# normalization vs. standardization
# https://www.statisticshowto.datasciencecentral.com/normalized/
# bamCoverage: normalizeUsing
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
# R normalizing depth scoverage: has multiplying by a scale factor
# https://stackoverflow.com/questions/33892579/normalizing-depth-coverage-among-samples
# quantile normalization: use or no?
# https://en.wikipedia.org/wiki/Quantile_normalization
# github rpkm and rpgc formulas:
# https://github.com/deeptools/deepTools/wiki/Normalizations
# galaxy formulas and variable names for bamCoverage normalizeUsing
# https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=d0d9710a573c56e6&render_repository_actions_for=tool_shed&tool_config=%2Fsrv%2Ftoolshed%2Fmain%2Fvar%2Fdata%2Frepos%2F002%2Frepo_2134%2FbamCoverage.xml&changeset_revision=b61d934cc1e8

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

parser = argparse.ArgumentParser(
	description='CNV caller.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--controls', required=False, type=str, default='../samples/down_syndrome/final_controls/subsetChr21/',
	metavar='<path>', help='path to a directory of bam control files')
parser.add_argument('--experimental', required=False, type=str, default='../samples/down_syndrome/final_experimentals/subsetChr21/',
	metavar='<path>', help='path to a directory of bam files')
parser.add_argument('--window', required=False, type=int, default=5000,
	metavar='<int>', help='size of the window')
arg = parser.parse_args()

###
sys.stdout = open(datetime.now().strftime('%I:%M%p_Dec%d') + '_chr21_cov_rpgc.out', 'w')

# Chromosome of interest
chromosomes = []
for i in range (1, 23): #1..22
        chromosomes.append('chr'+str(i))
chromosomes.extend(['chrX', 'chrY'])
# validChrms = chromosomes
chromosomes = ['chr21'] 

# Get bam files from specified paths
def getBamfiles(controlPath, expPath):
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


def getReadsAndLength(bamfilePaths):
	stats = {} 
	for bamfile in bamfilePaths:
		print("bamfile = %s" % bamfile)
		samfile = pysam.AlignmentFile(bamfile, "rb")	
		stats[bamfile] = {"totalReads": 0, "totalLen": 0}
		for chrm in validChrms:
			stats[bamfile]["totalLen"] += samfile.get_reference_length(chrm)
			stats[bamfile]["totalReads"] += samfile.count(chrm)
	print(stats)
# normalized coverage = coverage in bin / avg number of reads in all bins = avg number of reads in entire chrm or samp?

# Read bam files and store coverage for each chrm:bin:sample
def initCoverage(bamfilePaths):
	# coverage: dictionary that stores for each chrm:binkey:sample the 
	# 1) rawLen: sum of read lengths
	# 2) rawCov: raw coverage = rawLen / 5000
	# 2) normSamp: raw coverage / total reads in sample
	# 3) normChrm: raw coverage / total reads in sample:chrm
	coverage = {chrm: {} for chrm in chromosomes}
	readStats = {} #key: 0 for key in chromosomes.append('total')}
	readStatsKeys = [*chromosomes, 'totalReads', 'totalLen']
	print(chromosomes)
	print(readStatsKeys)	
	bamfiles = []
	samfileStats = []
	
	for bamfile in bamfilePaths:	
		bamfiles.append(bamfile)
		print('bamfile = %s' % bamfile)
		tstart = datetime.now()
		samfile = pysam.AlignmentFile(bamfile, 'rb') 
		#print(samfile.mapped)
		#print(samfile.lengths)
		#print(samfile.references)
		
		#print(samfile.get_reference_length('chr1'))
		#print(samfile.get_index_statistics()[1])# find indices for chr1-22, X, Y, chr1 at index 0
					
		readStats[bamfile] = {key: 0 for key in readStatsKeys}	
		for read in samfile:
			if read.is_duplicate == False:

				chrm = read.reference_name
				readLength = read.reference_length
				pos = read.reference_start
				if chrm in coverage:
					# increment readStats
					readStats[bamfile][chrm] += 1
					readStats[bamfile]['totalReads'] += 1				
					readStats[bamfile]['totalLen'] += readLength		
			
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
								'reads': 0,\
								'cov': 0,\
								'rpgc': 0
							}})
					else: 
						# add new binkey and sample with:
						# [0] coverage normalized with total reads in sample = 1*read length for now
						# [1] coverage normalized with total reads for bin = 1*read length for now  
						# [2] total reads per bin per sample = 1 for initialization, increment
						# [3] if ds: exp / average control for coverage1 = 0 for now
						# [4] if ds: exp / average control for coverage2 = 0 for now
	
						coverage[chrm].update({binkey: {bamfile: {\
							'rawReadsNum': 1,\
							'rawReadsLen': 1*readLength,\
							'reads': 0,\
							'cov': 0,\
							'rpgc':0
						}}})
				
				
		print('time for bamfile = %s' % (datetime.now() - tstart))
	return (coverage, bamfiles, readStats)

def finalizeCoverage(coverage, readStats, numControl):
	effGenomeSize = 3088269832

	# control samples
	if "JLKD" in bamfile:
		totReads = 119091944
	if "6978" in bamfile:
		totReads = 114536755
	if "6980" in bamfile:
		totReads = 117716353
	if "7015" in bamfile:
		totReads = 82000514 
	if "7016" in bamfile:
		totReads = 80885569
	
	# experimental samples
	if "7005" in bamfile:
		totReads = 98612299
	if "7006" in bamfile:
		totReads = 108168785
	if "7007" in bamfile:
		totReads = 120873577
	if "7008" in bamfile:
		totReads = 27859548


	for chrm in coverage:
		for binkey in coverage[chrm]:
			normfactor = 0
			for bamfile in coverage[chrm][binkey]:
				if "control" in bamfile:
					if bamfile in coverage[chrm][binkey].keys():
						#normfactor = normfactor + coverage[chrm][binkey][bamfile]['rawReadsNum']/readStats[bamfile]['totalReads']
						normfactorReads += coverage[chrm][binkey][bamfile]['rawReadsNum']/totReads
						
						# using coverage / total reads
						normfactorCov += ( (coverage[chrm][binkey][bamfile]['rawReadsNum'] * coverage[chrm][binkey]['rawReadsLen']) / arg.window ) / totReads 
						normfactorRpgc += (coverage[chrm][binkey][bamfile]['rawReadsNum'] * effGenomeSize) / (totReads * coverage[chrm][binkey][bamfile]['rawReadsLen']) 
								
			normfactor = (0.5 * normfactor) / numControl
			for bamfile in coverage[chrm][binkey]:

				if bamfile in coverage[chrm][binkey].keys():
					if normfactor != 0:
						#coverage[chrm][binkey][bamfile]['copyNum'] = \
						#	(coverage[chrm][binkey][bamfile]['rawReadsNum'] / readStats[bamfile]['totalReads']) / normfactor 			
						
						coverage[chrm][binkey][bamfile]['reads'] = \
							(coverage[chrm][binkey][bamfile]['rawReadsNum'] / totReads) / normfactorReads 			
				
						coverage[chrm][binkey][bamfile]['cov'] = \
							( (coverage[chrm][binkey][bamfile]['rawReadsNum'] * coverage[chrm][binkey]['rawReadsLen']) / arg.window ) / \
							totReads / normfactorCov
				
						coverage[chrm][binkey][bamfile]['rpgc'] = \
							(coverage[chrm][binkey][bamfile]['rawReadsNum'] * effGenomeSize) / \
							(totReads * coverage[chrm][binkey][bamfile]['rawReadsLen']) / normfactorRpgc 
	return coverage

def printCoverage(coverage, samples):
	with open(datetime.now().strftime('%I:%M%p_Dec%d') + '_chr21_cov_rpgc.txt', 'w', newline='') as full:
		wfull = csv.writer(full, delimiter='\t')
		samples.sort()
		wfull.writerow(['chr', 'start', 'end', *samples])
        
		for chrm in coverage:
			for binkey in coverage[chrm]:
				readsList = []
				covList = []
				rpgcList = []
				for sample in sorted (coverage[chrm][binkey]):
					readsList.append(coverage[chrm][binkey][sample]['reads'])
					covList.append(coverage[chrm][binkey][sample]['cov'])
					rpgcList.append(coverage[chrm][binkey][sample]['rpgc'])			
				readsList = [round(elem, 3) for elem in readsList]
				covList = [round(elem, 3) for elem in covList]
				rpgcList = [round(elem, 3) for elem in rpgcList]

				if len(coverage[chrm][binkey]) == len(samples):
					wfull.writerow([
						chrm,
						binkey*arg.window,
						(binkey+1)*arg.window,
						*readsList,
						*covList,
						*rpgcList,	
					]) 

# print output for full bins too

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

#	getReadsAndLength(bamfilePaths)

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
