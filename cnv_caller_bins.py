#!/usr/bin/env python3
import sys
from datetime import datetime
import pysam
import pandas as pd
import csv
import multiprocessing as mp
import math
import numpy as np
from itertools import repeat
import os

wdir = os.getcwd()
outDir = wdir + '/cnv_caller_output/'

sys.stdout = open(datetime.now().strftime('%I:%M%p_Nov%d') + '_prints.txt', 'w')

samfiles = []
# Multiprocessing shared objects
manager = mp.Manager()
positionsAll = manager.dict()
coverageAll = manager.dict()

chrmOfInterest = []

for i in range(1, 23):
	chrmOfInterest.append('chr'+str(i))
chrmOfInterest.extend(['chrX', 'chrY', 'chrM'])	
#chrmOfInterest = ['chrM', 'chrY']
binSize = 5000
	
#for chrmKey in chrmOfInterest:
#	coverageAll[chrmKey] = {}	
	
def pileupCountsLoop(sampleIndex):
	for chrmKey in chrmOfInterest:
		if chrmKey not in coverageAll:
			coverageAll[chrmKey] = manager.dict()

		startBin = positionsAll[sampleIndex][chrmKey][0]
		endBin = positionsAll[sampleIndex][chrmKey][1]
		print('pileupLoops(sampleIndex = %s) %s startBin:endBin = %s:%s' % (sampleIndex, chrmKey, startBin, endBin))
		firstBin = math.floor(startBin/binSize) # *binSize
		lastBin = math.ceil(endBin/binSize) # *binSize
		
		for binKey in range(firstBin, lastBin):
			binCov = 0
			for col in samfiles[sampleIndex].pileup(chrmKey, binKey*binSize, (binKey+1)*binSize):
				binCov += col.nsegments
			binCov /= binSize
			
			if binCov > 1:	
				if binKey not in coverageAll[chrmKey]:
					coverageAll[chrmKey].update({binKey: manager.dict()})			
				coverageAll[chrmKey][binKey].update({sampleIndex: binCov})
	print('end of pileupCountsLoop(%s) = %s' % (sampleIndex, datetime.now() - startTime))

def initPositionsAll():
	for sampleIndex in range(len(samfiles)):
		positionsAll[sampleIndex] = manager.dict()
		# positionsAll[sampleIndex].update({chrm:manager.list() for chrm in chrmOfInterest})
		positionsAll[sampleIndex].update({chrm:[] for chrm in chrmOfInterest})
		# works without manager.list(), why?		

def findChrmBoundaries(sampleIdx, chrm):
	posList = []
	for read in samfiles[sampleIdx].fetch(chrm, multiple_iterators=True):
		posList.append(read.reference_start)
	startBin = posList[0] 
	endBin = posList[-1] 
	print('findChrmBoundaries(sampleIdx = %s, chrm = %s) = (startBin, endBin) = (%s, %s)' % (sampleIdx, chrm, startBin, endBin))
	positionsAll[sampleIdx][chrm] = [startBin, endBin]

# 4 hours to write to file
# split by chrm - to temp files for each chrm
# then merge files and delete temp files
#def outputToFile():
#	with open(datetime.now().strftime('%I:%M%p_Nov%d')+'_pileup_cnv.csv', 'w', newline='') as f:
#		w = csv.writer(f, delimiter='\t')
#		samples = [inputfiles[i].split('/')[-1] for i in range(len(inputfiles))]
#		sample_names = [samples[i].split('_')[0] for i in range(len(samples))]
#		print(sample_names)
#		w.writerow(['chr', 'start', 'end', *sample_names])
#		
#		for ckey in coverageAll:
#			print(ckey)
#			# error: for bkey in coverageAll[ckey]: iterator, string, why?
#			for bkey in coverageAll[ckey].keys():
#				if len(coverageAll[ckey][bkey]) == len(samfiles):
#					print(bkey)
#					print(coverageAll[ckey][bkey])
#					w.writerow([ckey, bkey*binSize, (bkey+1)*binSize,
#					*[round(coverageAll[ckey][bkey][idx],2) for idx in sorted (coverageAll[ckey][bkey].keys())]
#                                        ])
				# ensure that order is correct
def outputToFile(ckey):
	with open(outDir + datetime.now().strftime('%I:%M%p_Nov%d') + '_pileup_cnv.%s' % (ckey), 'w', newline='') as f:
		w = csv.writer(f, delimiter='\t')
		
		if ckey == 'chr1':
			samples = [inputfiles[i].split('/')[-1] for i in range(len(inputfiles))]
			sample_names = [samples[i].split('_')[0] for i in range(len(samples))]
	#		print(sample_names)
			w.writerow(['chr', 'start', 'end', *sample_names])
		
		#for ckey in coverageAll:
		#	print(ckey)
			# error: for bkey in coverageAll[ckey]: iterator, string, why?
		for bkey in coverageAll[ckey].keys():
			if len(coverageAll[ckey][bkey]) == len(samfiles):
				#print(bkey)
				#print(coverageAll[ckey][bkey])
				w.writerow([ckey, bkey*binSize, (bkey+1)*binSize,
					*[round(coverageAll[ckey][bkey][idx],2) for idx in sorted (coverageAll[ckey][bkey].keys())]])
	

def initOutputDir():
	try:
		os.makedirs(outDir)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

inputfiles = []
def main():
	initOutputDir()
	control_samples = [
		'../samples/down_syndrome/JLKD062_filtered_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3536978_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3536980_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3537015_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3537016_trimmed_bismark_bt2.deduplicated.sorted.bam'
	]
	ds_samples = [
		'../samples/down_syndrome/SRR3537005_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3537006_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3537007_trimmed_bismark_bt2.deduplicated.sorted.bam',
		'../samples/down_syndrome/SRR3537008_trimmed_bismark_bt2.deduplicated.sorted.bam'	
	]		
	inputfiles.extend(control_samples + ds_samples) 
	#inputfiles.extend([control_samples[4], ds_samples[3]])
	print(inputfiles)
	print(len(inputfiles))
	for i in range(len(inputfiles)):
		samfiles.append(pysam.AlignmentFile(inputfiles[i], 'rb'))
		print('\nsamfiles[%s] = %s' % (i, inputfiles[i]))
	
	initPositionsAll()
	
	for sampleIdx in range(len(samfiles)):
		poolBounds = mp.Pool(processes = 25)
		poolBounds.starmap(findChrmBoundaries, zip(repeat(sampleIdx),chrmOfInterest))
		poolBounds.close()
		poolBounds.join()
	print('end of poolBounds = %s' % (datetime.now() - startTime))
	
	pool = mp.Pool(processes = 25)
	pool.map(pileupCountsLoop, range(len(samfiles)))
#	pool.map(pileupCountsLoopSamples, chrmOfInterest)
# 	gives pileup read errors
	pool.close()
	pool.join()
	print('end of pool = %s' % (datetime.now() - startTime))
	
	startFileTime = datetime.now()	
	print('startOutputToFile = %s' % (startFileTime))
	#outputToFile()
	poolFile = mp.Pool(processes = 25)
	poolFile.map(outputToFile, chrmOfInterest)
	poolFile.close()
	poolFile.join()
	print('endOutputToFile = %s' % (datetime.now() - startFileTime))

# to check for EOF:
# https://www.biostars.org/p/247903/

# dynamically create list and name it in python
# https://stackoverflow.com/questions/1373164/how-do-i-create-a-variable-number-of-variables
# https://stackoverflow.com/questions/31882211/how-to-create-list-with-index-in-name-in-python

# https://stackoverflow.com/questions/22487296/multiprocessing-in-python-sharing-large-object-e-g-pandas-dataframe-between
# share pandas dataframe between processes: create data_handler child process, or below:
# https://stackoverflow.com/questions/10721915/shared-memory-objects-in-multiprocessing/58514648#58514648

if __name__ == "__main__":
# main 
	startTime = datetime.now()
	print('startTime = %s' % (startTime))
	main()
	for i in range(len(samfiles)):
	       	samfiles[i].close()
	print('endTime = %s' % (datetime.now() - startTime))

# 4421M VIRT, 1484M RES, CPU ~82%, TIME 19h+
# 4421M *2 or 3
# sorted merged output 
# ls *.chr* -1v | xargs cat > cnv_caller_output_Nov08.txt		
