#!usr/bin/env python3

import sys
import os
import pysam
from datetime import datetime
from multiprocessing import Pool
import csv

# Write print statement outputs to file
sys.stdout = open(datetime.now().strftime('%I:%M%p_%b%d') + '_sex_caller.out', 'w')

sexChrm = ["chrX", "chrY"]

def getBamfiles(path):
	bamfilePaths = []
	for bamfile in os.listdir(path):
		if bamfile.endswith('.bam'):
			bamfilePaths.append(path + bamfile)
	return (bamfilePaths)

def sexCaller(bamfile):
	print(bamfile)
	with open("cffDNA_sex_output", 'a', newline='') as outfile:
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
	
if __name__ == '__main__':
	timestart = datetime.now()	
	print('timestart = %s' % timestart)

	bamfilePaths = getBamfiles("../../cffDNA/bam/")
	end_1 = datetime.now()
	print(*bamfilePaths, sep='\n')
	print(len(bamfilePaths))
	print('getBamfiles() %s' % (end_1 - timestart))

	with open("cffDNA_sex_output", 'w', newline='') as outfile:
		outfile = csv.writer(outfile, delimiter='\t')
		outfile.writerow(['Sample_name', 'ChrX_reads', 'ChrY_reads', 'ChrY:ChrX_ratio', 'ChrY:ChrX_%'])		
		
	pool = Pool(processes = 24)
	pool.map(sexCaller, bamfilePaths)
	pool.close()
	pool.join()
	
	print('timeend = %s' % (datetime.now() - timestart))
