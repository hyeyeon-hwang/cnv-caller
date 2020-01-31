#!usr/bin/env python3

import argparse
import sys
import os
import pysam
from datetime import datetime
from multiprocessing import Pool
import csv
import pandas as pd
from sklearn.cluster import KMeans

extended_help = """
Input the path to the directory containing the bam files.

To run:
python3 sex_caller.py --samples /share/lasallelab/Hyeyeon/projects/cnv-caller/allSamplesASD/ --cores 56 --sampleInfo /share/lasallelab/Hyeyeon/projects/cnv-caller/ASD_sample_info.csv
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

parser.add_argument(
	'--cores',
	required=False,
	type=int,
	default = 1,
	metavar='<int>',
	help='number of cores to use for parallel processing')

parser.add_argument(
	'--sampleInfo',
	required=False,
	type=str,
	default = None,
	metavar='<path>',
	help = "path to csv file that contains sex info for each sample")
	
arg = parser.parse_args()

# Write print statement outputs to file
# sys.stdout = open(datetime.now().strftime('%I:%M%p_%b%d') + '_sex_caller.out', 'w')

sexChrm = ["chrX", "chrY"]

def getBamfiles(path):
	bamfilePaths = []
	for bamfile in os.listdir(path):
		if bamfile.endswith('.bam'):
			bamfilePaths.append(path + '/' + bamfile)
	return (bamfilePaths)

def countSexChrmReads(bamfile):
	with open('tempfile.txt', 'a', newline='') as tempfile:
		tempfile = csv.writer(tempfile, delimeter='\t')
		stats{}
		
		samfile = pysam.AlignmentFile(bamfile, 'rb')
		stats[bamfile] = {'chrX': 0, 'chrY': 0}
		for read in samfile:
			if read.reference_name in sexChrm and read.is_duplicate == False:
				stats[bamfile][read.reference_name] += 1
		
		filename = bamfile.split("/")[-1]
		samplename = filename.split("_")[0]

		tempfile.writerow([
			samplename,
			stats[bamfile]['chrX'],
			stats[bamfile]['chrY'],
			stats[bamfile]['chrY']/stats[bamfile]['chrX'],
			])	
	
def predictSex():
	stats = pd.read_csv("tempfile.txt", sep='\t', engine='python')
	
	kmeans = KMeans(n_clusters = 2)
	kmeansData = stats[['ChrY:ChrX_ratio']] # add in placeholder column
	predictions = kmeans.fit_predict(kmeansData)
	centers = kmeans.cluster_centers_

	print(centers)
	print(predictions)

	if centers[0][1] > centers[1][1]:
		index0 = 'Male'
		index1 = 'Female'
	else:
		index0 = 'Female'
		index1 = 'Male'	

		
	
# split into predictSex
# output to file
# compute sex ratio
def sexCaller(bamfile):
	print(bamfile)
	with open("sex_caller_output_asd_Jan28.txt", 'a', newline='') as outfile:
		outfile = csv.writer(outfile, delimiter='\t')
		stats = {}
		print('bamfile = %s' % bamfile)
		tstart = datetime.now()
		
		samfile = pysam.AlignmentFile(bamfile, 'rb') 
		stats[bamfile] = {"chrX": 0, "chrY":0, "sex": None}		
		for read in samfile:
			if read.reference_name in sexChrm and read.is_duplicate == False:
				stats[bamfile][read.reference_name] += 1
		sexChrmRatio = stats[bamfile]["chrY"]/stats[bamfile]["chrX"]
			
		
		
		if stats[bamfile]["chrY"]/stats[bamfile]["chrX"] > 0.06: 
		# ratio for humans, different for cffDNA, make cutoff value an argument
			stats[bamfile]["sex"] = "Male"
		else: 
			stats[bamfile]["sex"] = "Female"
		
		
		
		filename = bamfile.split("/")[-1]
		samplename = filename.split("_")[0]

		if arg.sampleInfo == None:
	
			outfile.writerow([
				samplename,
				stats[bamfile]["chrX"],
				stats[bamfile]["chrY"],
				sexChrmRatio,
				round(sexChrmRatio * 100, 5),
				stats[bamfile]["sex"]
			])	
		else:
			sampleInfo = getSampleInfo(arg.sampleInfo)
			print(sampleInfo.loc[samplename, 'Sex'])

			outfile.writerow([
				samplename,
				stats[bamfile]["chrX"],
				stats[bamfile]["chrY"],
				sexChrmRatio,
				round(sexChrmRatio * 100, 5),
				stats[bamfile]["sex"],
				sampleInfo.loc[samplename, 'Sex']
			])	
			
		
		print('time for bamfile = %s' % (datetime.now() - tstart))

# sort output file by Name?
# def sortOutput(outputFile):

# infoFile = "ASD_sample_info.csv"
def getSampleInfo(sampleInfoFile):
	sampleInfo = pd.read_csv(sampleInfoFile)

	if 'Sex' in sampleInfo.columns and 'Name' in sampleInfo.columns:
		sampleInfo.index = list(sampleInfo['Name']) # set row names as sample names
		return sampleInfo			
	else:
		print("Sample info csv file must contain the following columns: 'Name', 'Sex'")
		# read first line of csv - should contain header
		# required columns are "Name" and "Sex"
			  	
	return sampleInfo
	
if __name__ == '__main__':
	timestart = datetime.now()	
	print('timestart = %s' % timestart)

	# "../../cffDNA/bam/"
	# arg.samples = /share/lasallelab/Hyeyeon/projects/cnv-caller/allSamplesASD
	bamfilePaths = getBamfiles(arg.samples)
	end_1 = datetime.now()
	print(*bamfilePaths, sep='\n')
	print(len(bamfilePaths))
	print('getBamfiles() %s' % (end_1 - timestart))

#	with open("sex_caller_output_asd_Jan28.txt", 'w', newline='') as outfile:
#		outfile = csv.writer(outfile, delimiter='\t')
#		if arg.sampleInfo == None:
#			outfile.writerow([	
#				'Sample_name', 'ChrX_reads', 'ChrY_reads', \
#				'ChrY:ChrX_ratio', 'ChrY:ChrX_percent', 'Sex'])		
#		else:
#			outfile.writerow([	
#				'Sample_name', 'ChrX_reads', 'ChrY_reads', \
#				'ChrY:ChrX_ratio', 'ChrY:ChrX_percent', 'Sex', 'Sample_info_sex'])		
	with open('tempfile.txt', 'w', newline='') as tempfile:
		tempfile = csv.writer(tempfile, delimeter='\t')
		tempfile.writerow(['Sample_name', 'ChrX_reads', 'ChrY_reads', 'ChrY:ChrX_ratio'])	

	
	# arg.cores = 56
	pool = Pool(processes = arg.cores) 
	pool.map(sexCaller, bamfilePaths)
	pool.close()
	pool.join()
	
	print('timeend = %s' % (datetime.now() - timestart))
