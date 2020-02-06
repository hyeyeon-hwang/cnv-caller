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
python3 sex_caller.py --samples /share/lasallelab/Hyeyeon/projects/cnv-caller/allSamplesASD/ --cores 27 --sampleInfo /share/lasallelab/Hyeyeon/projects/cnv-caller/ASD_sample_info.csv
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

def getBamfiles(path):
	bamfilePaths = []
	for bamfile in os.listdir(path):
		if bamfile.endswith('.bam'):
			bamfilePaths.append(path + '/' + bamfile)
	return (bamfilePaths)

def countSexChrmReads(bamfile):
	with open('tempfile.txt', 'a', newline='') as tempfile:
		tempfile = csv.writer(tempfile, delimiter='\t')
		stats = {}
		samfile = pysam.AlignmentFile(bamfile, 'rb')
		stats[bamfile] = {'chrX': 0, 'chrY': 0}
		for read in samfile:
			if read.reference_name in ['chrX', 'chrY'] and read.is_duplicate == False:
				stats[bamfile][read.reference_name] += 1
		
		filename = bamfile.split("/")[-1]
		samplename = filename.split("_")[0]
			
		sampleInfo = getSampleInfo(arg.sampleInfo)
			
		# started 2.51 G, max at 2.74 G
		sexChrmRatio = stats[bamfile]['chrY']/stats[bamfile]['chrX']
		
		tempfile.writerow([
			samplename,
			stats[bamfile]['chrX'],
			stats[bamfile]['chrY'],
			sexChrmRatio,
			sexChrmRatio * 100,
			sampleInfo.loc[samplename, 'Sex']
			])	

def predictSex():
	stats = pd.read_csv("tempfile.txt", sep='\t', engine='python')
	# Delete tempfile after reading in
	# os.remove('tempfile.txt')

	# round here because odd rounding behavior occurs when reading in rounded values
	stats['ChrY:ChrX_percent'] = round(stats['ChrY:ChrX_percent'], 5)
	# stats['ChrY:ChrX_percent'] = float('{0:.5f}'.format(stats['ChrY:ChrX_percent']))
		
	kmeansData = stats[['ChrY:ChrX_ratio']].copy()
	# Add placeholder column because kmeans.fit_predict() requires 2D data
	kmeansData['placeholder'] = [0] * len(kmeansData.index)
	kmeans = KMeans(n_clusters = 2)
	predictions = kmeans.fit_predict(kmeansData)
	centers = kmeans.cluster_centers_

	# divide two center values
	# choose smaller or larger as the numerator or denominator
	# if close to 1, then all samples are of one sex
	# if greater than 0.5? most will be greater than 0.9, very close to 1
	# 0.004/0.12 = 0.03, not even 0.1 
		
	index0 = None
	index1 = None
	if centers[0][0] > centers[1][0]:
		index0 = 'M'
		index1 = 'F'
	else:
		index0 = 'F'
		index1 = 'M'	

	predictedSex = []
	for predIdx in predictions:
		if predIdx == 0:
			predictedSex.append(index0)
		else:
			predictedSex.append(index1)
	
	stats['Predicted_sex'] = predictedSex

	# Add 'Sex_mismatch' column if there were mismatches in the 'Sample_info_sex' and 'Predicted_sex' columns
	# Rows with mismatches are labeled as 'Mismatch'
	# Rows with consistent sex labels are marked with '.'
	if stats['Sample_info_sex'].equals(stats['Predicted_sex']) == False:
		sexMismatchList = ['.'] * len(stats.index)
		for i in range(0, len(stats.index)):
			if stats['Sample_info_sex'][i] != stats['Predicted_sex'][i]:
				sexMismatchList[i] = "Mismatch"
	stats['Sex_mismatch'] = sexMismatchList
	stats = stats.sort_values(by = 'Sample_name')
	stats.to_csv(datetime.now().strftime('%I:%M%p_%b%d') + '_sex_caller_output.txt', sep = '\t', index = False)

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

	with open('tempfile.txt', 'w', newline='') as tempfile:
		tempfile = csv.writer(tempfile, delimiter='\t')
		tempfile.writerow([
			'Sample_name', 'ChrX_reads', 'ChrY_reads', 
			'ChrY:ChrX_ratio', 'ChrY:ChrX_percent',
			'Sample_info_sex'])

	
	pool = Pool(processes = arg.cores) 
	pool.map(countSexChrmReads, bamfilePaths)
	pool.close()
	pool.join()
	
	predictSex()
	print('timeend = %s' % (datetime.now() - timestart))
