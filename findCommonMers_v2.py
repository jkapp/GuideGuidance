#!/usr/bin/env python3
#Written: Josh Kapp (jkapp@ucsc.edu, joshkapp01@gmail.com)
#argparse usage modified from code written by David Bernick
###########
import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import argparse
###########
class CommandLine():
	def __init__(self, inOpts=None):
		self.parser = argparse.ArgumentParser()
		self.parser.add_argument('-r', '--r', action = 'store', help='reference file name')
		self.parser.add_argument('-t', '--t', action = 'store', help='mer count threshold')
		if inOpts is None :
			self.args = self.parser.parse_args()
		else :
			self.args = self.parser.parse_args(inOpts)
###########
class FindCommonMers:
	def __init__(self, inRef, inThresh):
		self.refFile = inRef #Reference genome in fasta format
		self.merThresh = int(inThresh) #User defined mer count threshold

		self.refSeq = '' #Load record into memory for fast parsing
		self.merDic = {} #Dictionary of 17mers
		self.outFile = '' #Name of output file

	def loadRef(self): 
		recordLen = 0  #how long is the current record
		currentPos = 0 #where are we as we inch down the current record
		currentMer = '' #whats the current 17mer
		merVal = 0 #How many instances of the current 17mer have we seen

		for record in SeqIO.parse(self.refFile, 'fasta'): ###Parse reference genome and throw into memory for now, one index for each header
			self.refSeq = str(record.seq) #load current record 
			recordLen = len(self.refSeq) #length of current record
			for inchFor in range(0,recordLen - 16): #start at index pos 0 and go all the down the current record, one base at a time until you hit record-16 
				currentMer = str(self.refSeq[currentPos:currentPos + 17]) #get current mer, 17 bases total
				currentMer = currentMer.upper() #uppercase standard incase theres a mixture in the record
				currentPos += 1 #increase current position by 1 for next interation
				if 'N' in currentMer: #disregard mers with any amount of Ns, mer would inflate dic with uninformative data
					continue
				if len(currentMer) < 17: #shouldnt happen but disregard mers shorter than 17
					continue
				if currentMer not in self.merDic: #add mer to list if its unique thus far
					self.merDic[currentMer] = 1 
				else:
					merVal = self.merDic[currentMer] #if already in the dic then count value up by 1
					self.merDic[currentMer] = merVal + 1
				currentMer = '' #reset mer for next iteration
				merVal = 0 #reset mer Val 

			currentPos = 0
			self.refSeq = Seq(self.refSeq, generic_dna)
			self.refSeq = self.refSeq.reverse_complement() #reverse comp the seq
			self.refSeq = str(self.refSeq) #make the seq an object again
			for inchRev in range(0, recordLen - 16): #do it all again but REVERSED
				currentMer = str(self.refSeq[currentPos:currentPos + 17]) #get current mer, 17 bases total
				currentMer = currentMer.upper() #uppercase standard incase theres a mixture in the record
				currentPos += 1 #increase current position by 1 for next interation
				if 'N' in currentMer: #disregard mers with any amount of Ns, mer would inflate dic with uninformative data
					continue
				if len(currentMer) < 17: #shouldnt happen but disregard mers shorter than 17
					continue
				if currentMer not in self.merDic: #add mer to list if its unique thus far
					self.merDic[currentMer] = 1 
				else:
					merVal = self.merDic[currentMer] #if already in the dic then count value up by 1
					self.merDic[currentMer] = merVal + 1
				currentMer = '' #reset mer for next iteration
				merVal = 0 #reset mer Val 

			recordLen = 0 #reset record Len 
			self.refSeq = '' #reset record


	def printMerOut(self): #Instead of writing to the command line this would populate a set() for fast parsing
		for key in self.merDic:
			if self.merDic[key] >= self.merThresh:
				print(key)
		print('Total Mers:', len(self.merDic))
###########
def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.  

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(myCommandLine)

    print (myCommandLine.args) # print the parsed argument string .. as there is nothing better to do

    inRefparam = myCommandLine.args.r
    inThreshparam = myCommandLine.args.t
    paramList = [inRefparam, inThreshparam]

    return paramList
###########
collectingParams = main()
#print(collectingParams)
inRef = collectingParams[0]
inThresh = collectingParams[1]
x = FindCommonMers(inRef, inThresh)
x.loadRef()
x.printMerOut()
