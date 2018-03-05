import itertools
import numpy as np
from sklearn.preprocessing import normalize

refFile = "data/reference.txt"
bases = ['A', 'T', 'C', 'G']
rcDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def kmerToIndex(string, k):
	assert(len(string) == k)
	index = 0
	for i in range(k):
		index += (4**i)*bases.index(string[k-1-i])
	return index

def indexToKmer(index, k):
	kmer = ''
	for i in range(k):
		kmer = bases[index%4] + kmer
		index //=4
	return kmer

def testConversion(k):
	k = 6
	kmerList = [''.join(p) for p in itertools.product(bases, repeat=k)]
	for i in range(len(kmerList)):
		assert(kmerToIndex(kmerList[i], k) == i)
		assert(indexToKmer(i, k) == kmerList[i])

def getRC(string):
	return ''.join([rcDict[b] for b in string])[::-1]

def parseRefData(k):
	output = list()
	kmerList = [''.join(p) for p in itertools.product(bases, repeat=k)]
	x = set()
	for kmer in kmerList:
		x.add(min(kmerToIndex(kmer, k), kmerToIndex(getRC(kmer), k)))
	x = list(x)
	with open(refFile) as f:
		for line in f:
			if line.startswith('>'):
				enhanced = line.strip()
				seq = ''
				continue
			if len(line.strip()) > 0:
				seq += line.strip().upper()
			elif len(seq) > 0:
				kmerVec = np.zeros(4**k)
				for i in range(len(seq) - k):
					subseq = seq[i:i + k]
					if set(subseq).issubset(set(bases)):
						kmerVec[min(kmerToIndex(subseq, k), kmerToIndex(getRC(subseq), k))] += 1
				output.append([kmerVec[x], enhanced])
				seq = ''
	if len(seq) > 0:
		kmerVec = np.zeros(4**k)
		for i in range(len(seq) - k):
			subseq = seq[i:i + k]
			if set(subseq).issubset(set(bases)):
				kmerVec[min(kmerToIndex(subseq, k), kmerToIndex(getRC(subseq), k))] += 1
		output.append([kmerVec[x], enhanced])

	return output

def processSVMData(data, region):
	outputX, outputY = zip(*data)
	outputX = normalize(outputX, norm='l2')
	outputY = np.array([int(region in y) for y in outputY])
	return [outputX, outputY]