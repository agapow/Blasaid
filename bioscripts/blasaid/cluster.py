"""
Various clustering functions.

Cluster functors are callable, and accept two sequences as arguments. They
return None if there is no match. If there is a match, they return the
preferred / better / longer sequence.

This imposes obvious limitations on the clustering process (and makes it order
dependent) but it serves our purposes.
"""

### IMPORTS

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from bioscripts.blasaid.filter import *


### CONSTANTS & DEFINES

### IMPLEMENTATION ###

def identity (seq_1, seq_2):
	if seq_1.seq.tostring() == seq_2.seq.tostring():
		return seq_1
	return None
	
	
def subsequence (seq_1, seq_2):
	str_1 = seq_1.seq.tostring()
	str_2 = seq_2.seq.tostring()
	if str_1 in str_2:
		return seq_2
	if str_2 in str_1:
		return seq_1
	return None
	

class PairwiseSimilarity (object):
	"""
	Cluster sequences above a certain level of similarity, as judged by pairwise alignment.
	"""
	
	def __init__ (self, threshold):
		self.threshold = threshold
		
	def __call__ (self, seq1, seq2):
		# sort by length
		sorted_seqs = sorted ([seq1, seq2],
			cmp=lambda x, y: cmp (len (x.seq), len (y.seq)))
		min_seq_len = len (sorted_seqs[0])
		from Bio import pairwise2 as prwise
		res = prwise.align.globalxx (*[s.seq for s in sorted_seqs], score_only=True) / min_seq_len
		if (self.threshold <= res):
			return sorted_seqs[1]
		else:
			return None



### END ######################################################################
