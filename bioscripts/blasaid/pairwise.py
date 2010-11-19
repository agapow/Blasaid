"""
Assess pairwise similarity using alignment.
"""

### IMPORTS


### CONSTANTS & DEFINES

# methods that directly return a similarity score
DIRECT_METHODS = [
	'identity',
]

# methods return return an alignment 
ALIGN_METHODS = [
	'pairwise',
]

ALL_METHODS = DIRECT_METHODS + ALIGN_METHODS


### IMPLEMENTATION ###

def similarity (meth, seq1, seq2, opts={}):
	## Preconditions:
	assert (meth in ALL_METHODS), "unrecognised similarity method '%s'" % meth
	if meth is None:
		meth = 'identity'
	## Main:
	if meth in DIRECT_METHODS:
		if meth == 'identity':
			sim = compare_identity (seq1, seq2)
		return sim
	elif meth in ALIGN_METHODS:
		if meth == 'identity':
			aln = align_pairwise (seq1, seq2)
		# do something with aln
	else:
		assert False, "unrecognised similarity method '%s'" % meth
		

def compare_identity (seq1, seq2):
	if seq1.data == seq2.data:
		return 1.0
	else:
		return 0.0
	

def align_pairwise (seq1, seq2):
	try:
		from Bio import cpairwise2 as prwise
	except:
		from Bio import pairwise2 as prwise
	res = prwise.align.globalxx (seq1.data, seq2.data)
	return max ([x[2]/x[4] for x in res])
	



### END #######################################################################
