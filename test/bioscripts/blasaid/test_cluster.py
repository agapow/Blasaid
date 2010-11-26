"""
Test fixtures for bioscripts.blasaid.qualutils.
"""

### IMPORTS

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from bioscripts.blasaid.cluster import *


### CONSTANTS & DEFINES

SHORT_SEQ = SeqRecord (Seq ('AACCGGTT'))
LONG_SEQ = SeqRecord (Seq ('CCGGTTTTTTAA'))
DUP_SEQ = SeqRecord (Seq ('AACCGGTT'))
SUB_SEQ = SeqRecord (Seq ('CCGGTT'))
LONG_SEQ2 = SeqRecord (Seq ('AAAATTTTTTAA'))


### TESTS ###

def test_identity():
	assert identity (SHORT_SEQ, DUP_SEQ)
	assert identity (SHORT_SEQ, LONG_SEQ) == None
	
	
def test_subsequence():
	assert (subsequence (SHORT_SEQ, DUP_SEQ))
	assert (subsequence (SUB_SEQ, LONG_SEQ) == LONG_SEQ)
	assert subsequence (SHORT_SEQ, LONG_SEQ) == None
	

class test_similarity():
	clusterer = PairwiseSimilarity (.75)
	assert clusterer (SHORT_SEQ, DUP_SEQ)
	assert clusterer (SUB_SEQ, LONG_SEQ2) == None
	assert clusterer (SUB_SEQ, SHORT_SEQ) == SHORT_SEQ
	

### END #######################################################################
	