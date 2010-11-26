"""
Test fixtures for bioscripts.blasaid.qualutils.
"""

### IMPORTS

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from bioscripts.blasaid.filter import *


### CONSTANTS & DEFINES

QUAL_ARR = [10, 20, 30, 40]
SUB_QUAL_ARR = [1, 2, 3, 4]

ANNOTATED_SEQ = SeqRecord (
	Seq ('ACGT'),
	letter_annotations={'phred_quality': QUAL_ARR},
)
UNANNOTATED_SEQ = SeqRecord (Seq ('ACGT'))
ANNOTATED_SEQ_WITH_OTHER_KEY = SeqRecord (
	Seq ('ACGT'),
	letter_annotations={'quality': QUAL_ARR},
)


### TESTS ###

def test_minlength():
	filter = MinLength(4)
	seq_len3 = SeqRecord (Seq ('ACG'))
	assert filter (seq_len3) == False
	seq_len4 = SeqRecord (Seq ('ACGT'))
	assert filter (seq_len4) == True
	seq_len5 = SeqRecord (Seq ('ACGTA'))
	assert filter (seq_len5) == True
	

def test_minbasequality():
	filter = MinBaseQuality(30)
	seq_len_a = SeqRecord (Seq ('ACGT'), letter_annotations={'phred_quality': [10, 20, 30, 40]})
	assert filter (seq_len_a) == False
	seq_len_b = SeqRecord (Seq ('ACGT'), letter_annotations={'phred_quality': [30, 30, 30, 40]})
	assert filter (seq_len_b) == True
	seq_len_c = SeqRecord (Seq ('ACGT'), letter_annotations={'phred_quality': [40, 40, 40, 40]})
	assert filter (seq_len_c) == True
	

def test_minbasequality():
	filter = MinAvgBaseQuality(30)
	seq_len_a = SeqRecord (Seq ('ACGT'), letter_annotations={'phred_quality': [10, 20, 30, 40]})
	assert filter (seq_len_a) == False
	seq_len_b = SeqRecord (Seq ('ACGT'), letter_annotations={'phred_quality': [30, 20, 30, 40]})
	assert filter (seq_len_b) == True
	seq_len_c = SeqRecord (Seq ('ACGT'), letter_annotations={'phred_quality': [10, 90, 40, 40]})
	assert filter (seq_len_c) == True


### END #######################################################################
	