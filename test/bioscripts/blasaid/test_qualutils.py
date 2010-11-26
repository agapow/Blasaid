"""
Test fixtures for bioscripts.blasaid.qualutils.
"""

### IMPORTS

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from bioscripts.blasaid.qualutils import *


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

def test_get_qual_key():
	assert get_qual_key (ANNOTATED_SEQ) == 'phred_quality'
	assert get_qual_key (UNANNOTATED_SEQ) == None
	assert get_qual_key (ANNOTATED_SEQ_WITH_OTHER_KEY) == 'quality'

	
def test_get_qual_arr():
	assert get_qual_arr (ANNOTATED_SEQ) == QUAL_ARR
	assert get_qual_arr (UNANNOTATED_SEQ) == None
	assert get_qual_arr (ANNOTATED_SEQ_WITH_OTHER_KEY, 'quality') == QUAL_ARR


def test_set_qual_arr():
	new_srec = SeqRecord (Seq ('ACGT'),
		letter_annotations={'phred_quality': QUAL_ARR})
	set_qual_arr (new_srec, SUB_QUAL_ARR)
	assert get_qual_arr (new_srec) == SUB_QUAL_ARR
	new_srec = SeqRecord (Seq ('ACGT'))
	set_qual_arr (new_srec, SUB_QUAL_ARR)
	assert get_qual_arr (new_srec) == SUB_QUAL_ARR
	

def test_add_dummy_quals():
	new_srec = SeqRecord (Seq ('ACGT'))
	assert get_qual_key (new_srec) == None
	assert get_qual_arr (new_srec) == None
	add_dummy_quals (new_srec)
	assert get_qual_key (new_srec) == 'phred_quality'
	assert get_qual_arr (new_srec)

	
	
### END #######################################################################
