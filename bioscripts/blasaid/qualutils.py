"""
Functions for manipulating quality scores on sequences.

We can usually expect the quality scores to be stored in `letter_annotations`,
in a dictionary under the key `phred_quality`. The utilities here should give
us the ability to cope with other quality scores types and missing scores.

"""

### IMPORTS

### CONSTANTS & DEFINES

DEFAULT_QUAL_KEY = 'phred_quality'


### IMPLEMENTATION ###

def get_qual_key (seqrec):
	if seqrec.letter_annotations.has_key (DEFAULT_QUAL_KEY):
		return DEFAULT_QUAL_KEY
	else:
		keys = seqrec.letter_annotations.keys()
		if keys:
			return keys[0]
		else:
			return None
	
	
def get_qual_arr (seqrec, key=None):
	key = key or DEFAULT_QUAL_KEY
	return seqrec.letter_annotations.get (key, None)


def set_qual_arr (seqrec, arr, key=None):
	key = key or DEFAULT_QUAL_KEY
	seqrec.letter_annotations[key] = arr


def add_dummy_quals (seqrec):
	qual_arr = [93] * len (seqrec)
	set_qual_arr (seqrec, qual_arr, key=DEFAULT_QUAL_KEY)
	

### END #######################################################################
