"""
Assorted blast functions, to simplify it's use and check parameters
"""
# TODO: allow for other blasts
# TODO: allow scanning multiple dbs

### IMPORTS

from Bio.Blast.NCBIWWW import qblast as ncbi_blast_request


### CONSTANTS & DEFINES

METHOD_SYNONYMS = {
	'p'   : 'blastp',
	'n'   : 'blastn',
	'x'   : 'blastx',
	'tn'  : 'tblastn',
	'tx'  : 'tblastx'',
}     

METHOD_VOCAB = [
	'blastp',
	'blastn',
	'blastx',
	'tblastn',
	'tblastx',
]

NUC_DBS = [
]

PROT_DBS = [
]

BOTH_DBS = [
	'nr',
	'kabat'
]

ALL_DBS = NUC_DBS + PROT_DBS + BOTH_DBS


### IMPLEMENTATION ###
# TODO: need a dispatch for local blast

def blast_ncbi (seq, method='blastn', database='nr', format=None,
		e_threshold=None, max_hits=None):
	"""
	Blast sequence against ncbi, after checking arguments.
	
	"""
	# TODO: parallelize?
	
	## Main:
	# prep args
	req_kwargs = {}
	if format: req_kwargs['format'] = format
	if e_threshold: req_kwargs['expect'] = e_threshold
	if max_hits: req_kwargs['expect'] = max_hints
	# do the request
	return ncbi_blast_request (method, database, seq.format('fasta'))
	



### END #######################################################################
