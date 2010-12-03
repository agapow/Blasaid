"""
Script for automated blasting of NGS data.

Aim
---

While next-generation sequencing (NGS) provides a bounty of data for analysis,
identifying sequence fragments without a reference genome is arduous. Each
fragment must be blasted seperately against a reference database, and the run
will include duplicates, largely overlapping sequences or those that are too
short or of too low quality. **blasaid** is a tool to helps reduce sequencing
runs to the more important or interesting reads, blasts each automatically and
synosizes the results for easy reading.

Function
--------

Each *blasaid* run consists of the following steps:

1. **Merge** all input sequences, incorporating any accompanying quality files.
Optionally **trim** 3' end of sequence to rid low quality trailing segments.

2. Optionally **filter** the input sequences, getting rid of any that fall below
thresholds in quality or length.

3. Optionally **cluster** the input sequences, reducing similar reads to a
single example

4. **Blast** the remaining reads against a chosen database and return results in
a given format.

5. The results are optionally **reduced** to the most interesting and
**summarized**. 

"""
# TODO: a blast summariser / reducer script

__docformat__ = 'restructuredtext en'
__author__ = 'Paul-Michael Agapow <paul-michael.agapow@hpa.org.uk>'


### IMPORTS ###

from os import path
import sys
from exceptions import BaseException, SystemExit
import traceback
from datetime import datetime
from itertools import islice

from Bio import AlignIO, SeqIO
from Bio.Blast import NCBIXML

from bioscripts.blasaid import filter, scriptlog, cluster
from bioscripts.blasaid.scriptlog import log
from bioscripts.blasaid.exseqreader import ExSeqReader
from bioscripts.blasaid.blast import blast_ncbi
from bioscripts.blasaid.summarywriter import SummaryWriter

from bioscripts.blasaid import fileutils, common


### CONSTANTS & DEFINES ###

try:
	from bioscripts.blasaid import __version__ as VERSION
except:
	VERSION = 'unknown'
	
INPUT_FORMATS = [
	'auto',
	'fasta',
]
DEFAULT_INPUT_FORMAT = INPUT_FORMATS[0]

OUTPUT_FORMATS = [
	'text',
	'csv',
	'xml'
]
DEFAULT_OUTPUT_FORMAT = OUTPUT_FORMATS[0]

# whaere and how to store intermediate sequences in
SCRATCH_DIR_PREFIX = datetime.now().strftime ("blasaid-scratch-%Y%m%dT%H%M-")
SCRATCH_FORMAT = 'fastq'
SCRATCH_EXT = '.fastq'

RESULTS_DIR_PREFIX = datetime.now().strftime ("blasaid-results-%Y%m%dT%H%M-")

BLAST_DIR_NAME = 'blast-results'

_DEV_MODE = True


### IMPLEMENTATION ###
	
### FILE UTILS

def create_intermediate_file (base, scratch_dir):
	#inter_path = path.join (scratch_dir, base + SCRATCH_EXT)
	#return open (inter_path, 'w')
	return fileutils.open_for_writing ([scratch_dir, base + SCRATCH_EXT])


### BLAST OPERATIONS

def make_blast_raw_dir (scratch_dir):
	"""
	Make a directory for the saving of Blast result files.
	"""
	res_path = path.join (scratch_dir, BLAST_DIR_NAME)
	log.debug ("Making raw blast result directory at '%s' ..." % res_path)
	mkdir (res_path)
	return res_path


def blast_seqs (in_file, scratch_dir, result_dir, e_threshold=None,
		max_hits=None):
	log.info ('Blasting reads ...')
	#
	# read in seqfile
	blast_cnt = 0
	rdr = ExSeqReader (in_file, fmt=SCRATCH_FORMAT, merge_quals=False)
	# for each sequence
	for seq in rdr.read():
		log.info ("Blasting sequence '%s' ..." % seq.id)
		# blast it and write down to scratch
		res = blast_ncbi (seq.format('fasta'),
			e_threshold=e_threshold,
			max_hits=max_hits,
		)
		blast_cnt += 1
		scratch_hndl = fileutils.write_to_file ([result_dir, '%s.xml' % seq.id],
			res.read())
		# reduce and write to results
		res.seek(0)
		result_hndl = fileutils.write_to_file ([result_dir, '%s.xml' % seq.id],
			res.read())
		# summarize
		
	## Postconditions & return:
	log.info ('%s sequences were blasted ...' % blast_cnt)
		
		
def dump_options (opts):
	log.debug ('Dumping options ...')
	print ("Options:")
	for attr in [x for x in dir (opts) if x not in ['read_file', 'read_module', 'ensure_value'] and not x.startswith ('_')]:
		print ("- %s: %s" % (attr, getattr (opts, attr)))


def merge_and_trim_seqs (input_files, scratch_dir, input_format=None,
		merge_quals=True, trim_right=None):
	"""
	Combine input sequence files, with merging with quality data and trimming.
	
	These reads and merges all the inputs into a single file in the designated
	intermediate format, combining them with any quality data available and
	trimming the 3' end if required. It is possible that this stage may be
	unnecessary (if the data is in a single file in the right format already and
	no trimming is required) but this will almost never be the case.
	
	"""
	# TOOD: check for rare conditions where we don't need to do this
	
	## Main:
	if trim_right:
		log.info ('Merging & trimming sequences ...')
	else:
		log.info ('Merging sequences ...')
		
	# create & open file for filtered seqs
	merged_out_hndl = create_intermediate_file ('merged_and_trimmed', scratch_dir)
	
	# read in seqfiles
	seq_cnt = 0
	for f in input_files:
		log.info ("Reading '%s' ..." % f)
		# check file exists
		assert (path.exists (f)), "the sequence file '%s' does not exist" % f
		# set default format
		if (input_format in [None, 'auto']):
			fmt = None
		else:
			fmt = input_format
		# make reader & read
		rdr = ExSeqReader (f, fmt=fmt, merge_quals=merge_quals)
		for seq in rdr.read():
			log.debug ("Reading sequence '%s' ..." % seq.id)
			# trim the sequnece if requested
			if trim_right:
				seq.seq = seq.seq[:-trim_right]
			# write it out
			SeqIO.write ([seq], merged_out_hndl, SCRATCH_FORMAT)
			seq_cnt += 1
	merged_out_hndl.close()
	## Postconditions & return:
	log.info ('%s sequences merged ...' % seq_cnt)
	return merged_out_hndl.name


def make_clusterer (identity=None, subsequence=None, similarity=None):
	log.debug ('Checking clusterer ...')
	if identity:
		return cluster.identity
	if subsequence:
		return cluster.subsequence
	if similarity:
		return PairwiseSimilarity (similarity)
	return None
	

def make_filters (length=None, base_quality=None, avg_quality=None):
	log.debug ('Checking filters ...')
	filter_list = []
	if (length):
		filter_list.append (filter.MinLength (float(length)))
	if (base_quality):
		filter_list.append (filter.MinBaseQuality (float(base_quality)))
	if (avg_quality):
		filter_list.append (filter.MinAvgBaseQuality (float(avg_quality)))
	return filter_list


def filter_seqs (infile, scratch_dir, filters):
	## Main:
	log.info ('Filtering sequences ...')
	# create & open file for filtered seqs
	filtered_out_hndl = create_intermediate_file ('filtered', scratch_dir)
	# read in seqfile
	filter_cnt = 0
	log.info ("Reading '%s' ..." % infile)
	# make reader & read
	rdr = ExSeqReader (infile, fmt=SCRATCH_FORMAT, merge_quals=False)
	for seq in rdr.read():
		log.debug ("Reading sequence '%s' ..." % seq.id)
		# if it passes all filters
		if all (filters):
			log.debug ("Accepting '%s' ..." % seq.id)
			SeqIO.write ([seq], filtered_out_hndl, SCRATCH_FORMAT)
			filter_cnt += 1
		else:
			log.debug ("Rejecting '%s' ..." % seq.id)
	filtered_out_hndl.close()
	## Postconditions & return:
	log.debug ('%s sequences remain after filtering ...' % filter_cnt)
	return filtered_out_hndl.name


def cluster_seqs (infile, scratch_dir, cluster_fn):
	# NOTE: the logic of this is quite hairy. We wish to traverse the number of
	# times through the collection looking for similar seqs. At the same time,
	# we don't wish to load the entirity of sequences into memory. So we open the
	# infile file and read it one by one. For each sequence read, we read the file
	# again and compare it to every sequence after it. The comparison (clustering)
	# function returns None if the two sequences do not cluster. Otherwise it
	# returns the preferred / better sequence of the two, which is then used for
	# subsequent comparisons on this loop. We save on unnecessary comparsions by
	# storing any sucessful matches in a "don't check" dict. We use a single file
	# handle to searh one, rewinding as need be, to avoid opening and closing
	# thousands.
	
	## Main:
	log.info ('Clustering sequences ...')
	# create & open file for results, & open handl for searching / comparing
	clustered_out_hndl = create_intermediate_file ('clustered', scratch_dir)
	search_hndl = open (infile, 'r')
	# read in seqfile and check seqs one-by-one
	seq_cnt = 0
	already_tested = {}
	log.info ("Reading '%s' ..." % infile)
	rdr = ExSeqReader (infile, fmt=SCRATCH_FORMAT, merge_quals=False)
	for i, seq_1 in enumerate (rdr.read()):
		log.debug ("Reading sequence '%s' ..." % seq_1.id)
		# if this seq hasn't previously been clustered
		if seq_1.id not in already_tested:
			# move to start of search file & start reading
			search_hndl.seek (0)
			rdr_2 = ExSeqReader (search_hndl, fmt=SCRATCH_FORMAT, merge_quals=False)
			# for every seq beyond the current one
			for seq_2 in islice (rdr_2.read(), i+1, None):
				# if it hasn't previously been clustered
				if seq_2.id not in already_tested:
					cluster_seq = cluster_fn (seq_1, seq_2)
					# if it clusters, place in "done" dict and update search term
					if cluster_seq:
						already_tested[seq_2.id] = True
						seq_1 = cluster_seq
			# save the surviving search term
			SeqIO.write ([seq_1], clustered_out_hndl, SCRATCH_FORMAT)
			seq_cnt += 1
	clustered_out_hndl.close()
	search_hndl.close()
	## Postconditions & return:
	log.debug ('%s sequences remain after clustering ...' % seq_cnt)
	return clustered_out_hndl.name


def parse_args():	
	# Construct the option parser.
	optparser = common.OptionParser (
		prog          = 'blasaid',
		description   = 'Reduce and blast NGS data.',
		epilog        = 'See the README or program documentation for further details',
	)
	
	# TODO: filter by range
	
	optparser.add_option ('--input-format', '-i',
		dest="input_format",
		choices=INPUT_FORMATS,
		default=DEFAULT_INPUT_FORMAT,
		help='''The format of the input sequence files. If not supplied, this will be
			inferred from the extension of the files.''',
		metavar='SEQ_FORMAT',
	)
	
	optparser.add_option ('--ignore-qual-files',
		dest="ignore_qual_files",
		action='store_true',
		default=False,
		help='''Do not merge qual files into sequence files.'''
	)
	
	optparser.add_option ('--trim-right',
		dest="trim_right",
		type=int,
		default=None,
		help='''Trim the 3' end of sequences by this many bases''',
		metavar='LENGTH',
	)
	
	optparser.add_option ('--filter-by-length',
		dest="filter_length",
		type=int,
		default=None,
		help='''Sequences that are shorter than this are discarded''',
		metavar='LENGTH',
	)
	
	optparser.add_option ('--filter-by-base-quality',
		dest="filter_base_qual_threshold",
		type=float,
		default=None,
		help='''Sequences that contain any nucleotide with a quality lower than
			this threshold are discarded''',
		metavar='QUALITY',
	)

	optparser.add_option ('--filter-by-average-quality',
		dest="filter_avg_qual_threshold",
		type=float,
		default=None,
		help='''Sequences with an average quality lower than
			this threshold are discarded''',
		metavar='QUALITY',
	)
	
	optparser.add_option ('--cluster-by-identity',
		dest="cluster_identity",
		action='store_true',
		default=False,
		help='''Sequences that are identical are reduced to a single example''',
	)
	
	optparser.add_option ('--cluster-by-subsequence',
		dest="cluster_subsequence",
		action='store_true',
		default=False,
		help='''Sequences that are subsequences of others are eliminated''',
	)
	
	optparser.add_option ('--cluster-by-similarity',
		dest="cluster_similar",
		type=float,
		default=None,
		help='''Sequences that are more similar than this threshold are clustered
			into a single query''',
		metavar='PERCENT',
	)
	
	optparser.add_option ('--blast-max-e',
		dest="blast_expect",
		type=float,
		default=None,
		help='''Blast with this e-threshold (expectation) cut-off.''',
		metavar='PERCENT',
	)	

	optparser.add_option ('--blast-max-hits',
		dest="blast_max_hits",
		type=int,
		default=None,
		help='''Blast with this e-threshold (expectation) cut-off.''',
		metavar='PERCENT',
	)
	
	# TODO:
	#blast-max-gaps
	#blast-min-alignment
	#blast-min-score
	#blast-min-bit-score
	#blast-min-identity
	
	# TODO: call this results format?
	#optparser.add_option ('--output-format', '-o',
	#	dest="output_format",
	#	choices=OUTPUT_FORMATS,
	#	default=DEFAULT_OUTPUT_FORMAT,
	#	help='''The format of the output blast result files. If not supplied, this
	#		will default to %s.''',
	#	metavar='BLAST_FORMAT',
	#)
	
	optparser.add_results_path_option()
	optparser.add_scratch_path_option()
	optparser.add_debug_options()
	
	# parse 
	options, infiles = optparser.parse_args()
	
	## Postconditions & return:
	# can only use one clustering method
	assert (bool(options.cluster_identity) + bool(options.cluster_subsequence) +
		bool(options.cluster_similar) <= 1), "multiple clustering methods set"
	
	# return
	return infiles, options


def main():
	input_files, options = parse_args()

	try:
		# auxillary and help functions
		if options.dump_options:
			dump_options (options)
			
		# setup
		scriptlog.init_logger (
			format="%(message)s",
			verbosity=options.verbosity,
		)
		
		scratch_dir = fileutils.create_scratch_dir (SCRATCH_DIR_PREFIX,
			options.scratch_path)
		log.info ("Making temporary files at '%s' ..." % scratch_dir)
		
		# read in. merge and trim seq files
		work_file = merge_and_trim_seqs (input_files, scratch_dir,
			merge_quals=not options.ignore_qual_files,
			trim_right=options.trim_right,
			input_format=options.input_format
		)

		# filter stuff, if requested
		filters = make_filters (
			length = options.filter_length,
			base_quality = options.filter_base_qual_threshold,
			avg_quality = options.filter_avg_qual_threshold,
		)
		if (filters):
			work_file = filter_seqs (work_file, scratch_dir, filters)
		
		# cluster stuff, if requested
		clusterer = make_clusterer (
			options.cluster_identity,
			options.cluster_subsequence,
			options.cluster_similar,
		)
		if clusterer: 
			work_file = cluster_seqs (work_file, scratch_dir, clusterer)
			
		# blast stuff, if not a dryrun
		if not options.dryrun:
			# directory for raw results straight from server
			blast_out_dir = make_blast_raw_dir (scratch_dir)
			# directory for distilled results
			results_dir = fileutils.make_output_dir (RESULTS_DIR_PREFIX,
				options.results_path)
			log.info ("Making final results directory at '%s' ..." % results_dir)
		
			result_dir = blast_seqs (work_file,
				scratch_dir=blast_out_dir,
				result_dir=results_dir,
				e_threshold=options.blast_expect,
				max_hits=options.blast_max_hits,
			)
			# and combine, reduce & summarize if requested

	except BaseException, err:
		if options.traceback:
			traceback.print_exc()
		log.critical ("A problem: %s" % err)
		log.info ('Fatal error, terminating.')
		sys.exit(1)
		return
		
	log.info ('Finished.')
	sys.exit(0)
	

if __name__ == '__main__':
	main()


### END #######################################################################
