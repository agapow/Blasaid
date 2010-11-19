"""
Script for automated blasting of NGS data.
"""

__docformat__ = 'restructuredtext en'
__author__ = 'Paul-Michael Agapow <paul-michael.agapow@hpa.org.uk>'


### IMPORTS ###

from os import path
import sys
from exceptions import BaseException, SystemExit
import traceback
import tempfile
from datetime import datetime

from Bio import AlignIO, SeqIO

from bioscripts.blasaid import filter, scriptlog
from bioscripts.blasaid.scriptlog import log
from bioscripts.blasaid.exseqreader import ExSeqReader


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

# what form to store intermediate sequences in
SCRATCH_FORMAT = 'fastq'
SCRATCH_EXT = '.fastq'

_DEV_MODE = True


### IMPLEMENTATION ###
	
def get_scratch_dir (path):
	"""
	Create a directory for scratch files.
	"""
	prefix = datetime.now().strftime ("blasaid-%Y%m%dT%H%M-")
	if (path):
		return tempfile.mkdtemp (prefix=prefix, dir=path)
	else:
		return tempfile.mkdtemp (prefix=prefix)
		
	
def dump_options (opts):
	log.debug ('Dumping options ...')
	print ("Options:")
	for attr in [x for x in dir (opts) if x not in ['read_file', 'read_module', 'ensure_value'] and not x.startswith ('_')]:
		print ("- %s: %s" % (attr, getattr (opts, attr)))


def merge_and_trim_seqs (input_files, scratch_dir, input_format=None,
		merge_quals=True, trim_right=options.trim_right):
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
		log.debug ('Merging & trimming sequences ...')
	else:
		log.debug ('Merging sequences ...')
		
	# create & open file for filtered seqs
	merged_out_path = path.join (scratch_dir, 'merged_and_trimmed' + SCRATCH_EXT)
	merged_out_hndl = open (merged_out_path, 'w')
	
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
	return merged_out_path


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
	log.debug ('Filtering sequences ...')
	# create & open file for filtered seqs
	filtered_out_path = path.join (scratch_dir, 'filtered.fasta')
	filtered_out_hndl = open (filtered_out_path, 'w')
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
			SeqIO.write ([seq], filtered_out_hndl, "fasta")
			filter_cnt += 1
		else:
			log.debug ("Rejecting '%s' ..." % seq.id)
	filtered_out_hndl.close()
	## Postconditions & return:
	log.debug ('%s sequences remain after filtering ...' % filter_cnt)
	return filtered_out_path


def cluster_seqs (infile, scratch_dir, cluster_fn):
	## Main:
	log.debug ('Filtering sequences ...')
	# create & open file for filtered seqs
	filtered_out_path = path.join (scratch_dir, 'filtered.fasta')
	filtered_out_hndl = open (filtered_out_path, 'w')
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
			SeqIO.write ([seq], filtered_out_hndl, "fasta")
			filter_cnt += 1
		else:
			log.debug ("Rejecting '%s' ..." % seq.id)
	filtered_out_hndl.close()
	## Postconditions & return:
	log.debug ('%s sequences remain after filtering ...' % filter_cnt)
	return filtered_out_path


def parse_args():	
	# Construct the option parser.
	from optparse import OptionParser
	optparser = OptionParser (
		prog          = 'blasaid',
		usage         = '%prog [options] SEQFILES ...',
		version       = 'version %s' %  VERSION,
		description   = "Reduce and blast NGS data.",
		# epilog='FORMAT must be one of %s.\n',
		epilog        = '',
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
		type=store_true,
		default=False,
		help='''Sequences that are identical are reduced to a single example''',
	)
	
	optparser.add_option ('--cluster-by-subsequence',
		dest="cluster_subsequence",
		type=store_true,
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
	
	optparser.add_option ('--output-format', '-o',
		dest="output_format",
		choices=OUTPUT_FORMATS,
		default=DEFAULT_OUTPUT_FORMAT,
		help='''The format of the output blast result files. If not supplied, this
			will default to %s.''',
		metavar='BLAST_FORMAT',
	)
	
	optparser.add_option ('--intermediate-files',
		dest="intermediate_files",
		type=str,
		default=None,
		help='''Where to store intermediate files. By default these will be kept
			with the system temporary files (e.g. /tmp).''',
		metavar='NAME',
	)
	
	optparser.add_option ('--dryrun',
		dest="dryrun",
		action='store_true',
		default=False,
		help='''Parse and filter the sequences but don't run the blast search'''
	)
	
	optparser.add_option ('--dump-options',
		dest="dump_options",
		action='store_true',
		default=False,
		help='''Print the current value of all program options'''
	)
	
	optparser.add_option ('--verbosity', '-v',
		dest="verbosity",
		default='4',
		help='''How much output to generate.'''
	)
	
	optparser.add_option ('--traceback',
		action='store_true',
		dest="traceback",
		default=False,
		help='''Exceptions are logged.''' # TODO
	)
	
	# parse 
	options, infiles = optparser.parse_args()
	
	### Postconditions & return:
	if (not infiles):
		optparser.error ('No input files specified')
		
	return infiles, options


def main():
	input_files, options = parse_args()

	try:
		scriptlog.init_logger (
			format="%(message)s",
			verbosity=options.verbosity,
		)
		
		if options.dump_options:
			dump_options (options)
		
		scratch_dir = get_scratch_dir (options.intermediate_files)
		log.info ("Making temporary files at '%s' ..." % scratch_dir)
		
		merge_file = merge_and_trim_seqs (input_files, scratch_dir,
			merge_quals=not options.ignore_qual_files,
			trim_right=options.trim_right,
			input_format=options.input_format
		)

		filter_file = filter_seqs (merge_file, scratch_dir,
			make_filters (
				length = options.filter_length,
				base_quality = options.filter_base_qual_threshold,
				avg_quality = options.filter_avg_qual_threshold,
			)
		)
		
		log.debug ('Clustering sequences ...')
		cluster_file = cluster_seqs (filter_file, scratch_dir,
			merge_quals=not options.ignore_qual_files,
			input_format=options.input_format)
		
	
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
