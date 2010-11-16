"""
Script for automated blasting of NGS data.
"""

__docformat__ = 'restructuredtext en'
__author__ = 'Paul-Michael Agapow <paul-michael.agapow@hpa.org.uk>'


### IMPORTS ###

from os import path
from exceptions import BaseException
import logging

from Bio import AlignIO

from bioscripts.blasaid import filter


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

_DEV_MODE = True


### IMPLEMENTATION ###

def init_logger ()
	global (logger)
	logger = logging.getLogger("main")
	handler = logging.StreamHandler()
	handler.setLevel(logging.DEBUG)
	
def make_filters (length=None, base_quality=None, avg_quality=None):
	filter_list = []
	if (length):
		filter_list.append (filter.MinLength (float(length)))
	if (base_quality):
		filter_list.append (filter.MinBaseQuality (float(base_quality)))
	if (avg_quality):
		filter_list.append (filter.MinAvgBaseQuality (float(avg_quality)))
	return filter_list
	
	
def parse_args (arg_arr):
	init_logger()
	
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
	# TODO: filter by range
	
	optparser.add_option ('--input-format', '-i',
		dest="input_format",
		choices=INPUT_FORMATS,
		default=DEFAULT_INPUT_FORMAT,
		help='''The format of the input sequence files. If not supplied, this will be
			inferred from the extension of the files.''',
		metavar='SEQ_FORMAT',
	)
	
	optparser.add_option ('--output-format', '-o',
		dest="output_format",
		choices=OUTPUT_FORMATS,
		default=DEFAULT_OUTPUT_FORMAT,
		help='''The format of the output blast result files. If not supplied, this
			will default to %s.''',
		metavar='BLAST_FORMAT',
	)
	
	optparser.add_option ('--filter-by-similarity', '',
		dest="filter_similar",
		type=float,
		default=None,
		help='''Sequences that are more similar than this threshold are clustered
			into a single query''',
		metavar='BLAST_FORMAT',
	)
	
	optparser.add_option ('--filter-by-length', '',
		dest="filter_length",
		type=int,
		default=None,
		help='''Sequences that are shorter than this are discarded''',
		metavar='LENGTH',
	)
	
	optparser.add_option ('--filter-by-base-quality', '',
		dest="filter_base_qual_threshold",
		type=float,
		default=None,
		help='''Sequences that contain any nucleotide with a quality lower than
			this threshold are discarded''',
		metavar='QUALITY',
	)

	optparser.add_option ('--filter-by-average-quality', '',
		dest="filter_avg_qual_threshold",
		type=float,
		default=None,
		help='''Sequences with an average quality lower than
			this threshold are discarded''',
		metavar='QUALITY',
	)

	optparser.add_option ('--intermediate-file', '',
		dest="intermediate_file",
		type=str,
		default='filtered-seqs.fasta',
		help='''The file to save filtered sequences into''',
		metavar='NAME',
	)
	
	
	optparser.add_option ('--dryrun',
		dest="dryrun",
		help='''Parse and filter the sequences but don't run the blast search'''
	)
	
	optparser.add_option ('--verbose', '-v',
		dest="verbosity",
		help='''How much output to generate.'''
	)
	
	optparser.add_option ('--traceback',
		action='store_true',
		dest="traceback",
		default=False,
		help='''Exceptions are logged.''' # TODO
	)
	
	options, pargs = optparser.parse_args (args=arg_arr)
	
	infiles = pargs[1:]
	if (not infiles):
		optparser.error ('No input files specified')
	
	return infiles, options



def main():
	import sys
	input_files, options = parse_args (sys.argv[1:])
	filters = make_filters (
		length = options.filter_length,
		base_quality = options.filter_base_qual_threshold,
		avg_quality = options.filter_avg_qual_threshold,
	)
	for f in input_files:
		if (options.format == 'auto'):
			rdr = ExSeqReader (f)
		else:
			rdr = ExSeqReader (f, options.format)
		for seq in rdr:
			if (all ([f(seq) for f in filters])):
				
			
		
		
	

if __name__ == '__main__':
	main()


### END #######################################################################
