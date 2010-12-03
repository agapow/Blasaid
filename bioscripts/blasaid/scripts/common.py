"""
Common utilities for scripts.
"""

### IMPORTS ###

import optparse




from os import path, mkdir
import sys
from exceptions import BaseException, SystemExit
import traceback
import tempfile
from datetime import datetime
from itertools import islice

from Bio import AlignIO, SeqIO
from Bio.Blast import NCBIXML

from bioscripts.blasaid import filter, scriptlog, cluster
from bioscripts.blasaid.scriptlog import log
from bioscripts.blasaid.exseqreader import ExSeqReader
from bioscripts.blasaid.blast import blast_ncbi
from bioscripts.blasaid.summarywriter import SummaryWriter


### CONSTANTS & DEFINES ###

try:
	from bioscripts.blasaid import __version__ as VERSION
except:
	VERSION = 'unknown'


### IMPLEMENTATION ###
	
class OptionParser (optparse.OptionParser):
	"""
	An option parser with common options for all blasaid-family scripts.
	"""
	def __init__ (self, **kwargs):
		kwargs['usage']    = '%prog [options] SEQFILES ...',
		kwargs['version']  = 'version %s' %  VERSION,		

	def add_results_path_option (self):
		# TODO: call this output path?
		self.add_option ('--results-path',
			dest="results_path",
			type=str,
			default=None,
			help='''Where blast results and summaries will be saved. By default this
				will be in the current directory.''',
			metavar='NAME',
		)
	
	def add_scratch_path_option (self):
		self.add_option ('--scratch-path',
			dest="scratch_path",
			type=str,
			default=None,
			help='''Where to store intermediate files. By default these will be kept
				with the system temporary files (e.g. /tmp).''',
			metavar='NAME',
		)
	
	def add_dryrun_option (self):
		self.add_option ('--dryrun',
			dest="dryrun",
			action='store_true',
			default=False,
			help='''Parse and filter the sequences but don't run the blast search'''
		)

	def add_dump_option (self):
		self.add_option ('--dump-options',
			dest="dump_options",
			action='store_true',
			default=False,
			help='''Print the current value of all program options'''
		)
	
	def add_verbosity_option (self):
		self.add_option ('--verbosity', '-v',
			dest="verbosity",
			default='4',
			help='''How much output to generate.'''
		)

	def add_traceback_option (self):
		self.add_option ('--traceback',
			action='store_true',
			dest="traceback",
			default=False,
			help='''Exceptions are logged.''' # TODO
		)
	
	def add_debug_options (self):
		self.add_dryrun_option()
		self.add_dump_option()
		self.add_verbosity_option()
		self.add_traceback_option()
		
	def parse_args (self):
		options, infiles = optparse.OptionParser.parse_args (self)
		validate_args (options, infiles)
		return options, infiles
	
	def validate_args (self, options, infiles):
		# need at least some infiles
		if (not infiles):
			optparser.error ('No input files specified')


### END #######################################################################
