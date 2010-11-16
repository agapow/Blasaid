"""
A reader for sequence files.

Just does a few extra things on top on BioPython standard routines:

* opens file handle
* auto-detects format from extension
* merges quality file if present

"""

### IMPORTS

from os.path import splitext, exists, join
from itertools import izip

from Bio import SeqIO

import formats


### CONSTANTS & DEFINES

### IMPLEMENTATION

class ExSeqReader (object):
	def __init__ (self, path_or_hndl, fmt=None):
		# open any passed filepaths
		if (isinstance (path_or_hndl, basestring)):
			self.hndl = open (path_or_hndl, 'r')
			self.opened_hndl = True
		else:
			self.hndl = path_or_hndl
			self.opened_hndl = False
		# set format, guess from extension if need be
		if (fmt):
			self.fmt = fmt
		else:
			ext = splitext (self.hndl.name)[1][1:]
			assert (ext), "can't deduce format without extension on '%s'" % hndl.name
			self.fmt = formats.get_format_from_extension (ext)
		# if there's a qual file, read it
		self.qual_hndl = None
		if (not format_has_quality (self.fmt)):
			qual_path = join (splitext (self.hndl.name)[0], '.qual')
			if exists (qual_path):
				self.qual_hndl = open (qual_path, 'r')
	
	def __iter__ (self):
		return self
	
	def next():
		if (self.qual_hndl):
			for s, q in izip (SeqIO.parse (self.hndl, self.fmt), SeqIO.parse (self.qual_hndl, 'qual')):
				s.letter_annotations = q.letter_annotations
		else:
			for s in SeqIO.parse (self.hndl, self.fmt):
				yield s			

			



### END #######################################################################
