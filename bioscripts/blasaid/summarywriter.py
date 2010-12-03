"""
Writes a user-friendly summary of blast results.
"""
# TODO: other formats?

### IMPORTS

from datetime import datetime

from bioscripts.blasaid import __version__


### CONSTANTS & DEFINES

INDENT = '   '


### IMPLEMENTATION ###

class SummaryWriter (object):
	def __init__ (self, pth):
		self.hndl = open (pth, 'w')
		self.write_header()
		
	def __del__ (self):
		self.write_footer()
		self.hndl.close()
		
	def write_header (self):
		self.hndl.write ("# Blast summary written by blasaid v%s, %s\n" % (
			__version__,
			datetime.now().strftime ("%Y%m%dT%H%M")
		))			
	
	def write_result (self, hit):
		self.hndl.write ("---\n")
		self._write_params (hit)
		self.hndl.write ('hits:\n')
		for d, a in zip (hit.descriptions, hit.alignments):
			self._write_hit (d, a)
			
	def _write_params (self, hit):
		self.hndl.write ('application:\n')
		self.hndl.write (INDENT + 'name: %s\n' % hit.application)
		self.hndl.write (INDENT + 'version: %s\n' % hit.version)
		self.hndl.write (INDENT + 'query: %s\n' % hit.query)
		self.hndl.write (INDENT + 'database: %s\n' % hit.database)
	
		self.hndl.write ('parameters:\n')
		self.hndl.write (INDENT + 'gap-penalties: %s\n' % str (hit.gap_penalties))
		self.hndl.write (INDENT + 'blast-cutoff: %s\n' % str (hit.blast_cutoff))
		self.hndl.write (INDENT + 'window-size: %s\n' % hit.window_size)
		self.hndl.write (INDENT + 'threshold: %s\n' % hit.threshold)
	
	def _write_hit (self, desc, aln):
		self.hndl.write (' - ' + 'title: %s\n' % desc.title)
		self.hndl.write (INDENT + 'score: %s\n' % desc.score)
		self.hndl.write (INDENT + 'expect: %s\n' % desc.e)
		self.hndl.write (INDENT + 'hsps:\n')
		for h in aln.hsps:
			self._write_hsp (h)		
	
	def _write_hsp (self, hsp):
		self.hndl.write (INDENT + ' - ' + 'score: %s\n' % hsp.score )
		self.hndl.write (INDENT + INDENT + 'bits: %s\n' % hsp.bits )
		self.hndl.write (INDENT + INDENT + 'expect: %s\n' % hsp.expect )
		self.hndl.write (INDENT + INDENT + 'identities: %s\n' % hsp.identities )
		self.hndl.write (INDENT + INDENT + 'positives: %s\n' % hsp.positives )
		self.hndl.write (INDENT + INDENT + 'gaps: %s\n' % hsp.gaps )
		self.hndl.write (INDENT + INDENT + 'query-alignment: %s-%s\n' % (hsp.query_start, hsp.query_end))
		self.hndl.write (INDENT + INDENT + 'hit-alignment: %s-%s\n' % (hsp.sbjct_start, hsp.sbjct_end))
		self.hndl.write (INDENT + INDENT + 'alignment: |\n')
		self.hndl.write (INDENT * 3 + '%s   %s\n' % ('query'.ljust(16), hsp.query))
		self.hndl.write (INDENT * 3 + '%s   %s\n' % (''.ljust(16), hsp.match))
		self.hndl.write (INDENT * 3 + '%s   %s\n' % ('subject'.ljust(16), hsp.sbjct))
			

	def write_footer (self):
		self.hndl.write ("...\n")
	
			
		


### END #######################################################################
