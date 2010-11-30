"""
Parse a Blast results expressed as XML, allow deletion of hits and output.
"""

### IMPORTS

from cStringIO import StringIO
from xml.etree import cElementTree as elemtree

									 
### CONSTANTS & DEFINES

### IMPLEMENTATION ###

class BlastReader (object):
	def __init__ (self, xml):
		# if it's a string, stuff into a buffer because etree only reads xml
		# documents from file-likes
		if isinstance (xml, basestring):
			xml = StringIO (xml)
		self.doc_tree = elemtree.parse (xml)
		self.doc_root = self.doc_tree.getroot()
		self.params = self._init_header()
		# we assume there's only a single iteration in this tree
		self.hits_elem = r.find ('.//Iteration_hits')
		
	def _init_header (self):
		params = BlastParams()
		params.program = self.doc_root.findtext ('BlastOutput_program').strip()
		params.version = self.doc_root.findtext ('BlastOutput_version').strip()
		params.db = self.doc_root.findtext ('BlastOutput_db').strip()
		params.query_len = self.doc_root.findtext ('BlastOutput_query-len').strip()
		params_node = self.doc_root.find ('.//Parameters').strip()
		params.params = dict ([(n.tag.split('_', 1)[1], n.text) for n in params_node])
		return params
			
	
	# TODO: for Python 3.0
	# __next__ = next
	
	def read_hits (self):
		"""
		Returns the node in the tree
		"""
		for n in self.hits_elem.getchildren():
			yield n
		
	def __len__ (self):
		# return number of hits
		return len (self.hits_elem)


class BlastParams (object):
	"""
	Just a container for Blast call details"
	"""
	def __init__ (self):
		self.program = None
		self.version = None
		self.db = None
		self.query_len = None
		self.params = None


class BlastHit (object):
	"""
	Just a container for Blast hit details"
	"""
	def __init__ (self):
		self.bit_score = None
		self.score = None
		self.e_value = None
		self.query_align = None
		self.hit_align = None
		self.align_len = None
		self.identity = None
		self.positive = None
		self.gaps = None
		self.query_seq = None
		self.hit_seq = None


class BlastHsp (object):
	"""
	Just a container for Blast hit details"
	"""
	def __init__ (self):
		self.bit_score = None
		self.score = None
		self.e_value = None
		self.query_align = None
		self.hit_align = None
		self.align_len = None
		self.identity = None
		self.positive = None
		self.gaps = None
		self.query_seq = None
		self.hit_seq = None
		
		
		
def hit_node_to_obj (node):
	
	



### END #######################################################################
