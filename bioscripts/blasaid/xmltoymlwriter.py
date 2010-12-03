"""
A writer that takes XML and maps it to (easier to read) YAML.
"""

### IMPORTS

from xml.etree import cElementTree as elemtree



### CONSTANTS & DEFINES

### IMPLEMENTATION ###

class XmlToYamlWriter (object):
	def __init__ (self, path_or_hndl):
		if isinstance (path_or_hndl, basestring):
			self.hndl_opened = True
			hndl_or_path = open (path_or_hndl, 'w')
		else:
			self.hndl_opened = False
		self.hndl = hndl_or_path
			
	def __del__ (self):
		if self.hndl_opened:
			self.hndl.close()		
			
	def write_xml_doc (self, file_hndl_or_xml):
		if isinstance (file_hndl_or_xml, basestring):
			hndl_opened = True
			file_hndl_or_xml = open (file_hndl_or_xml, 'r')
		else:
			hndl_opened = False
		if hasattr (file_hndl_or_xml, 'read'):
			file_hndl_or_xml = open (file_hndl_or_xml, 'r')		
		doc_tree = elemtree.parse (file_hndl_or_xml)
		doc_root = doc_tree.getroot()
		self._write_doc (self, doc_root)
		if hndl_opened:
			file_hndl_or_xml.close()
			
	def _write_doc (self, doc_root):
		self.hndl.write ('--- !%s\n' % doc_root.tag)
		self.hndl.write (yaml.dump (
			[self._convert_node (n) for n in doc_root.get_children()]
		))
		self.hndl.write ('...\n')
		
	def _convert_node (self, node):
		meth_name = '_convert_%s_node' % node.tag
		if hasattr (self, meth_name):
			return getattr (self, meth_name)(node)
		else:

	
### END #######################################################################
