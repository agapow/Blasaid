"""
Assorted format related information.
"""


### IMPORTS

### CONSTS & DEFINES

FORMAT_LIST = {
	'clustal':     ['aln', 'clustal'],
	'fasta':       ['fas'],
	'fastq':       ['fastq-sanger'],
	'genbank':     ['gb'],
	'nexus':       ['nex', 'paup'],
	'stockholm':   ['sth'],

#	['tab'],
#	['csv'],
}

# what extensions map to what formats
FORMAT_MAP = {}
EXT_MAP = {}
for format, ext_list in FORMAT_LIST.iteritems():
	if (format not in ext_list):
		ext_list = [format] + ext_list
	for e in ext_list:
		EXT_MAP[e] = format
	FORMAT_MAP[format] = ext_list[0]


### IMPLEMENTATION ###

def get_format_from_extension (ext):
	"""
	Map an extension to a file format.
	
	:Returns:
		the file format in lowercase (e.g. 'fasta')
	
	Note that if the format is unknown, the extension is taken as a format.
	
	"""
	canon_ext = ext.lower()
	return EXT_MAP.get (canon_ext, canon_ext)
	

def get_extension_for_format (fmt):
	"""
	Get the canonical (preferred) extension for a format.
	
	:Returns:
		the extension in lower case (e.g. 'aln')
	
	Note that if the format is unknown, the format is taken as an extension.

	"""
	canon_fmt = fmt.lower()
	return FORMAT_MAP.get (canon_fmt, canon_fmt)


def format_has_quality (fmt):
	return fmt.startswith ('fastq')


### END #######################################################################
