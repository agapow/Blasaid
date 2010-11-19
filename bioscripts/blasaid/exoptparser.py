"""
An commandline option parser that wraps and extends the standard implementation.

Although the option parser in the standard Python library is fine, there are
several shortcomings:

* `optparse` is deprecated as of Python 2.7 to be replaced by `argparse`. While
  the latter is an improvement, we cannot yet rely on having 2.7. At the same
  time, we would like to be ready for Python 3.0.

* Some of the common argument types are verbose to define and arguably fail to
  convey well what they are doing.
  
* Some obvious flexibility (e.g. synonymous commands) is missing.

`ExOptParser` wraps the best available option parser and adds missing functionality
  
  
"""
# ??? allow just names to be added ('d' 'dest') and handled appropriately?
# ??? positionals
# ??? synonyms for options
# ??? synonyms for choices
# ??? postvalidation & auto postvalidation
# ??? auto postvalidation that most options are used twice?
# only error out when options are unsaitisfiable or contradictory, when things
# are merely pointless or nonsensical (a save name when you aren't saving, or a
# filter that select no sequences) be quiet
# ???: log as yaml?
# sensible default args: empty list for append

# seqtool filter rename count transform edit convert sort unique merge info

# --defaults: prints out default values and then quits

### IMPORTS

from optparse import OptionParser

__all__ = [
	'ExOptionParser',
]


### CONSTANTS & DEFINES

BOOLEAN_SYNONYMS = {
	'true'    : True,
	't'       : True,
	'yes'     : True,
	'y'       : True,
	'1'       : True,
	'false'   : False,
	'f'       : False,
	'no'      : False,
	'n'       : False,
	'0'       : False,
	'none'    : False,
	'null'    : False,
	'nil'     : False,
}


### IMPLEMENTATION ###

class CallableDict (dict):
	def __call__ (self, key):
		return self.__get__ (key)
		
		
class ExOptionParser (OptionParser):

	def __init__ (self, prog_name=None, desc=None, long_desc=None, usage=None):
		self.prog_name = prog_name
		pass
	
	def add_option (self, name, options, default):
		OptionParser.add_option (*_format_options_arg (options),
			dest=name,
			help=desc,
			default=default
		)
	
	def add_flag_option (self, name, options, default=False, store_val=True, 
			desc=None, metavar=None):
		"""
		Add an option that sets a value when seen.
		"""
		pass

	def add_value_option (self, name, options, default=None, convert=None,
			desc=None, metavar=None):
		"""
		Add an option that stores a passed value.
		
		:Parameters:
			convert
				A callable, or list of callables, that when passed the option
				value validates and/or converts it to a more suitable form.
				
		
		"""
		pass
	
	def add_append_value_option (self, name, options, default=[], convert=None,
			desc=None, metavar=None):
		"""
		Add an option that stores a passed value, and can be called multiple times.
		"""
		pass	
		
	def add_choice_option (self, name, options, choices, default=None,
			convert=None, desc=None, metavar=None):
		"""
		Add an option that is a value from a list. 
		"""
		pass

	def add_append_choice_option (self, name, options, choices, default=[],
			convert=None, desc=None, metavar=None):
		"""
		Add an option that is a value from a list, and can be called multiple times.
		"""
		pass
		# ??? use append, or use comma
		
	def add_boolean_option (self, name, options, default, desc, metavar):
		"""
		Add an option that can accept 'true' or 'false' values or various synonyms.
		
		This options accepts (in a case insensitive way):
		
		* true-false or t-f
		* yes-no or y-n
		* on-off
		* 1-0
		* ?-none/nil/null
		
		"""
		convert = lambda x: BOOLEAN_SYNONYMS[x.lower()]
		return self.add_choice_option (self, options, [True, False], default=default,
			convert=convert, desc=desc, metavar=metavar)
	
	def add_toggle_option (self, name, options, default, desc, metavar):
		"""
		Add two options to switch something on or off.
		
		Generates two options "--foo" and "--nofoo". Suggest that single letters
		are ruled out.
		"""
		pass
		
		
	def add_save_as_option (self, options=['-s', '--save-as'],
			dest="save_as",
			default="{dir}{stem}-out{dir}",
			help="name output files in this form",
			metavar='STR',
			tmpl="{dir}{substem}-{progname}{ext}",
			subs=None,
			interpolator=None,
		):
		pass
		
	def add_output_format_option (self, choices, options=['-o', '--output-format'], 
			dest="output_format",
			action="store",
			default=None,
			help="save output in this format",
			metavar='STR',
		):
		default = default or choices[0]
	
	def add_input_format_option (self, choices, options=['-i', '--input-format'], 
			dest="input_format",
			action="store",
			default=None,
			help="expect input to be in this format",
			metavar='STR',
		):
		default = default or choices[0]
		
	def add_debug_option (self, options=['--debug'],
			dest="debug",
			default=False,
			help="raise a full traceback on errors",
		):
		pass

	def add_verbosity_option (self, options=['v', 'verbosity'],
			dest="input_format",
			action="store",
			choices=[],
			default=None,
			desc="expect input to be in this format",
			metavar='STR',
		):
		default = default or choices[0]
	
	def add_positional_argument (self, number=1, default=None, desc=None, metavar=None):
		pass

	def add_trailing_argument (self, default=None, desc=None, metavar=None):
		pass	
		
	def post_validate (self):
		pass
	
	
### UTILITIES

def _make_sequence (x):
	"""
	If the argument is not a sequence, put it in one.
	
	:Returns:
		The argument in a list, unless it is a list or sequence already.
		
	This is just a covenience function for method signatures, such that a value
	or list of values may be used interchangeably.
	
	For example:
	
		>>> _make_sequence (1)
		[1]
		>>> _make_sequence ("foo")
		["foo"]
		>>> _make_sequence ([1])
		[1]
		>>> _make_sequence (["foo", "bar"])
		["foo", "bar"]
		
	"""
	if isinstance (x, (list, tuple)):
		return x
	else:
		return [x]


def _format_option_string (x):
	"""
	Convert options to the correct (leading hyphen) form, if need be.
	
	:Returns:
		An option string in the canonical (leading hyphen) form.
		
	This is a convienience function, that allows options to to be be described
	as single letters or words and the appropriate number of hyphens tagged on
	the front. For backwards compatiability, options that already begin with a
	hyphen are left alone.
	
	For example:
	
		>>> _format_option_string ('f')
		'-f'
		>>> _format_option_string ('foo')
		'--foo'
		>>> _format_option_string ('--f')
		'--f'
		
	"""
	if (x.startswith('-')):
		return x
	else:
		if (len(x) == 1):
			return "-%s" % x
		else:
			return "--%s" % x


def _format_options_arg (x):
	"""
	Convert an options definition to the correct form, if need be.
	
	:Returns:
		A list of options in the canonical (leading hyphen) form.
		
	This is a convienience function, that allows options to to be be described
	as a single value or a list, with leading hyphens or without.
	
	For example:
	
		>>> _format_options_arg ("f")
		['-f']
		>>> _format_options_arg (["foo"])
		['--foo']
		>>> _format_options_arg (["-foo", "bar"])
		['-foo', '--bar']
		
	"""
	return [_format_option_string (a) for a in _make_sequence (x)]
	

	
### END #######################################################################
