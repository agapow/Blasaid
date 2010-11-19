"""
Simple log utilities and constants for the package.
"""

### IMPORTS

import logging as log


### CONSTANTS & DEFINES
 
# the various verbosity (log) levels and their names
VERBOSITY_LVLS_AND_VOCAB = {
	0                  : ['6', 'all'],
	log.DEBUG          : ['5', 'debug'],
	log.INFO           : ['4', 'info'],
	log.WARNING        : ['3', 'warning'],
	log.ERROR          : ['2', 'error'],
	log.CRITICAL       : ['1', 'critical'],
	log.CRITICAL+100   : ['0', 'none'],
}

# mapping verbosity vocab to levels
VERBOSITY_VOCAB_TO_LVL = {}
for lvl, vocab in VERBOSITY_LVLS_AND_VOCAB.iteritems():
	for vocab_item in vocab:
		VERBOSITY_VOCAB_TO_LVL[vocab_item] = lvl


### IMPLMENTATION ###

def init_logger (**kwargs):
	# NOTE: any calling or setting to the defautl logger will instantiate and
	# 'fix' it, so that future settings silently have no effect. Thus, you get
	# one shot at configuring the default logger.
	if kwargs.has_key ('verbosity'):
		kwargs['level'] = VERBOSITY_VOCAB_TO_LVL[kwargs['verbosity']]
	log.basicConfig (**kwargs)



### END #######################################################################
