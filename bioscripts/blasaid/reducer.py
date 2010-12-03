"""
Evaluate and filter Blast hits and HSPs.
"""


class BaseReducer (object):
	def __call__ (self, x):
		return True

class HitReducer (BaseReducer):
	pass

class HspReducer (BaseReducer):
	pass