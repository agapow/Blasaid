"""
Methods for selecting sequences based on different criteria.
"""


class MinLength (object):
	"""
	Reject sequences falling below a certain length.
	"""
	
	def __init__ (self, threshold):
		self.threshold = threshold
		
	def __call__ (self, seq):
		return self.threshold <= len (seq)
		

class MinBaseQuality (object):
	"""
	Reject sequences with any bases below a certain quality.
	
	Sequences without a quality are let through.
	"""
	
	def __init__ (self, threshold):
		self.threshold = threshold
		
	def __call__ (self, seq):
		# if there's a quality, check it
		quals = seq.letter_annotations
		if quals:
			# check phred quality or else whatever they give
			scores = quals.get ('phred_quality') or quals.get (quals.keys[0])
			for s in scores:
				if (s < self.threshold):
					return False
		return True


class MinAvgBaseQuality (object):
	"""
	Reject sequences with an average quality below a given threshold.
	
	Sequences without a quality are let through.
	"""
	
	def __init__ (self, threshold):
		self.threshold = threshold
		
	def __call__ (self, seq):
		# if there's a quality, check it
		quals = seq.letter_annotations
		if quals:
			# check phred quality or else whatever they give
			scores = quals.get ('phred_quality') or quals.get (quals.keys[0])
			avg_qual = sum (scores) / float (len (scores))
			if (avg_qual < self.threshold):
				return False
		return True
		