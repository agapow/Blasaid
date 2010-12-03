"""
Assorted file utilities.
"""

### IMPORTS ###

from os import path
import tempfile


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###

def create_scratch_dir (name, pth=None):
	"""
	Create a directory for scratch files.
	
	Do this in temporary files unless otherwise prompted. The name will be used
	as the prefix of the directory name (i.e. a random suffix will be appended
	to prevent name collisions).
	"""
	if (pth):
		return tempfile.mkdtemp (prefix=name, dir=pth)
	else:
		return tempfile.mkdtemp (prefix=name)


def make_output_dir (name, pth=None):
	"""
	Create a directory for saving results.
	
	Do this in the current directory unless otherwise prompted. A random suffix
	will be appended to prevent name collisions.
	"""
	return create_scratch_dir (name, pth or '.')
	

def open_file_hndl (paths, mode):
	"""
	Join the path components, open the file there and return the handle.
	
	A convienience method.
	"""
	return open (path.join (*paths), mode)
	

def open_for_writing (paths, mode=''):
	"""
	Join the path components, open the file there and return the handle.
	
	A convienience method.
	"""
	return open_file_hndl (paths, 'w' + mode)
	

def open_for_reading (paths):
	"""
	Join the path components, open the file there and return the handle.
	
	A convienience method.
	"""
	return open_file_hndl (paths, 'r' + mode)


def write_to_file (paths, data, mode=''):
	hndl = open_for_writing (paths, mode)
	hndl.write (data)
	hndl.close()
	return hndl.name



### END #######################################################################
