Introduction
============

Processing next-generation sequencing (NGS) data can be arduous. When the data comes from an unknown or obscure organism, each of the plethora of reads has to be assessed and blasted so as to ascertain its identity. **Blasaid** (*BLAH-sade*) is a simple commandline program to:

* filter NGS data by length, quality and/or similarity
* automatically blast the remaining data against a database
* synopsize the results in a number of user-friendly formats


Installation
============

*Blasaid* [#homepage]_ can be installed in a number of ways. First, it is
made available as a Python package bioscripts.blasaid [#pypi_blasaid]_
First, it is made available as a Python package bioscripts.blasaid
[#pypi_blasaid]_, and may be installed via a number of automated methods using
setuptools [#setuptools]_, which will install all prerequsites. A manual installation will
also suffices, but prerequsiites must be installed beforehand .

Via setuptools / easy_install
-----------------------------

From the commandline call::

	% easy_install bioscripts.blasaid
	
Superuser privileges may be required. 


Via setup.py
------------

Download a source tarball, unpack it and call setup.py to
install::

	% tar zxvf bioscripts.blasaid.tgz
	% cd bioscripts.blasaid
	% python setup.py install
	
Superuser privileges may be required. 

Manual
------

Download and unpack the tarball as above. Ensure Biopython is available. Copy
the scripts in bioscripts.blasaid to a location they can be called from.


Usage
=====

Depending on your platform, scripts may be installed as ``.py`` scripts,
or some form of executable, or both.

The script is called::

	blasaid [options] SEQFILES ...

Options include:

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i SEQ_FORMAT, --input-format=SEQ_FORMAT
                        The format of the input sequence files. If not
                        supplied, this will be
                        inferred from the extension of the files.
  -o BLAST_FORMAT, --output-format=BLAST_FORMAT
                        The format of the output blast result files. If not
                        supplied, this                         will default to
                        %s.
  --filter-by-similarity=BLAST_FORMAT
                        Sequences that are more similar than this threshold
                        are clustered                         into a single
                        query
  --filter-by-length=LENGTH
                        Sequences that are shorter than this are discarded
  --filter-by-base-quality=QUALITY
                        Sequences that contain any nucleotide with a quality
                        lower than                         this threshold are
                        discarded
  --filter-by-average-quality=QUALITY
                        Sequences with an average quality lower than
                        this threshold are discarded
  --intermediate-file=NAME
                        The file to save filtered sequences into
  -v VERBOSITY, --verbose=VERBOSITY
                        How much output to generate.
  --traceback           Exceptions are logged.


Developer notes
===============

This module isn't intended primarily for importing, but the setuptools packaging and
infrastructure make for simple distribution of scripts, allowing the checking
of prerequisites, consistent installation and updating.

The ``bioscripts`` namespace was chosen as a convenient place to "keep" these
scripts and is open to other developers.


References
==========

.. [#homepage] `bioscripts.blasaid homepage <http://www.agapow/net/software/bioscripts.convert>`__

.. [#biopython] `Biopython homepage <http://www.biopython.org>`__

.. [#setuptools] `Installing setuptools <http://peak.telecommunity.com/DevCenter/setuptools#installing-setuptools>`__



