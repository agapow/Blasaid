Introduction
============

While next-generation sequencing (NGS) provides a bounty of data for analysis,
identifying sequence fragments without a reference genome is arduous. Each
fragment must be blasted seperately against a reference database, and the run
will include duplicates, largely overlapping sequences or those that are too
short or of too low quality. **blasaid** is a tool to helps reduce sequencing
runs to the more important or interesting reads, blasts each automatically and
synosizes the results for easy reading.

Each *blasaid* run consists of the following steps:

1. **Merge** all input sequences, incorporating any accompanying quality files.
Optionally **trim** 3' end of sequence to rid low quality trailing segments.

2. Optionally **filter** the input sequences, getting rid of any that fall below
thresholds in quality or length.

3. Optionally **cluster** the input sequences, reducing similar reads to a
single example

4. **Blast** the remaining reads against a chosen database and return results in
a given format.

5. The results are optionally **reduced** to the most interesting and
**summarized**.


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
  --ignore-qual-files   Do not merge qual files into sequence files.
  
  --trim-right=LENGTH   Trim the 3' end of sequences by this many bases
  
  --filter-by-length=LENGTH
                        Sequences that are shorter than this are discarded
  --filter-by-base-quality=QUALITY
                        Sequences that contain any nucleotide with a quality
                        lower than                         this threshold are
                        discarded
  --filter-by-average-quality=QUALITY
                        Sequences with an average quality lower than
                        this threshold are discarded
								
  --cluster-by-identity
                        Sequences that are identical are reduced to a single
                        example
  --cluster-by-subsequence
                        Sequences that are subsequences of others are
                        eliminated
  --cluster-by-similarity=PERCENT
                        Sequences that are more similar than this threshold
                        are clustered                         into a single
                        query
								
  -o BLAST_FORMAT, --output-format=BLAST_FORMAT
                        The format of the output blast result files. If not
                        supplied, this                         will default to
                        %s.
  --intermediate-files=NAME
                        Where to store intermediate files. By default these
                        will be kept                         with the system
                        temporary files (e.g. /tmp).
								
  --dryrun              Parse and filter the sequences but don't run the blast
                        search
  --dump-options        Print the current value of all program options
  -v VERBOSITY, --verbosity=VERBOSITY
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



