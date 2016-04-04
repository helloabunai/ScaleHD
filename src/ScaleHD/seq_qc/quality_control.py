##
## Generals
import glob
import os
import logging as log

##
## Backend junk
from ..backpack import Colour as clr

class SeqQC:

	def __init__(self, input_data, instance_rundir, stage):
		self.input_data = input_data
		self.instance_rundir = instance_rundir

		if stage == 'valid':
			self.verify_input()
		if stage == 'dmpx':
			self.execute_demultiplex()
		if stage == 'trim':
			self.execute_trimming()

	def verify_input(self, raise_exception=True):
		for fqfile in glob.glob(os.path.join(self.input_data, '*')):
			if fqfile.endswith('.fq') or fqfile.endswith('.fastq'):
				return True
		if raise_exception:
			log.error('{}{}{}{}'.format(clr.red,'shd__ ',clr.end,'I/O: Non-FQ/FastQ file found in input directory.'))
		return False

	def execute_demultiplex(self):
		print 'demultiplex here', self.input_data

	def execute_trimming(self):
		print 'trimming here', self.input_data




