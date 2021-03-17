#!/usr/bin/python

'''
Script to generate toy data for use with the
visualization tools.

Circumvents the need to use real IceCube data
in the visualizition plots

'''
from numba import jit

class PseudoPulse:
	def init(self, charge, time):
		self.charge = charge
		self.time = time

	@property
	def charge(self):
		return self._charge

	@property
	def time(self):
		return self._time

	@time.setter
	def time(self, value):
		assert value>0, 'ERROR: time must be positive!'
	

def to_hdf5(filename, input_dict, group='/', tablename='toydata'):
	'''
	Convert a python dict into hdf5 table

	inputs
	---

	filename: str (name of hdf5 name)
	
	input_dict: dict (the data in the form of a python dict)

	group: str (name of pytable group where table is stored)

	tablename: str (name of hdf5 table)

	'''

	import tables


def generate_pseudo_data(livetime, spe_peak, start_date):
	'''
	Create a multidimensional array of data 
	'''
	data_dict = {}

	# Generate a random trigger time in mjd

	# draw a reco energy from some distribution

	# draw a reco vertex from some distribution

	# draw charge from an SPE-like distribution

	# apply discriminator threshold

	# return a dict
	return data_dict


@jit
def discriminator_threshold(input_charge, discriminator_threshold=0.2):
	'''
	mimicks the role of the discriminator threshold
	by removing pulses 
	'''
	assert isinstance(input_data, dict), 'ERROR: input data must be a dict'

	for pulse in input_charge:
		if pulse.charge


if __name__=='__main__':

	import argparse
	parser = argparse.ArgumentParser('Generate toy data')
	parser.add_argument('-l', '--livetime', help='Livetime of data requested, in years', default=1.0)
	parser.add_argument('-spe', '--spe-peak', help='Define the location of the simulated SPE peak', default=1.0)
	parser.add_argument('-o', '--output-file', help='name of output .hd5 file', default='test_pseudo_data.hd5')

	args = parser.parse_args()

