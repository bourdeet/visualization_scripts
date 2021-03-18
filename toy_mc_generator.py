#!/usr/bin/python

'''
Script to generate toy data for use with the
visualization tools.

Circumvents the need to use real IceCube data
in the visualizition plots

'''
import numpy as np
import scipy as scp
from numba import jit
from astropy.time import Time, TimeDelta


# Baseline event rate in Hz
BASELINE_RATE = 1.0

# Set a default start date in mjd
DEFAULT_START = Time(2010., format='decimalyear').mjd

# default SPE peak location
DEFAULT_SPE = 1.029

def spe_functional_form(x, spe_peak=DEFAULT_SPE):

	# 
	# values taken from the example shown in 2002.00997
	# location of spe peak is made to vary
	# 
	P1 = 0.186
	w1 = 0.027
	P2 = 0.255
	w2 = 0.424
	mu = spe_peak
	sigma = 0.280

	expo_1 = P1/w1*np.exp(-x/w1)

	expo_2 = P2/w2*np.exp(-x/w2)

	gauss = (1-P1-P2)/sigma/np.sqrt(np.pi2)/scp.special.erfc(-mu/sigma/np.sqrt(2))*np.exp(-(x-mu)**2./(2*sigma**2.))

	return expo_1+expo_2+gauss

	

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


def generate_pseudo_data(livetime, spe_peak, start_date, rate=BASELINE_RATE):
	'''
	Create a multidimensional array of data 

	livetime: float (data livetime in seconds)

	rate: float (rate in Hz)
	'''
	data_dict = {}

	# Generate a random trigger time in mjd
	data_dict['time_mjd'] = draw_trigger_times(livetime, rate)
	n_events = data_dict['time_mjd'].size[0]

	# draw a reco energy from some distribution
	data_dict['reco_energy'] = draw_energy_distribution(n_events)

	# draw a reco vertex from some distribution
	x, y, z = draw_vertex_locations(n_events)
	data_dict['vertex_x'] = x
	data_dict['vertex_y'] = y
	data_dict['vertex_z'] = z

	# draw charge from an SPE-like distribution
	data_dict['charge'] = draw_charge_distribution(n_events, spe_peak=spe_peak)

	# apply discriminator threshold
	data_dict['mask'] = discriminator_threshold(data_dict['charge'])

	# return a dict
	return data_dict

#@jit
def draw_trigger_times(livetime, rate, mjd_offset=DEFAULT_START):
	'''
	livetime: float (livetime desired in seconds)

	rate: float (desired baseline rate in Hz)

	mjd_offset: float (mark the time of the beginning of the data period, in mjd)
	'''

	# Determine a number of events to draw, based on
	# Poisson expectation
	lamb=rate*livetime
	n_events = np.random.poisson(lam=lamb)

	# convert livetime into mjd Time Delta
	dT = TimeDelta(livetime, format='sec').mjd


	# Draw n_events number from a uniform random distribution
	time_series = []
	for i in range(n_events):
		t=np.random.uniform()
		time_series.append(mjd_offset+(dT*t))

	return np.sort(time_series)

#@jit
def draw_energy_distribution(n_events, mu=30., sigma= 10.):
	'''
	energy distribution is just a gaussian
	'''
	energy_array = []
	
	for i in range(n_events):
		E = np.random.normal(loc=mu, scale=sigma)
		energy_array.append(E)
	return energy_array

#@jit
def draw_vertex_locations(n_events, mu=0., sigma=120.):
	'''
	vertex locations drawn from a sphere
	of gaussian radial intensity
	'''

	x = []
	y = []
	z = []

	for i in range(n_events):
		r = np.random.normal(loc=mu, scale=sigma)
		theta = np.random.uniform(0.,np.pi)
		phi = np.random.uniform(0., 2*np.pi)

		x.append(r*np.sin(theta)*np.sin(phi))
		y.append(r*np.sin(theta)*np.cos(phi))
		z.append(r*np.cos(theta))
	
	return x, y, z

#@jit
def draw_charge_distribution(n_events, spe_peak=DEFAULT_SPE):
	'''

	'''
	integral = scp.integrate.quad(spe_functional_form, 0., 5.)
	max_value = spe_functional_form(0, spe_peak=spe_peak)/integral

	n=0
	charge = []
	while n<n_events:

		first_draw = np.random.uniform(0., 5.)
		second_draw = np.random.uniform(0., max_value)
		first_y = spe_functional_form(first_draw, spe_peak=spe_peak)

		if second_draw<=first_y:
			charge.append(first_draw)
			n+=0


	return charge




def discriminator_threshold(input_charge, discriminator_threshold=0.2):
	'''
	mimicks the role of the discriminator threshold
	by removing pulses 
	'''
	assert isinstance(input_data, dict), 'ERROR: input data must be a dict'

	mask = []
	for pulse in input_charge:
		if pulse<discriminator_threshold:
			mask.append(0)
		else:
			mask.append(1)
	return mask


if __name__=='__main__':

	import argparse
	parser = argparse.ArgumentParser('Generate toy data')
	parser.add_argument('-l', '--livetime', help='Livetime of data requested, in years', default=1.0)
	parser.add_argument('-spe', '--spe-peak', help='Define the location of the simulated SPE peak', default=1.0)
	parser.add_argument('-o', '--output-file', help='name of output .hd5 file', default='test_pseudo_data.hd5')

	args = parser.parse_args()

