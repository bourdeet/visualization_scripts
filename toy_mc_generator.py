#!/usr/bin/python

'''
Script to generate toy data for use with the
visualization tools.

Circumvents the need to use real IceCube data
in the visualizition plots

'''
from collections import OrderedDict
import numpy as np
import scipy as scp
import scipy.integrate as integrate
from numba import jit
from astropy.time import Time, TimeDelta


# Baseline event rate in Hz
DEFAULT_RATE = 0.001

# Set a default start date in mjd
DEFAULT_START = Time(2010., format='decimalyear')

# default SPE peak location
DEFAULT_SPE = 1.029

def get_dtype(data):
    '''
    Given a dict, generate a nested numpy dtype
    Stolen from: https://gist.github.com/blink1073/feb3ac5920653ae222c2
    '''
    fields = []
    for (key, value) in list(data.items()):
        # make strings go to the next 64 character boundary
        # pytables requires an 8 character boundary
        if isinstance(value, str):
            value += ' ' * (64 - (len(value) % 64))
            # pytables does not support unicode
            if isinstance(value, str):
                value = value.encode('utf-8')
        elif isinstance(value, str):
            value += ' ' * (64 - (len(value) % 64))

        if isinstance(value, dict):
            fields.append((key, get_dtype(value)))
        else:
            value = np.array(value)
            fields.append((key, '%s%s' % (value.shape, value.dtype)))
    return np.dtype(fields)

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

    gauss = (1-P1-P2)/sigma/np.sqrt(np.pi/2)/scp.special.erfc(-mu/sigma/np.sqrt(2))*np.exp(-(x-mu)**2./(2*sigma**2.))

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


    first_row_dict = OrderedDict([ (k,v[0]) for k,v in list(input_dict.items()) ])
    dtype = get_dtype(first_row_dict)


    file = tables.open_file(filename, mode='w')
    new_table = file.create_table(where=group, name=tablename, description=dtype)
    
    n_rows = len(list(input_dict.values())[0])
    for i in range(n_rows):
        r = new_table.row 
        for k,v in input_dict.items():
            r[k] = v[i] 
        r.append()
        new_table.flush()

    file.close()
    print('table stored in ',filename)


def generate_pseudo_data(livetime, spe_peak, start_date, rate=DEFAULT_RATE, debug=False):
    '''
    Create a multidimensional array of data 

    livetime: float (data livetime in seconds)

    rate: float (rate in Hz)
    '''
    data_dict = {}

    # Generate a random trigger time in mjd
    if debug:
        print('drawing trigger times...')
    data_dict['time_mjd'] = draw_trigger_times(livetime, rate)
    n_events = len(data_dict['time_mjd'])

    assert n_events>0, 'ERROR: your chosen livetime and rate yields no data'

    # draw a reco energy from some distribution
    if debug:
        print('drawing energies...')
    data_dict['reco_energy'] = draw_energy_distribution(n_events)

    # draw a reco vertex from some distribution
    if debug:
        print('drawing vertex locations...')
    x, y, z = draw_vertex_locations(n_events)
    data_dict['vertex_x'] = x
    data_dict['vertex_y'] = y
    data_dict['vertex_z'] = z

    # draw charge from an SPE-like distribution
    if debug:
        print('drawing charge...')
    data_dict['charge'] = draw_charge_distribution(n_events, spe_peak=spe_peak)

    # apply discriminator threshold
    if debug:
        print('applying discriminator threshold...')
    #data_dict['mask'] = discriminator_threshold(data_dict['charge'])

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
    dT = TimeDelta(livetime, format='sec')


    # Draw n_events number from a uniform random distribution
    time_series = []
    for i in range(n_events):
        t=np.random.uniform()
        mjd_t = mjd_offset+TimeDelta(dT.sec*t,format='sec')
        time_series.append(mjd_t.value)

    return list(np.sort(time_series))

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
    integral = integrate.quad(spe_functional_form, 0., 5.)[0]
    max_value = spe_functional_form(0, spe_peak=spe_peak)/integral


    n=0
    charge = []
    while n<n_events:

        first_draw = np.random.uniform(0., 5.)
        second_draw = np.random.uniform(0., max_value)
        first_y = spe_functional_form(first_draw, spe_peak=spe_peak)

        if second_draw<=first_y:
            charge.append(first_draw)
            n+=1


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
    parser.add_argument('-l', '--livetime', help='Livetime of data requested, in years', type=float, default=1.0)
    parser.add_argument('-sd','--start-day', help='provide a reference start date in mjd', type=float, default=DEFAULT_START)
    parser.add_argument('-spe', '--spe-peak', help='Define the location of the simulated SPE peak', type=float, default=DEFAULT_SPE)
    parser.add_argument('-r', '--rate', help='define a raw event rate priot to any trigger', type=float, default=DEFAULT_RATE)
    parser.add_argument('-o', '--output-file', help='name of output .hd5 file', default='test_pseudo_data.hd5')

    parser.add_argument('--debug', help='debug mode', action='store_true')

    args = parser.parse_args()

    # Approximate conversion from years to seconds
    if args.debug:
        livetime_s = np.pi*(10**7)*0.005
    else:
        livetime_s = np.pi*(10**7)*args.livetime

    pseudo_data = generate_pseudo_data(livetime=livetime_s,
                                       spe_peak=args.spe_peak,
                                       start_date=args.start_day,
                                       rate=args.rate,
                                       debug=args.debug)

    to_hdf5(filename=args.output_file,
            input_dict=pseudo_data)

    print('Done')
