'''
Script to perform a comparison of seasonal
data.


Etienne Bourbeau
'''

# Standard Python stuff
import os, sys, collections, glob, datetime
import numbers
import numpy as np
from astropy.time import Time
import json
import pickle
from collections import OrderedDict



# Code borrowed from the icecube oscillation-wg fridge
# Will be replaced at a later time
from utils.plotting.mesh import add_text_values_to_mesh

#
# Standard plotting style
#
import matplotlib
matplotlib.use('agg')
from matplotlib import rcParams
FONTSIZE=22
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 15
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
plt.rcParams.update({'font.size': 22})
N_YEARS = 10
#DATA_COLOR_SCALE = ColorScale("Set1",N_YEARS)



# Set a default start date in mjd
DEFAULT_START = Time(2010., format='decimalyear')

# Constant to approximate a year in seconds
YEAR_IN_SEC = np.pi*1e7

# Set a default json config file containing livetime information 
# of our pseudo-datasets
DEFAULT_CONFIG = 'season_config.json'

# Define the binning we want to use on each variables
# we are plotting
#
TOY_MC_VARIABLES = OrderedDict()
TOY_MC_VARIABLES['time_mjd'] = np.linspace(55196.,58849.,121)
TOY_MC_VARIABLES['reco_energy'] = np.linspace(0, 100, 31)
TOY_MC_VARIABLES['vertex_x'] = np.linspace(-200, 200, 31)
TOY_MC_VARIABLES['vertex_y'] = np.linspace(-200, 200, 31)
TOY_MC_VARIABLES['vertex_z'] = np.linspace(-200, 200, 31)
TOY_MC_VARIABLES['charge'] = np.linspace(0., 3., 31)

# Define a smaller list of variables when debugging
#
DEBUG_VARIABLES = OrderedDict()
DEBUG_VARIABLES['reco_energy'] = TOY_MC_VARIABLES['reco_energy']

########################################################################
# Bespoke function to read the toy MC hdf5 file format
# TODO: generalize in a separate script

def hdf5_table_to_dict(filename, group='/', tablename='tablename', list_of_keys=None):
    '''

    '''
    import tables

    F = tables.open_file(filename)
    T = F.get_node('{0}{1}'.format(group,tablename))

    if list_of_keys is None:
        list_of_keys = T.colnames
    else:
        assert isinstance(list_of_keys,list), 'list_of_keys must be a list if provided'

    output_dict = {}

    for k in list_of_keys:
        output_dict[k] = T.col(k)

    F.close()

    return output_dict

########################################################################
# Plotting subfunctions
#

def create_uberfigure(variable=None):
    '''
    Create an Uberfigure that will contain all
    axes needed to plot all comparison plots 
    for one variable
    '''
    assert not variable is None,'Error: No variable specified for plotting'

    # LateX-proof the name of the variable
    varkey = ' '.join(variable.split('_'))

    fig = plt.figure(figsize=(30,45))
    fig.suptitle(varkey, size=40)
    gs = gridspec.GridSpec(3, 2,hspace=0.3,wspace=0.3)

    uberfigure = {'fig':fig,
                  'grid':gs}

    return uberfigure

def fraction(numerator,denominator) :

    #Case 1 : Inputs are numbers
    if isinstance(numerator,numbers.Number) and isinstance(denominator,numbers.Number) :
        return float(numerator) / float(denominator) if denominator > 0 else 0.

    #Case 2 : Inputs are numpy arrays
    elif isinstance(numerator,np.ndarray) and isinstance(denominator,np.ndarray) :
        if numerator.dtype == np.integer : numerator = numerator.astype(float)
        if denominator.dtype == np.integer : denominator = denominator.astype(float)
        return numerator / denominator

    else :
        raise Exception("Cannot calculate fraction of types %s / %s" % (type(numerator),type(denominator)) )


class ColorScale(object):
    '''
    Define a scale of colors that can be scanned through
    '''
    def __init__(self, cmap, n=100):
        self.cmap = matplotlib.cm.get_cmap(cmap)
        self.num_steps = n
        self.counter = 0

    def get(self, i):
        return self[i]  # returns a RGBA tuple

    def __getitem__(self, i):
        if i >= self.num_steps:
            raise Exception(
                "Cannot get step %i from color scale : Only %i step(s) are defined"
                % (i, self.num_steps)
            )
        return self.cmap(fraction(i, self.num_steps))

    def get_next(self):
        i = self.counter
        self.counter += 1
        return self[i]

    def get_all(self):  # TODO Can probably make this more efficient
        return [self[i] for i in range(0, self.num_steps)]

    def reset(self):
        self.counter = 0

DATA_COLOR_SCALE = ColorScale("tab10", 10)
######################################################################
# Base class 
#

class SeasonalVariationPlots():
    '''
    Make plots of the processing level variables and other observables
    '''
    def __init__(self, mHz=False, pdfname=None, variables_to_plot=None):
        '''
        initialize the class

        '''
        self.mHz = mHz
        self.data_sets = collections.OrderedDict()
        self.histogram_sets = collections.OrderedDict()
        self.histogram_cdf_sets = collections.OrderedDict()
        self.histogram_cdf_uncertainties_sets = collections.OrderedDict()
        self.histogram_uncertainties_sets = collections.OrderedDict()
        self.livetimes = collections.OrderedDict()

        
        self.rate_time_binning = np.linspace(55700,58900,103)
        self.monthly_rates=collections.OrderedDict({'binning':self.rate_time_binning,
                                                    'counts':np.zeros(self.rate_time_binning.shape[0]-1),
                                                    'livetime':np.zeros(self.rate_time_binning.shape[0]-1)})
        
        # Define the variables to plot and load
        self.variables_to_plot = variables_to_plot
        self.variables_in_data = []
        self.variables_to_load = list(variables_to_plot.keys()) #+ ['data_livetime']


        self.pdf_handle = None

        if pdfname is not None:
            self.setup_pdf(pdfname=pdfname)

        # Init unit handling
        # Unit handling
        self.rate_unit = "mHz" if self.mHz else "Hz"
        self.rate_scaling = 1.e3 if self.mHz else 1.



    def load_data(self, season_dict=None, apply_cut=None, use_pisa = False, debug=False):

        if use_pisa:
            # Load the real IceCube data
            self._load_pisa_data(season_dict, apply_cut)

        else:
            # Load the toy data
            self._load_toy_data(season_dict, apply_cut, debug=debug)


    def _load_pisa_data(self, season_dict=None, apply_cut=None):

        for season, data in season_dict.items():

            self.livetimes[season] = data['livetime']

            new_season_dict = collections.OrderedDict()

            for k,v in season_dict[season]['container'].array_data.items():
                new_season_dict[k] = v.get('host')
                if k in list(self.variables_to_plot.keys()):
                    self.variables_in_data.append(k)

            self.data_sets[season] = new_season_dict


    def _load_toy_data(self, season_dict=None, apply_cut=None, debug=False):
        '''
        load the toy data
        '''
        assert isinstance(season_dict, str), 'ERROR: toy data loader expects a str for season_dict'
        assert season_dict.endswith('.json'), 'ERROR: season_dict must be the name of a JSON file.'

        with open(season_dict) as json_file:
            seasonal_data = json.load(json_file)

            for season, metadata in seasonal_data.items():

                table_data = hdf5_table_to_dict(metadata['event_file'],
                                                group='/',
                                                tablename='toydata',
                                                list_of_keys=self.variables_to_load)

                if debug:
                    print('fields contained in this table:')
                    for k in table_data.keys():
                        print(k)
                
                if apply_cut:
                    raise Exception('No cuts implemented yet.')


                self.data_sets[season] = table_data
                self.livetimes[season] = metadata['livetime_yr']*YEAR_IN_SEC

                self.variables_in_data = []
                for k in table_data.keys():
                    if k in list(self.variables_to_plot.keys()):
                        self.variables_in_data.append(k)





    def histogram_the_data(self):
        '''
        histogram data form the containers based on the binning that was provided


        '''
        assert self.variables_to_plot is not None,'ERROR: no variables to plot'
        assert self.data_sets is not None, 'ERROR: no data provided'


        for variable,bins in self.variables_to_plot.items():

            # skip if the variable is not in data
            if variable not in self.variables_in_data:
                print('WARNING: Could not find variable {} in the data'.format(variable))
                continue

            if bins is None:
                bins=31

            self.histogram_sets[variable] = collections.OrderedDict()
            self.histogram_uncertainties_sets[variable] = collections.OrderedDict()
            self.histogram_cdf_sets[variable] = collections.OrderedDict()
            self.histogram_cdf_uncertainties_sets[variable] = collections.OrderedDict()

            for season in self.data_sets.keys():

                data = self.data_sets[season][variable]
                c,bin_edges = np.histogram(data,bins=bins)

                # Overwrite the binning of that variable
                self.variables_to_plot[variable] = bin_edges

                self.histogram_sets[variable][season] = c.astype(float)
                self.histogram_uncertainties_sets[variable][season] = np.sqrt(c.astype(float))

                # Compute the cumulative distribution function
                cumul = np.cumsum(c.astype(float),axis=0)
                cdf = np.divide(cumul,sum(c),out=np.zeros_like(cumul),where=cumul.sum()!=0)
                self.histogram_cdf_sets[variable][season] = cdf
                self.histogram_cdf_uncertainties_sets[variable][season] = np.divide(np.sqrt(cumul),cumul.sum(),out=np.zeros_like(cumul),where=cumul.sum()!=0)


                
    ################################################################################
    # Function that computes seasonal variation metrics
    #

    def compute_chi2(self):
        '''
        Compute the ratio of each year in the sample.

        To do so, divides each histogram by the year's 
        livetime, and compute the ratio of any combination 
        of years

        Computes the chi2 between two years as :

                chi2 = (observed-expected)**2/expected

        Also records the chi2 value, uncertainties on the ratio and numbers of degrees of freedom
        '''
        self.chi2_sets = collections.OrderedDict()
        self.ratio_sets = collections.OrderedDict()
        self.ratio_unc_sets = collections.OrderedDict()
        self.chi2_ndof_sets = collections.OrderedDict()

        for variable in sorted(self.variables_in_data):

            current_chi2 = collections.OrderedDict()
            current_ratio = collections.OrderedDict()
            current_unc = collections.OrderedDict()
            current_ndof = collections.OrderedDict()

            for ref in self.data_sets.keys():

                count_denom = self.histogram_sets[variable][ref]
                unc_denom   = self.histogram_uncertainties_sets[variable][ref]
                cst_denom = self.livetimes[ref]
                assert self.livetimes[ref]>0, 'ERROR: No valid livetime for season {}'.format(year)



                current_chi2[ref]=collections.OrderedDict()
                current_ratio[ref] = collections.OrderedDict()
                current_unc[ref] = collections.OrderedDict()
                current_ndof[ref] = collections.OrderedDict()

                for year in self.data_sets.keys():


                    
                    cst_num = self.livetimes[year]
                    # Break if there is no livetime
                    assert self.livetimes[year]>0, 'ERROR: No valid livetime for season {}'.format(year)

                    count_num = self.histogram_sets[variable][year]/cst_num*cst_denom
                    unc_num = self.histogram_uncertainties_sets[variable][year]/cst_num*cst_denom


                    #
                    # Ratio of counts (compare absolute count w.r.t livetime of the reference year)
                    #
                    ratio = np.divide(count_num,count_denom,out=np.ones_like(count_num),where=count_denom!=0.)

                    #
                    # Uncertainty on that ratio
                    #
                    unc_component_numerator = np.divide(unc_num**2,count_denom**2,out=np.ones_like(count_denom),where=count_denom!=0)
                    unc_component_denominator=np.divide(unc_denom**2*count_num**2.,count_denom**4.0,out=np.ones_like(count_denom),where=count_denom!=0)
                    ratio_unc = np.sqrt(unc_component_denominator+unc_component_numerator)

                    # Remove low-statistics bins from the chi2 computation
                    mask = (count_num>10)*(count_denom!=0)
                    if sum(mask)==0:
                        chi2=np.NaN
                    else:
                        chi2 = sum((count_num[mask]-count_denom[mask])**2./count_denom[mask])
 
                    current_chi2[ref][year]  = chi2
                    current_ratio[ref][year] = ratio
                    current_unc[ref][year] = ratio_unc
                    current_ndof[ref][year] = sum(mask)


            self.chi2_sets[variable] = current_chi2
            self.ratio_sets[variable] = current_ratio
            self.ratio_unc_sets[variable] = current_unc
            self.chi2_ndof_sets[variable] = current_ndof



    ##############
    # KS-tests
    #
    def compute_KS_test(self):
        '''
        Compute KS test value from CDF sets

        the return values will be tuples of (KS test value, p-value)
        '''
        import scipy.stats as scp
        from scipy.stats import ks_2samp

        self.ks_test_sets = collections.OrderedDict()
        for variable in sorted(self.variables_in_data):

            current_ks_set = collections.OrderedDict()

            for ref in self.data_sets.keys():

                current_ks_set[ref] = collections.OrderedDict()
                #cdf_ref = self.histogram_cdf_sets[variable][ref]

                reference_data_array = self.data_sets[ref][variable]

                for year in self.data_sets.keys():

                    comparison_data_array = self.data_sets[year][variable]

                    #cdf_compare  = self.histogram_cdf_sets[variable][year]
                    #current_ks_set[ref][year] = max(np.abs(cdf_compare-cdf_ref))

                    TS, pval = ks_2samp(reference_data_array, comparison_data_array)
                    current_ks_set[ref][year] = (TS, pval)

            self.ks_test_sets[variable] = current_ks_set


    ###############
    # Monthly Rates
    #
    def compute_monthly_rates(self, toy_data=True, debug=False):
        '''
        Bin event times into pre-defined bins
        '''
        if toy_data:
            self._monthly_rates_toy(debug=debug)
        else:
            self._monthly_rates_icecube(debug=debug)


    def _monthly_rates_toy(self, season_config=DEFAULT_CONFIG, debug=False):
        '''

        '''

    def _monthly_rates_icecube(self, debug=False):
        '''
        Bin IceCube event times into a predefined
        monthly-ish time intervals. Handle cases where
        a run data overlaps between two bins

        '''
        for season in self.data_sets.keys():

            times_mjd = self.data_sets[season]['I3EventHeader.start_time.mod_julian_day_double']

            xtra_counts,_ = np.histogram(times_mjd,bins=self.monthly_rates['binning'])
            self.monthly_rates['counts']+=xtra_counts.astype(float)


            # compute the fraction of a season's livetime that falls into each bin 
            from processing.runs.run_database import RunDatabase

            D = RunDatabase(data_type='data',
                            database_file=os.environ['FRIDGE_DIR']+'/processing/samples/oscNext/processing/oscNext_run_db.json',
                            pass_num=2,
                            dataset=season[-2:])

            List_of_runs = D.get_run_ids(stage='level7')

            n=0
            for run in List_of_runs:
                #print('run: ',run)
                meta = D.get_run_metadata(run_id=run)
                #print(self.monthly_rates['livetime'])
                for subrun in meta:
                    st = subrun['first_event_time'].mod_julian_day_double # event start time 
                    lt = subrun['last_event_time'].mod_julian_day_double  # event end time
                    L  = subrun['livetime']                               # livetime

                    #print('subrun start time: ',st,' \t end time: ',lt,' \t livetime: ',L)
                    bin_first = np.digitize(st, bins=self.rate_time_binning)-1
                    bin_last  = np.digitize(lt, bins=self.rate_time_binning)-1

                    #print('first_bin: ',bin_first,' \t last_bin:', bin_last)

                    # Easy case: the subrun is contained in a single bin
                    if bin_first==bin_last:
                        self.monthly_rates['livetime'][bin_first]+=L

                    # Second case: subrun overlaps a border of the histogram
                    elif bin_last-bin_first==1:

                        # top edge of the bin in which the subrun starts
                        top_edge = self.rate_time_binning[bin_first+1]
                        assert top_edge-st>0,'ERROR: top_edge-starting time should be positive.'
                        assert lt-top_edge>0,'ERROR: end time - top edge should be positive'
                        first_bin_livetime = (top_edge-st)*24*3600 # MJD must be converted into seconds
                        second_bin_livetime = (lt-top_edge)*24*3600

                        self.monthly_rates['livetime'][bin_first]+=first_bin_livetime
                        self.monthly_rates['livetime'][bin_last]+=second_bin_livetime

                        if debug:
                            print('Catching an edge case: ')
                            print('start time: ',st)
                            print('top_edge: ',top_edge)
                            print('end time: ',lt)
                            print('first livetime added: ',(top_edge-st)*24*3600)
                            print('second livetime added: ',(lt-top_edge)*24*3600 )
                    else:
                        if debug:
                            print('starting_time: ',st)
                            print('ending_time: ',lt)
                            print('livetime: ',L)
                            print('1st bin time: ',first_bin_livetime)
                            print('2nd bin time: ',second_bin_livetime)
                            print('bin edges: ',self.rate_time_binning)
                        raise Exception('Weird livetime situation')



    ################################################################################
    # Plotting Functions (copied from season_comparison)
    #

    def plot_yearly_chi2(self,hist_key=None,fig=None, grid_elements=None):
        '''
        Plot the ratios and chi2 maps for a variable
        '''
        assert not hist_key is None,'Error: No hist_key specified for plotting'

        return_fig=False
        if fig is None:
            fig = plt.figure(figsize=(8,7))
            #fig.suptitle(hist_key, size=15)
            gs = gridspec.GridSpec(1, 1,hspace=0.4,wspace=0.4)
            chi2_plot = fig.add_subplot(gs[0])
            return_fig=True

        else:
            assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'
            chi2_plot = fig.add_subplot(grid_elements)


        #
        #  Plot the chi2 map
        #
        oned_axis = np.arange(len(list(self.chi2_sets[hist_key].keys()))+1)

        reference_axis, comparison_axis = np.meshgrid(oned_axis,oned_axis)

        chi2_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])
        i=0
        j=0
        for r in list(self.chi2_sets[hist_key].keys()):

            for y in list(self.chi2_sets[hist_key][r].keys()):

                ndof = float(self.chi2_ndof_sets[hist_key][r][y])
                C2 = self.chi2_sets[hist_key][r][y]
                chi2_map[i][j] = np.divide(C2, ndof,out = np.zeros_like(C2),where=ndof!=0)

                if np.isclose(chi2_map[i][j],0., atol=1e-5):
                    chi2_map[i][j] = 0.

                j+=1
            j=0
            i+=1

        chi2_mesh_handle = chi2_plot.pcolormesh(reference_axis, comparison_axis, chi2_map.T, cmap='Greens',)
        chi2_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
        chi2_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
        chi2_plot.set_xticklabels(list(self.histogram_sets[hist_key].keys()), minor=False, rotation=45)
        chi2_plot.set_yticklabels(list(self.chi2_sets[hist_key]['2012'].keys()), minor=False)
        plot_title = hist_key.replace('_', ' ')
        chi2_plot.set_ylabel(r'\textbf{Comparison Year}', fontweight='bold')
        chi2_plot.set_xlabel(r'\textbf{Reference year}', fontweight='bold')

        chi2_plot.set_title(r'\textbf{Reduced $\chi^{2}$ - '+plot_title+'}', fontweight='bold', pad=20)
        add_text_values_to_mesh(chi2_mesh_handle, fmt="%.2g", finiteOnly=True, fontweight='bold')

        if return_fig:
            return fig, chi2_plot

    def plot_yearly_ks_test(self, hist_key=None, fig=None, grid_elements=None):
        '''
        plot the ks_tests
        '''
        import matplotlib as mpl
        assert not hist_key is None,'Error: No hist_key specified for plotting'

        return_fig=False
        if fig is None:
            fig = plt.figure(figsize=(8,7))
            #fig.suptitle(hist_key, size=15)
            gs = gridspec.GridSpec(1, 1,hspace=0.4,wspace=0.4)
            KS_plot = fig.add_subplot(gs[0])
            return_fig=True

        else:
            assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'
            KS_plot = fig.add_subplot(grid_elements)


        oned_axis = np.arange(len(list(self.ks_test_sets[hist_key].keys()))+1)
        reference_axis,comparison_axis = np.meshgrid(oned_axis,oned_axis)
        KS_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])


        i=0
        j=0
        for r in list(self.ks_test_sets[hist_key].keys()):

            for y in list(self.ks_test_sets[hist_key][r].keys()):

                TS, pval = self.ks_test_sets[hist_key][r][y]
                KS_map[i][j] = pval

                #KS_plot.text(i+0.2,j+0.4, '%.2g'%pval, fontsize=20, fontweight='bold')
                #KS_plot.text(i+0.2,j+0.1, '(TS = %.03f)'%TS)
                j+=1
            j=0
            i+=1

        ks_map_handle = KS_plot.pcolormesh(reference_axis, comparison_axis, KS_map.T, cmap='Blues', norm=mpl.colors.LogNorm(vmax=1.0, vmin=1.e-3))
        KS_plot.set_xticks(oned_axis[:-1]+0.5, minor=False)
        KS_plot.set_yticks(oned_axis[:-1]+0.5, minor=False)
        KS_plot.set_xticklabels(list(self.ks_test_sets[hist_key].keys()), minor=False, rotation=45)
        KS_plot.set_yticklabels(list(self.ks_test_sets[hist_key]['2012'].keys()), minor=False)
        KS_plot.set_ylabel(r'\textbf{Comparison Year}', fontweight='bold')
        KS_plot.set_xlabel(r'\textbf{Reference year}', fontweight='bold')
        plot_title = hist_key.replace('_', ' ')
        KS_plot.set_title(r'\textbf{KS Test - '+plot_title+'}' , fontweight='bold', pad=20)
        #KS_plot.set_title(r'KS test Value Map', fontweight='bold', pad=20)
        add_text_values_to_mesh(ks_map_handle, fmt="%.1g", finiteOnly=True, fontweight='bold')

        cbar = fig.colorbar(ks_map_handle, ax=KS_plot)
        cbar.set_label('p-value')

        if return_fig:
            return fig, KS_plot


    def plot_monthly_rates(self, rate_mode=False, rates_only=False):
        '''
        Plot the monthly variations of the rate, for each season
        '''

        if rates_only:
            F, ax = plt.subplots(figsize=(15,5))
            F.suptitle('Monthly rates and counts')
        else:
            F = Figure(4, 1, figsize=(20,30), title='Monthly-ish rates and counts')

        C = self.monthly_rates['counts'].astype(float)
        T = self.monthly_rates['livetime'].astype(float)
        x = self.monthly_rates['binning'][:-1]+0.5*(self.monthly_rates['binning'][1:]-self.monthly_rates['binning'][:-1])

        R = np.divide(C, T, out = np.zeros_like(C), where=T!=0)/1.e-3
        err = np.divide(np.sqrt(C),T,out = np.zeros_like(C),where=T!=0)*1000.

        legend_lines = collections.OrderedDict()
        for i, qty, error, label in zip(range(3), [C,T,R],[np.sqrt(C),None,err],['Event Count','Livetime','Rate (mHz)']):
            
            if rates_only and i!=2:
                continue

            if not rates_only:
                ax = F.get_ax(i=i)
            ax.errorbar(x,qty,yerr= error, drawstyle='steps-mid',color='k',fmt='o', markersize=5)
            ax.set_ylabel(label)


            
            o = ax.axvline(55694, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(0), zorder=-1)
            legend_lines['2011'] = o
            o = ax.axvline(56062, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(1), zorder=-1)
            legend_lines['2012'] = o
            o = ax.axvline(56414, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(2), zorder=-1)
            legend_lines['2013'] = o
            o = ax.axvline(56783, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(3), zorder=-1)
            legend_lines['2014'] = o
            o = ax.axvline(57160, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(4), zorder=-1)
            legend_lines['2015'] = o
            o = ax.axvline(57528, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(5), zorder=-1)
            legend_lines['2016'] = o
            o = ax.axvline(57891, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(6), zorder=-1)
            legend_lines['2017'] = o
            o = ax.axvline(58309, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(7), zorder=-1)
            legend_lines['2018'] = o
            o = ax.axvline(58682, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(8), zorder=-1)
            legend_lines['2019'] = o
            o = ax.axvline(58998, 0., 1., linewidth=2.0, color=DATA_COLOR_SCALE.get(9), zorder=-1)
            legend_lines['2020'] = o

            ax.set_xlabel('MJD date')

            year_jan_1st = [55927, 56293, 56658, 57023, 57388, 57754, 58119, 58484, 58849]

            for k,y in enumerate(sorted(year_jan_1st)):
                if k==0:
                    o = ax.axvline(y, 0.,1., linewidth=1.0, linestyle='--', color='k', zorder=-1)
                    legend_lines['1st January'] = o
                else:
                    ax.axvline(y, 0.,1., linewidth=1.0, linestyle='--', color='k', zorder=-1)

            

            if label=='Rate (mHz)':
                # constrain  axis limits and fit a line
                ax.set_ylim([0.9,1.3])
                mask = (qty>0)*(x>56062)
                chi2_obj = Chi2Regression(linear_fit, x=x[mask],y=qty[mask], sy=error[mask])
                fit = Minuit(chi2_obj,  pedantic=False, print_level=0, a=0., b=np.mean(qty[mask]))
                fit.migrad()
                slope = fit.values['a']
                slope_err = fit.errors['a']
                intercept = fit.values['b']
                intercept_err = fit.errors['b']
                chi2 = fit.fval
                ndof = qty.shape[0]-2
                o = ax.plot(x,linear_fit(x,a=slope, b=intercept), linestyle='-', color='grey', zorder=-1)
                legend_lines['Fit'] = o
                txtbox=ax.text(0.7,0.7,'Linear Fit:\n---\n slope: {0:.3g} +/- {1:.1g}\n intercept: {2:.3g} +/- {3:.1g}\n chi2 (dof = {4}): {5:.3g}'.format(slope,
                                                                                                                                 slope_err,
                                                                                                                                 intercept,
                                                                                                                                 intercept_err,
                                                                                                                                 ndof,
                                                                                                                                 chi2),
                transform=ax.transAxes,
                alpha=1.,
                backgroundcolor='w'
                )
                txtbox.set_bbox(dict(facecolor='w', alpha=1., edgecolor='k'))
                

        #
        # add a shaded region indicating the discarded 2011 data
        #
        ax.fill_between([55694, 56061], y1=0, y2=2, alpha=0.5, color='gray', hatch='/')

        #
        # Wrap seasonal plots into a yearly time frame
        #
        if not rates_only:

            #
            # Position one legend on the side of each subplot
            #
            F.get_ax(i=2).legend(legend_lines.values(),              # List of the line objects
                   labels=legend_lines.keys(),       # The labels for each line
                   loc="center right",)
            F.tight_layout()
            return F.fig

        else:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width*0.93, box.height])
            ax.legend([x for x in legend_lines.values()],              # List of the line objects
                   labels=[l for l in legend_lines.keys()],       # The labels for each line
                   bbox_to_anchor=(1.02, 1.0),
                   fontsize='small',
                   )
                   #loc="center right",)
            #F.tight_layout()

            return F


    def plot_nan_counts(self):
        '''
        Plot a chart reporting the number
        of NaN items in each variable harvested
        '''

        fig,ax = plt.subplots(figsize=(20,10),gridspec_kw={'hspace':0.5})

        labels = list(self.numbers_of_nans.keys())
        counts = list(self.numbers_of_nans.values())
        x_position = np.arange(len(labels))


        ax.bar(x_position,counts,tick_label=labels,color='b')
        ax.set_ylabel("Number of NaN's",fontweight='bold')

        for tick in ax.get_xticklabels():
            tick.set_rotation(90)

        return fig


    def plot_histograms(self, hist_key=None, fig = None, grid_elements=None, plot_cdf=True):
        '''
        Plot histograms and CDF of the variable
        '''

        if fig is None:
            
            if plot_cdf:
                fig = plt.figure(figsize=(5*2,int(5*45/35.)))
                grid_elements = gridspec.GridSpec(2, 2, hspace=0.2, wspace=0.2, height_ratios=[5,1])
            else:
                fig = plt.figure(figsize=(7,int(5*45/35.)))
                grid_elements = gridspec.GridSpec(2, 1, hspace=0.05, wspace=0.2, height_ratios=[5,2])
        else:
            assert grid_elements is not None,'ERROR: if you provide a figure, you must provide a gridspec element.'

        # get rid of the underscore plotting issue with latex
        hist_key_label= ' '.join(hist_key.split('_'))

        ax1 = fig.add_subplot(grid_elements[0])
        ax1.set_ylabel('Rates (mHz)',fontweight='bold')
        ax1.set_xlabel(hist_key_label,fontweight='bold')

        if plot_cdf:
            ax2 = fig.add_subplot(grid_elements[1])
            ax2.set_ylabel('1D Cumulative Distribution',fontweight='bold')
            ax2.set_xlabel(hist_key_label,fontweight='bold')
        else:
            ax2=None

        # Distribution histograms
        legend_lines = OrderedDict()
        for i, year in enumerate(self.histogram_sets[hist_key].keys()) :
            
            # retrieve the histogram bin edges
            B = self.variables_to_plot[hist_key]
            if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
                x_values = np.log10(10**(B[:-1])+0.5*(10**(B[1:])-10**(B[:-1])))
            else:
                x_values = B[:-1]+0.5*(B[1:]-B[:-1])
            
            # Retrieve the histogram and associated uncertainties
            h   = self.histogram_sets[hist_key][year]
            unc = self.histogram_uncertainties_sets[hist_key][year]
            L   = self.livetimes[year]
            o=ax1.errorbar(x_values,h/L*1000.,yerr=unc/L*1000.,drawstyle='steps-mid',color = DATA_COLOR_SCALE.get(i), linewidth=2.0)
            legend_lines[year] = o

            if ax2 is not None:
                # Retrieve the cdf information
                cdf     = self.histogram_cdf_sets[hist_key][year]
                cdf_unc = self.histogram_cdf_uncertainties_sets[hist_key][year]
                ax2.errorbar(x_values,cdf,yerr=cdf_unc,drawstyle='steps-mid', label=year,color = DATA_COLOR_SCALE.get(i), linewidth=2.0)




        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width*0.9, box.height])
        ax1.legend([x for x in legend_lines.values()],              # List of the line objects
               labels=[l for l in legend_lines.keys()],       # The labels for each line
               bbox_to_anchor=(1.02, 0.7),
               fontsize='small',
               )

        if ax2 is not None:
            ax2.legend()


        if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
            ax1.set_xscale('log')
            if ax2 is not None:
                ax2.set_xscale('log')



        return fig, grid_elements, ax1

    def plot_single_ratio(self, hist_key, fig=None, ref_year='2014', grid_elements=None):
        '''
        bespoke plotter to return a single ratio instead of all ratios
        '''
        import matplotlib.gridspec as gridspec
        assert not hist_key is None,'Error: No variable specified for plotting'

        if fig is None:
            fig = plt.figure(figsize=(5,3))
            fig.suptitle(hist_key, size=20)
            grid_elements = gridspec.GridSpec(1, 1,hspace=0.4,wspace=0.4)

        else:
            assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'


        hist_key_label= ' '.join(hist_key.split('_'))

        #
        # Get parameters related to the hist_key to plot
        #
        xbins = self.variables_to_plot[hist_key]
        if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
            X = np.log10(10**(xbins[:-1])+0.5*(10**(xbins[1:])-10**(xbins[:-1])))
        else:
            X = xbins[:-1]+0.5*(xbins[1:]-xbins[:-1])
        
        #
        # Set yearly colorscale
        #
        color_year_mapping = collections.OrderedDict()
        i=0
        for year in self.histogram_sets[hist_key].keys():
            color_year_mapping[year] = DATA_COLOR_SCALE.get(i)
            i+=1

        #
        # Plot Ratios
        #           
        n=0
        ref = self.ratio_sets[hist_key][ref_year]

        ax = fig.add_subplot(grid_elements[1])
        ax.set_ylabel('Ratio to {}'.format(ref_year), fontweight='bold', multialignment='center')

        for year in list(self.ratio_sets[hist_key][ref_year].keys()):
                
            if year!=ref:
                A = self.ratio_sets[hist_key][ref_year][year]
                B = self.ratio_unc_sets[hist_key][ref_year][year]
                C = color_year_mapping[year]
                l = ax.errorbar(X, A, yerr = B, marker = '.', color=C, label=year)

        ax.set_ylim([0.8,1.2])
        ax.axhline(y=1.0, xmin=0., xmax=1., color='k', zorder=-1000)
        ax.set_xlabel(hist_key_label)
        if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
            ax.set_xscale('log')

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.9, box.height])




    def plot_ratios(self,hist_key=None,fig=None,grid_elements=None,):
        '''
        Plot Ratios of rates, for each year

        hist_key: string (name of the variable we plot)

        fig: pyplot.figure object (figure handle)

        grid_elements: GridSpec object (location or list of locations within a gridspec)

        '''
        import matplotlib.gridspec as gridspec

        assert not hist_key is None,'Error: No variable specified for plotting'

        if fig is None:
            fig = plt.figure(figsize=(15,30))
            fig.suptitle(hist_key, size=20)
            gs = gridspec.GridSpec(N_YEARS, 1,hspace=0.4,wspace=0.4)

        else:
            assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'
            gs = gridspec.GridSpecFromSubplotSpec(N_YEARS, 1,subplot_spec=grid_elements,hspace=0.4,wspace=0.4)


        #
        # Get parameters related to the hist_key to plot
        #
        xbins = self.variables_to_plot[hist_key]
        if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
            X = np.log10(10**(xbins[:-1])+0.5*(10**(xbins[1:])-10**(xbins[:-1])))
        else:
            X = xbins[:-1]+0.5*(xbins[1:]-xbins[:-1])
        
        line_labels = []
        line_handles =[]

        #
        # Set yearly colorscale
        #
        color_year_mapping = collections.OrderedDict()
        i=0
        for year in self.histogram_sets[hist_key].keys():
            color_year_mapping[year] = DATA_COLOR_SCALE.get(i)
            i+=1

        #
        # Plot Ratios
        #           
        n=0
        for ref in list(self.ratio_sets[hist_key].keys()):

            ax = fig.add_subplot(gs[n,0])
            ax.set_ylabel('Ratio to \n %s'%ref, fontweight='bold', labelpad=30, multialignment='center')

            for year in list(self.ratio_sets[hist_key][ref].keys()):
                    
                if year!=ref:
                    A = self.ratio_sets[hist_key][ref][year]
                    B = self.ratio_unc_sets[hist_key][ref][year]
                    C = color_year_mapping[year]
                    l = ax.errorbar(X, A, yerr = B, marker = '.', color=C, label=year)
                    if year not in line_labels:
                        line_handles.append(l)
                        line_labels.append(year)

            ax.set_ylim([0.8,1.2])
            ax.axhline(y=1.0, xmin=0., xmax=1., color='k', zorder=-1000)

            if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
                ax.set_xscale('log')

            n+=1

        #fig.legend(line_handles,line_labels,loc = (.5,0.5))




    def plot_nan_counts(self):
        '''
        Plot a chart reporting the number
        of NaN items in each variable harvested
        '''

        fig,ax = plt.subplots(figsize=(20,10),gridspec_kw={'hspace':0.5})

        labels = list(self.numbers_of_nans.keys())
        counts = list(self.numbers_of_nans.values())
        x_position = np.arange(len(labels))


        ax.bar(x_position,counts,tick_label=labels,color='b')
        ax.set_ylabel("Number of NaN's",fontweight='bold')

        for tick in ax.get_xticklabels():
            tick.set_rotation(90)

        return fig






    def plot_passing_rates(self):
        '''
        print out the number of events in each year
        and the number of events passing the particular
        processing level boolean.

        plot the corresponding chi2 map for both numbers
        '''


        fig = plt.figure(figsize=(30,20))
        fig.suptitle('Passing Fractions', size=20)
        gs = gridspec.GridSpec(1, 2,hspace=0.3,wspace=0.4)

        uncut_evts_plot = fig.add_subplot(gs[:,0])
        passing_evts_plot = fig.add_subplot(gs[:,1])

        oned_axis = np.arange(len(list(self.metadata.keys()))+1)
        reference_axis,comparison_axis = np.meshgrid(oned_axis,oned_axis)

        pass_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])
        uncut_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])


        i=0
        j=0
        pass_frac = collections.OrderedDict()

        for ref_year in list(self.metadata.keys()):

            ref_events_passing = float(self.metadata[ref_year]['passing_events'])
            ref_events_uncut   = self.metadata[ref_year]['uncut_events']
            pass_frac[ref_year] = ref_events_passing/ref_events_uncut


            for y in list(self.metadata.keys()):

                comp_events_passing = float(self.metadata[y]['passing_events'])
                comp_events_uncut   = self.metadata[y]['uncut_events']

                pass_map[i][j] = np.divide((comp_events_passing-ref_events_passing)**2.0,comp_events_passing,out = np.zeros(1),where=comp_events_passing!=0)
                uncut_map[i][j] = np.divide((comp_events_uncut-ref_events_uncut)**2.0,comp_events_uncut,out = np.zeros(1),where=comp_events_uncut!=0)


                passing_evts_plot.text(i+0.5,j+0.5,'%.02f'%(pass_map[i][j]))
                uncut_evts_plot.text(i+0.5,j+0.5,'%.02f'%(uncut_map[i][j]))

                j+=1
            j=0
            i+=1

        passing_evts_plot.pcolormesh(reference_axis,comparison_axis,pass_map.T,cmap='BuPu')
        passing_evts_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
        passing_evts_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
        passing_evts_plot.set_xticklabels(list(self.metadata.keys()),minor=False)
        passing_evts_plot.set_yticklabels(list(self.metadata.keys()),minor=False)

        passing_evts_plot.set_ylabel('Comparison Year')
        passing_evts_plot.set_xlabel('Reference year')
        passing_evts_plot.set_title(r'# Events passing level cut')

        uncut_evts_plot.pcolormesh(reference_axis,comparison_axis,uncut_map.T,cmap='BuPu')
        uncut_evts_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
        uncut_evts_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
        uncut_evts_plot.set_xticklabels(list(self.metadata.keys()),minor=False)
        uncut_evts_plot.set_yticklabels(list(self.metadata.keys()),minor=False)

        uncut_evts_plot.set_ylabel('Comparison Year')
        uncut_evts_plot.set_xlabel('Reference year')
        uncut_evts_plot.set_title(r'# Events in files')

        #
        # Just print the numbers
        #
        for y in list(self.metadata.keys()):
            print((y,' passing fraction: ',pass_frac[y]))



        return fig



    def plot_passing_rates(self):
        '''
        print out the number of events in each year
        and the number of events passing the particular
        processing level boolean.

        plot the corresponding chi2 map for both numbers
        '''


        fig = plt.figure(figsize=(30,20))
        fig.suptitle('Passing Fractions', size=20)
        gs = gridspec.GridSpec(1, 2,hspace=0.3,wspace=0.4)

        uncut_evts_plot = fig.add_subplot(gs[:,0])
        passing_evts_plot = fig.add_subplot(gs[:,1])

        oned_axis = np.arange(len(list(self.metadata.keys()))+1)
        reference_axis,comparison_axis = np.meshgrid(oned_axis,oned_axis)

        pass_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])
        uncut_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])


        i=0
        j=0
        pass_frac = collections.OrderedDict()

        for ref_year in list(self.metadata.keys()):

            ref_events_passing = float(self.metadata[ref_year]['passing_events'])
            ref_events_uncut   = self.metadata[ref_year]['uncut_events']
            pass_frac[ref_year] = ref_events_passing/ref_events_uncut


            for y in list(self.metadata.keys()):

                comp_events_passing = float(self.metadata[y]['passing_events'])
                comp_events_uncut   = self.metadata[y]['uncut_events']

                pass_map[i][j] = np.divide((comp_events_passing-ref_events_passing)**2.0,comp_events_passing,out = np.zeros(1),where=comp_events_passing!=0)
                uncut_map[i][j] = np.divide((comp_events_uncut-ref_events_uncut)**2.0,comp_events_uncut,out = np.zeros(1),where=comp_events_uncut!=0)


                passing_evts_plot.text(i+0.5,j+0.5,'%.02f'%(pass_map[i][j]))
                uncut_evts_plot.text(i+0.5,j+0.5,'%.02f'%(uncut_map[i][j]))

                j+=1
            j=0
            i+=1

        passing_evts_plot.pcolormesh(reference_axis,comparison_axis,pass_map.T,cmap='BuPu')
        passing_evts_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
        passing_evts_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
        passing_evts_plot.set_xticklabels(list(self.metadata.keys()),minor=False)
        passing_evts_plot.set_yticklabels(list(self.metadata.keys()),minor=False)

        passing_evts_plot.set_ylabel('Comparison Year')
        passing_evts_plot.set_xlabel('Reference year')
        passing_evts_plot.set_title(r'# Events passing level cut')

        uncut_evts_plot.pcolormesh(reference_axis,comparison_axis,uncut_map.T,cmap='BuPu')
        uncut_evts_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
        uncut_evts_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
        uncut_evts_plot.set_xticklabels(list(self.metadata.keys()),minor=False)
        uncut_evts_plot.set_yticklabels(list(self.metadata.keys()),minor=False)

        uncut_evts_plot.set_ylabel('Comparison Year')
        uncut_evts_plot.set_xlabel('Reference year')
        uncut_evts_plot.set_title(r'# Events in files')

        #
        # Just print the numbers
        #
        for y in list(self.metadata.keys()):
            print((y,' passing fraction: ',pass_frac[y]))


        return fig


    ########################################
    ########################################
    # Main plotting function

    def plot_func(self) :

        #
        # Passing rates
        #
        #passing_rate_fig = self.plot_passing_rates()
        #self.pdf_handle.savefig(passing_rate_fig)
        #plt.close(passing_rate_fig)

        #
        # Monthly rates
        #

        #month_fig = self.plot_monthly_rates(rate_mode=True)
        #self.pdf_handle.savefig(month_fig)
        #plt.close(month_fig)

        

        #
        # NaN counts
        #
        #fig = self.plot_nan_counts()
        #self.pdf_handle.savefig(fig)
        #plt.close(fig)
        
        #
        # Plot the histograms of individual variables
        #
        yscale = "linear"
        if self.mHz: rate_unit='(mHz)'
        else: rate_unit = '(Hz)'

        for variable_name in self.histogram_sets.keys():
            print(('plotting info about ',variable_name,'...'))

            # Create a figure handle 
            uberfig = create_uberfigure(variable=variable_name)

            self.plot_histograms(hist_key=variable_name, fig = uberfig['fig'], grid_elements = uberfig['grid'])
            self.plot_ratios( hist_key=variable_name, fig = uberfig['fig'],grid_elements = uberfig['grid'][1:,0])
            self.plot_yearly_chi2(hist_key=variable_name, fig = uberfig['fig'], grid_elements = uberfig['grid'][1,1])
            self.plot_yearly_ks_test(hist_key=variable_name, fig = uberfig['fig'], grid_elements = uberfig['grid'][2,1])

            self.pdf_handle.savefig(uberfig['fig'])


            plt.close(uberfig['fig'])
            del uberfig
        



    ################################################################################


    def setup_pdf(self,pdfname=None):

        self.pdf_handle = PdfPages(pdfname)

    def close_pdf(self):

        self.pdf_handle.close()
   



#######################################################################################
# Function to Harvest the desired keys

#######################################################################################
# Main executor
if __name__ == "__main__" :


    import argparse
    parser = argparse.ArgumentParser('Compare several seasons of data to each other.')

    # Debug-mode
    parser.add_argument('--debug', help='debug mode.', action='store_true')

    # input metadata config
    parser.add_argument('-m', '--metadata', help='path to a metadata file storing seasonal info.', default=DEFAULT_CONFIG)

    # output name
    parser.add_argument('-o', '--output-file', help='output pdf file name', default='test_season_comparison.pdf')
    parser.add_argument('-c', '--apply-cuts', help='Apply the cuts associated to processing level', action='store_true')

    # Select the type of comparison you want
    parser.add_argument('-yr', '--years', help='list the season you want to compare to', nargs='+', type=int, default= [12,13,14,17,18])
    args = parser.parse_args()



    # Steer harvesting
    if args.debug:
        years = [12,13]
        variables_to_harvest = DEBUG_VARIABLES
    else:
        variables_to_harvest = OrderedDict()
        for k,v in TOY_MC_VARIABLES.items():
            variables_to_harvest[k] = v

        years = args.years




    comparator = SeasonalVariationPlots(mHz=False, 
                                        pdfname=args.output_file,
                                        variables_to_plot=variables_to_harvest)

    comparator.load_data(season_dict=args.metadata, 
                         apply_cut=args.apply_cuts,
                         use_pisa=False,
                         debug=args.debug)
    
    # plot monthly rates
    comparator.compute_monthly_rates()

    # bin and compute quantities
    comparator.histogram_the_data()
    comparator.compute_chi2()
    comparator.compute_KS_test()

    # Plotting
    comparator.plot_func()
    comparator.close_pdf()

    print('pdf saved: ', args.output_file)

        
        
    #
    # Save the rate plot (TODO)
    #
    #fig = comparator.plot_monthly_rates(rates_only=True)
    #plt.savefig('thesis_plot_season_comparison_rates.png')


    # save a few selected plots as individual png figures
    for v in ['reco_energy', 'vertex_x']:

        print('plotting: ',v)
        fig, ax = comparator.plot_yearly_ks_test(hist_key=v)
        fig.set_size_inches([10.,9.5])
        plt.savefig('season_comparison_{}_ks.png'.format(v))

        fig, ax = comparator.plot_yearly_chi2(hist_key=v)
        plt.savefig('season_comparison_{}_chi2.png'.format(v))     
        fig, gs, ax = comparator.plot_histograms(hist_key=v, plot_cdf=False)
        ax.get_xaxis().set_ticks([])
        comparator.plot_single_ratio(hist_key=v, fig=fig, grid_elements=gs)
        s = fig.get_size_inches()
        fig.set_size_inches([10.,7.])

        plt.savefig('season_comparison_{}_hist.png'.format(v))



        print('saved some figures to PNG')