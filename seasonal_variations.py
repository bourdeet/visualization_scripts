'''
Script to perform a comparison of seasonal
data.


Etienne Bourbeau
'''

# Standard Python stuff
import os, sys, collections, glob, datetime
import numpy as np
import pickle
from collections import OrderedDict

# Standard plotting stuff
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams.update({'font.size': 22})
N_YEARS = 9
DATA_COLOR_SCALE = ColorScale("Set1",N_YEARS)


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

	fig = plt.figure(figsize=(30,45))
	fig.suptitle(variable, size=40)
	gs = GridSpec(3, 2,hspace=0.3,wspace=0.3)

	uberfigure = {'fig':fig,
				  'grid':gs}

	return uberfigure



######################################################################
# Base class 
#

class SeasonalVariationPlots():
	'''
	Make plots of the processing level variables and other observables
	'''
	def __init__(self,mHz=False,pdfname=None,variables_to_plot = None):
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
		self.variables_to_load = list(variables_to_plot.keys()) + ['data_livetime']


		self.pdf_handle = None

		if pdfname is not None:
			self.setup_pdf(pdfname=pdfname)

		# Init unit handling
		# Unit handling
		self.rate_unit = "mHz" if self.mHz else "Hz"
		self.rate_scaling = 1.e3 if self.mHz else 1.



	def load_data(self,season_dict=None,apply_cut=None, use_pisa = False):

		if use_pisa:
	
			for season, data in season_dict.items():

				self.livetimes[season] = data['livetime']

				new_season_dict = collections.OrderedDict()

				for k,v in season_dict[season]['container'].array_data.items():
					new_season_dict[k] = v.get('host')
					if k in list(self.variables_to_plot.keys()):
						self.variables_in_data.append(k)

				self.data_sets[season] = new_season_dict

		else:

			for season, filepath in season_dict.items():
				metadata,new_dict = read_hdf5(filepath,choose=self.variables_to_load)

				for k in new_dict.keys():
					print(k)
				# if a mask is selected:
				if apply_cut:
					M = new_dict['L7_MuonClassifier_ProbNu']>0.7
					new_new_dict = collections.OrderedDict()
					for k,v in new_dict.items():
						new_new_dict[k] = v[M]

					new_dict = new_new_dict

				self.data_sets[season] = new_dict
				self.livetimes[season] = sum(metadata['livetimes'])

				self.variables_in_data = []
				for k in new_dict.keys():
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
	def compute_monthly_rates(self):
		'''
		compute the monthly rates
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
						print('Catching an edge case: ')
						print('start time: ',st)
						print('top_edge: ',top_edge)
						print('end time: ',lt)
						print('first livetime added: ',(top_edge-st)*24*3600)
						print('second livetime added: ',(lt-top_edge)*24*3600 )

					else:
						print('starting_time: ',st)
						print('ending_time: ',lt)
						print('livetime: ',L)
						print('1st bin time: ',first_bin_livetime)
						print('2nd bin time: ',second_bin_livetime)
						print('bin edges: ',self.rate_time_binning)
						raise Exception('Weird livetime situation')

				#print(self.monthly_rates['livetime'])
				#print(self.rate_time_binning)


	################################################################################
	# PLotting Functions (copied from season_comparison)
	#
	def plot_ratios(self,hist_key=None,fig = None,grid_elements=None):
		'''
		Plot Ratios of rates, for each year
		'''

		assert not hist_key is None,'Error: No variable specified for plotting'

		if fig is None:
			fig = plt.figure(figsize=(15,30))
			fig.suptitle(hist_key, size=20)
			gs = GridSpec(N_YEARS, 1,hspace=0.4,wspace=0.4)

		else:
			assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'
			gs = gridspec.GridSpecFromSubplotSpec(N_YEARS, 1,subplot_spec=grid_elements,hspace=0.4,wspace=0.4)


		#
		# Get parameters related to the hist_key to plot
		#
		xbins = self.variables_to_plot[hist_key]
		X = xbins[:-1]+0.5*(xbins[1:]-xbins[:-1])
		line_labels = []
		line_handles =[]

		#
		# Set yearly colorscale
		#
		color_year_mapping = collections.OrderedDict()
		i=0
		for year in list(self.histogram_sets[hist_key].keys()):
			color_year_mapping[year] = DATA_COLOR_SCALE.get(i)
			i+=1

		#
		# Plot Ratios
		#			
		n=0
		for ref in list(self.ratio_sets[hist_key].keys()):

			ax = fig.add_subplot(gs[n,0])
			#ax.set_xlabel(hist_key)
			ax.set_ylabel('ratio to %s'%ref,fontweight='bold')

			for year in list(self.ratio_sets[hist_key][ref].keys()):
					
				if year!=ref:
					l = ax.errorbar(X,
									self.ratio_sets[hist_key][ref][year],
									yerr = self.ratio_unc_sets[hist_key][ref][year], 
									marker = '.',
									color=color_year_mapping[year],
									label=year)

					if year not in line_labels:
						line_handles.append(l)
						line_labels.append(year)

			ax.set_ylim([0.2,1.8])

			n+=1

		#fig.legend(line_handles,line_labels,loc = (.5,0.5))




	def plot_yearly_chi2(self,hist_key=None,fig=None, grid_elements=None):
		'''
		Plot the ratios and chi2 maps for a variable
		'''
		assert not hist_key is None,'Error: No hist_key specified for plotting'


		if fig is None:
			fig = plt.figure(figsize=(10,10))
			fig.suptitle(hist_key, size=15)
			gs = GridSpec(1, 1,hspace=0.4,wspace=0.4)
			chi2_plot = fig.add_subplot(gs[0])

		else:
			assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'
			chi2_plot = fig.add_subplot(grid_elements)


		#
		#  Plot the chi2 map
		#
		oned_axis = np.arange(len(list(self.chi2_sets[hist_key].keys()))+1)

		reference_axis,comparison_axis = np.meshgrid(oned_axis,oned_axis)

		chi2_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])
		i=0
		j=0
		for r in list(self.chi2_sets[hist_key].keys()):

			for y in list(self.chi2_sets[hist_key][r].keys()):

				ndof = float(self.chi2_ndof_sets[hist_key][r][y])
				C2 = self.chi2_sets[hist_key][r][y]
				chi2_map[i][j] = np.divide(C2,ndof,out = np.zeros_like(C2),where=ndof!=0)

				chi2_plot.text(i+0.3,j+0.4,'%.02f'%(chi2_map[i][j]))

				j+=1
			j=0
			i+=1

		chi2_plot.pcolormesh(reference_axis,comparison_axis,chi2_map.T,cmap='Greens',)
		chi2_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
		chi2_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
		chi2_plot.set_xticklabels(list(self.histogram_sets[hist_key].keys()),minor=False)
		chi2_plot.set_yticklabels(list(self.chi2_sets[hist_key]['2012'].keys()),minor=False)

		chi2_plot.set_ylabel('Comparison Year',fontweight='bold')
		chi2_plot.set_xlabel('Reference year',fontweight='bold')
		chi2_plot.set_title(r'Reduced $\chi^{2}$ Map ',fontweight='bold')



	def plot_yearly_ks_test(self, hist_key=None, fig=None, grid_elements=None):
		'''
		plot the ks_tests
		'''
		assert not hist_key is None,'Error: No hist_key specified for plotting'


		if fig is None:
			fig = plt.figure(figsize=(10,10))
			fig.suptitle(hist_key, size=15)
			gs = GridSpec(1, 1,hspace=0.4,wspace=0.4)
			KS_plot = fig.add_subplot(gs[0])

		else:
			assert not grid_elements is None,'Must specify the gridspec elements where to create the axes'
			KS_plot = fig.add_subplot(grid_elements)


		oned_axis = np.arange(len(list(self.chi2_sets[hist_key].keys()))+1)
		reference_axis,comparison_axis = np.meshgrid(oned_axis,oned_axis)
		KS_map = np.zeros([oned_axis.shape[0]-1,oned_axis.shape[0]-1])


		i=0
		j=0
		for r in list(self.ks_test_sets[hist_key].keys()):

			for y in list(self.ks_test_sets[hist_key][r].keys()):

				TS, pval = self.ks_test_sets[hist_key][r][y]
				KS_map[i][j] = pval
				KS_plot.text(i+0.3,j+0.4, '%.02f'%pval)
				KS_plot.text(i+0.2,j+0.1, '(TS = %02f)'%TS)
				j+=1
			j=0
			i+=1

		KS_plot.pcolormesh(reference_axis,comparison_axis,KS_map.T,cmap='Blues')#,vmax=0.05)
		KS_plot.set_xticks(oned_axis[:-1]+0.5,minor=False)
		KS_plot.set_yticks(oned_axis[:-1]+0.5,minor=False)
		KS_plot.set_xticklabels(list(self.ks_test_sets[hist_key].keys()),minor=False)
		KS_plot.set_yticklabels(list(self.ks_test_sets[hist_key]['2012'].keys()),minor=False)
		KS_plot.set_ylabel('Comparison Year',fontweight='bold')
		KS_plot.set_xlabel('Reference year',fontweight='bold')
		KS_plot.set_title(r'KS test Value Map',fontweight='bold')



	def plot_monthly_rates(self,rate_mode = False):
		'''
		Plot the monthly variations of the rate, for each season
		'''

		F = Figure(3,1,figsize=(20,30),title='Monthly-ish rates and counts')

		C = self.monthly_rates['counts'].astype(float)
		T = self.monthly_rates['livetime'].astype(float)
		x = self.monthly_rates['binning'][:-1]+0.5*(self.monthly_rates['binning'][1:]-self.monthly_rates['binning'][:-1])


		R = np.divide(C,T,out = np.zeros_like(C),where=T!=0)/1.e-3
		err = np.divide(np.sqrt(C),T,out = np.zeros_like(C),where=T!=0)*1000.

		for i, qty, error, label in zip(range(3), [C,T,R],[np.sqrt(C),None,err],['Event Count','Livetime','Rate (mHz)']):
			
			ax = F.get_ax(i=i)
			ax.errorbar(x,qty,yerr=error,drawstyle='steps-mid',color='k',fmt='o')
			ax.set_ylabel(label)
			if i==2:
				ax.set_ylim([2.5,3.5])

			legend_lines = collections.OrderedDict()
			o = ax.plot([55694,55694],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(0),label='2011')
			legend_lines['2011'] = o[0]
			o = ax.plot([56062,56062],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(1),label='2012')
			legend_lines['2012'] = o[0]
			o = ax.plot([56414,56414],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(2),label='2013')
			legend_lines['2013'] = o[0]
			o = ax.plot([56783,56783],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(3),label='2014')
			legend_lines['2014'] = o[0]
			o = ax.plot([57160,57160],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(4),label='2015')
			legend_lines['2015'] = o[0]
			o = ax.plot([57528,57528],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(5),label='2016')
			legend_lines['2016'] = o[0]
			o = ax.plot([57891,57891],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(6),label='2017')
			legend_lines['2017'] = o[0]
			o = ax.plot([58309,58309],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(7),label='2018')
			legend_lines['2018'] = o[0]
			o = ax.plot([58682,58682],[0.,1.15*max(qty)],linewidth=2.0,color=DATA_COLOR_SCALE.get(8),label='2019')
			legend_lines['2019'] = o[0]
		
			ax.set_xlabel('MJD date')

			year_jan_1st = [55927,56293,56658,57023,57388,57754,58119,58484]

			for y in sorted(year_jan_1st):
				ax.plot([y,y],[0.,1.15*max(R)],linewidth=1.0,linestyle='--',color='k')


		F.fig.legend(legend_lines.values(),              # List of the line objects
				   labels=legend_lines.keys(),       # The labels for each line
				   loc="center right",) 
	
		return F.fig


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


	def plot_histograms(self, hist_key=None, fig = None, grid_elements=None):
		'''
		Plot histograms and CDF of the variable
		'''


		ax1 = fig.add_subplot(grid_elements[0])
		ax1.set_ylabel('Rates (mHz)',fontweight='bold')
		ax1.set_xlabel(hist_key,fontweight='bold')

		ax2 = fig.add_subplot(grid_elements[1])
		ax2.set_ylabel('1D Cumulative Distribution',fontweight='bold')
		ax2.set_xlabel(hist_key,fontweight='bold')
		ratio_grid = grid_elements[3]

		# Distribution histograms
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
			ax1.errorbar(x_values,h/L*1000.,yerr=unc/L*1000.,drawstyle='steps-mid',label=year, color = DATA_COLOR_SCALE.get(i), linewidth=2.0)


			# Retrieve the cdf information
			cdf     = self.histogram_cdf_sets[hist_key][year]
			cdf_unc = self.histogram_cdf_uncertainties_sets[hist_key][year]
			ax2.errorbar(x_values,cdf,yerr=cdf_unc,drawstyle='steps-mid', label=year, color = DATA_COLOR_SCALE.get(i), linewidth=2.0)


		ax1.legend()
		ax2.legend()

		if 'total_energy' in hist_key or 'cascade_energy' in hist_key:
			ax1.set_xscale('log')
			ax2.set_xscale('log')






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

		# month_fig = self.plot_monthly_rates(rate_mode=True)
		# self.pdf_handle.savefig(month_fig)
		# plt.close(month_fig)

		

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

		





	def plot_passing_rates(self):
		'''
		print out the number of events in each year
		and the number of events passing the particular
		processing level boolean.

		plot the corresponding chi2 map for both numbers
		'''


		fig = plt.figure(figsize=(30,20))
		fig.suptitle('Passing Fractions', size=20)
		gs = GridSpec(1, 2,hspace=0.3,wspace=0.4)

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

	################################################################################


	def setup_pdf(self,pdfname=None):

		self.pdf_handle = PdfPages(pdfname)

	def close_pdf(self):

		self.pdf_handle.close()





#########################################################################################################
#
# Functions that plot things
#



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
			gs = GridSpec(N_YEARS, 1,hspace=0.4,wspace=0.4)

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
			ax.set_ylabel('ratio to %s'%ref,fontweight='bold')

			for year in list(self.ratio_sets[hist_key][ref].keys()):
					
				if year!=ref:
					A = self.ratio_sets[hist_key][ref][year]
					B = self.ratio_unc_sets[hist_key][ref][year]
					C = color_year_mapping[year]
					l = ax.errorbar(X,A,yerr = B,marker = '.',color=C,label=year)
					if year not in line_labels:
						line_handles.append(l)
						line_labels.append(year)

			ax.set_ylim([0.2,1.8])

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
		gs = GridSpec(1, 2,hspace=0.3,wspace=0.4)

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



#######################################################################################
# Function to Harvest the desired keys

#######################################################################################
# Main executor
if __name__ == "__main__" :


	import argparse
	parser = argparse.ArgumentParser('Compare several seasons of data to each other.')

	# Debug-mode
	parser.add_argument('--debug', help='debug mode.', action='store_true')

	# output name
	parser.add_argument('-os', '--output-stem', help='output file name. only supply the root, not the extension',default='test_season_comparison')
	parser.add_argument('-is', '--input-stem', help='input file stem name of the files you want to plot',required=True)

	parser.add_argument('-c', '--apply-cuts', help='Apply the cuts associated to processing level', action='store_true')
	parser.add_argument('-f', '--fixed-period', help='harvest only runs between August and december of dataset year',action='store_true')

	# Select the type of comparison you want
	parser.add_argument('-yr', '--years', help='list the season you want to compare to', nargs='+', type=int, default= [12,13,14,17,18])
	parser.add_argument('-ss', '--subsample', help='dataset to compare: test,burn or full', default='test')

	args = parser.parse_args()


	# Steer harvesting
	if args.debug:
		years = [12,13]
		variables_to_harvest = DEBUG_VARIABLES
	else:

		variables_to_harvest = OrderedDict()
		for k,v in ALL_OSCNEXT_VARIABLES.items():
			variables_to_harvest[k] = v

		for k,v in SPECIAL.items():
			variables_to_harvest[k] = v
		years = args.years



[	pdfname = args.output_stem+'.pdf'
	variables_to_harvest['I3EventHeader.start_time.mod_julian_day_double'] = np.linspace(0.,60000,31)

	comparator = SeasonalVariationPlots(mHz=False,pdfname=pdfname, variables_to_plot = variables_to_harvest)

	# format the input file names into a dict
	input_files = collections.OrderedDict()

	for year in years:
		filename = args.input_stem+'_data_20%0i.hdf5'%year
		if os.path.isfile(filename):
			input_files['20%0i'%year] = filename
		else:
			print('could not find a file for year ',year)



	comparator.load_data(season_dict=input_files,apply_cut=args.apply_cuts)
	#comparator.compute_monthly_rates()

	comparator.histogram_the_data()
	comparator.compute_chi2()
	comparator.compute_KS_test()

	# Plotting
	comparator.plot_func()
	comparator.close_pdf()

	print('pdf saved: ',pdfname)]


