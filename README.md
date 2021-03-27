# visualization_scripts
A collection of scripts created over the years to read and visualize data. 


Solarized dark             |  Solarized Ocean
:-------------------------:|:-------------------------:
![Energy distribution](/example_plots/season_comparison_reco_energy_hist.png) | ![KS Test](/example_plots/season_comparison_reco_energy_ks.png)

##  A. Data stability plots

These plots were designed to compare different data samples from IceCube taken during different detector seasons, and assess whether these samples are compatible with one another. This is done by looking at the 1D distribution of several control variables in the data.

The seasonal_variation.py script is a Python script that generates these plots and prints them out into a single PDF that summarizes the check on a variable in a single page. The plots features in a typical page are:

* A histogram of the 1D distribution
* A histogram of the CDF of the distribution
* A set of ratio plots, where we compare the 1D distributions of every season to a reference season
* A chi2 map, quantifying the overall level of agreement between the 1D distributions of any pair of season
* A KS test, quantifying the agreement between the shapes of any pair of distributions.

[This PDF](example_plots/seasonal_variation_guide.pdf) provides additional explanation regarding how to interpret the figures produced by this script. The same folder also contains an example PDF file containing plots generated using a set of 10 seasons of pseudo-data drawn from identical distributions.

### About the toy event generator

Since the original plots were generated using proprietary IceCube data, I have created a small script that generates pseudo-data to mimick the sort of information that is normally used in IceCube. The script can produce any specified livetime, using any user-specified baseline event rate. It will draw the following set of variables:

* An mjd trigger time from a uniform distribution
* A reconstructed energy value from an asymmetric chi2 distribution
* A set of vertex x/y/z coordinates using a uniform theta/phi distribution, but a normally distributed radius
* A charge value drawn from an SPE template based on the one featured in arxiv:2002.00997

When generating the data, the user can specify different livetime for different season, or choose to generate event charge according to a different value of SPE peak. The goal here is then to see how these changes between season can be detected in the subsequent plots.

### Pseudo Data generated for this plot

For a vanilla set of plots (where all seasons are identical), the pseudo data used for the plots features here were produced using the following set of commands:

```
python toy_mc_generator.py -sd 55197 -o pseudo_data_2010.hdf5
python toy_mc_generator.py -sd 55562 -o pseudo_data_2011.hdf5
python toy_mc_generator.py -sd 55927 -o pseudo_data_2012.hdf5
python toy_mc_generator.py -sd 56293 -o pseudo_data_2013.hdf5
python toy_mc_generator.py -sd 56658 -o pseudo_data_2014.hdf5
python toy_mc_generator.py -sd 57023 -o pseudo_data_2015.hdf5
python toy_mc_generator.py -sd 57388 -o pseudo_data_2016.hdf5
python toy_mc_generator.py -sd 57754 -o pseudo_data_2017.hdf5
python toy_mc_generator.py -sd 58119 -o pseudo_data_2018.hdf5
python toy_mc_generator.py -sd 58484 -o pseudo_data_2019.hdf5
```
Where we have generated 10 seasons of data, each with a livetime of exactly 1 year and an identical SPE distribution. The plots themselves can be obtained by running:

```
python seasonal_variations.py
```

In the pdf called seasonal_comparison_mismatch.pdf, I have created the same set of plots, but for a case in which a) the 2011 season has 5% smaller livetime than the other season and b) the 2015 season has a slightly different SPE peak (1.02 instead of 1.029).

```
python toy_mc_generator.py -sd 55562 -o pseudo_data_2011.hdf5 -l 0.95
python toy_mc_generator.py -sd 57023 -o pseudo_data_2015.hdf5 -spe 1.02
```

From that PDF it can be seen that the chi2 mapping is really good at finding changes in overall rates, whereas the KS test is king at spotting small changes in the shape of a distribution.

