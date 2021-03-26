# visualization_scripts
A collection of scripts created over the years to read and visualize data.


# Pseudo Data generated for this plot

The pseudo data used for the plots features here were produced using the following set of commands:

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

the plots themselves can be obtained by running:

```
python seasonal_variations.py
```