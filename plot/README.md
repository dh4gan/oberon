This collection of python scripts plots data from `dh4gan/oberon`

Data from this C++ code comes in several forms:

***N Body Output Data***

This is plotted using

`plot_orbit.py`
`plot_separations.py`
`plot_positions.py`

and animated (crudely) using 

`show_orbit.py`

***Latitudinal Energy Balance Model (EBM) data***

This is plotted using

`plot_EBMlog.py` - which plots the global properties stored in the `.log` file

`plot_EBMlog_multiplefiles.py` - which plots the global properties (`.log`) of multiple objects with EBM data
`plot_EBMlog_multiplefiles_anomaly.py` - which plots the deviation from the mean of global properties (`.log`) of multiple objects with EBM data

`plot_timeaverage_EBM.py` - which collects several latitudinal snapshots and plots time averaged curves

`construct_lat_vs_t.py` - which plots the latitudinal properties as a function of time on a single plot

`plot_EBMlat.py` - plots the temperature as a function of latitude and time, (`.lat`)

`EBM_movie.py` - create several latitudinal plots for the sake of animation

`plot_EBMlog_periodogram.py` - calculates periodograms of the `.log` data using Fast Fourier Transforms


The `filefinder` package (enclosed) is used to identify files with the correct format for plotting.  
In the case of `plot_multiple_EBMlog.py`, it is used to find all `.log` files in the directory

The `io_EBM` package handles I/O for the EBM data
The `io_nbody` package handles I/O for the N Body files
The `io_paramfile` package handles I/O for input parameter files (in development)
